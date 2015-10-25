#include "Dstar.h"
#include "SurfaceExplorer.h"

/* SurfaceExplorer()
 * ----------------------------------------------------------------------------
 * Set the distance per cell in the resulotion to start using the explorer.
 */
SurfaceExplorer::SurfaceExplorer(double r, double d2s, double zs, double zm) {
	resolution = r;
	distance2Surface = d2s;
	zStep = zs;
	zMax = zm;
}

/* init(Vector3d s, Vector3d o, vector<MatrixXd>& m)
 * ----------------------------------------------------------------------------
 * Initializing the surface explorer allows you to set the start point of the
 * robot and the sensor location (observer).  It also allows you to set a base
 * map, or to the best of your knowledge what the base map looks like. 
 * The map is checekd for sizing (i.e. each layer much have the same
 * dimensions).  
 */
bool SurfaceExplorer::init(Vector3d s, Vector3d o, vector<MatrixXd>& m) {
	start = s;
	observer = o;

	// Make sure all the layers of the map are the same size
	int mapXYsize = m[0].size();
	for(vector<int>::size_type i=0; i!=m.size(); i++) {
		if (m[i].size() != mapXYsize) {
			throw mapException;
		}
	}

	// Store the map
	map3D = m;

	// zMax handles planning limits, if the map isn't big enough, need 
	// to adjust
	if(map3D.size()<zMax) zMax = map3D.size();

	return true;
}

/* checkNeighbors(Vector3d location, int connectivity, int value)
 * ----------------------------------------------------------------------------
 * Checks the neigboring cells for a certain value.  Handles 8-connectivity
 * and 4 connectivity, 4 is default. Returns a truth matrix (3x3) where the
 * center is always zero.  
 */
Matrix3d SurfaceExplorer::checkNeighbors(Vector3d location, int connectivity, int value) {
	Matrix3d neighbors = Matrix3d::Zero();
	int x = location(0); int y = location(1); int z = location(2);
	if(x<map3D[z].cols()-1) 
		if(map3D[z](y,x+1)==value) 
			neighbors(1,2) = 1;
	if(y<map3D[z].rows()-1)
		if(map3D[z](y+1,x)==value)
			neighbors(2,1) = 1;
	if(x>0)
		if(map3D[z](y,x-1)==value)
			neighbors(1,0) = 1;
	if(y>0)
		if(map3D[z](y-1,x)==value)
			neighbors(0,1) = 1;
	if (connectivity==8) {
		if(x<map3D[z].cols()-1 && y<map3D[z].rows()-1) 
			if(map3D[z](y+1,x+1)==value) 
				neighbors(2,2) = 1;
		if(y<map3D[z].rows()-1 && x>0)
			if(map3D[z](y+1,x-1)==value)
				neighbors(2,0) = 1;
		if(x>0 && y>0)
			if(map3D[z](y-1,x-1)==value)
				neighbors(0,0) = 1;
		if(y>0 && x<map3D[z].cols()-1)
			if(map3D[z](y-1,x+1)==value)
				neighbors(0,2) = 1;
	}
	return neighbors;
}

/* findEndPoints()
 * ----------------------------------------------------------------------------
 * This finds start and end points on each layer.  The trick here is that these
 * points will occur on the border of occlusion and obstacle buffers.  Cells
 * labeled as BUFFER_OBS will have exactly one 4-connected neighbor who is 
 * OCCLUSION if it is a start or end point.
 */
MatrixXd SurfaceExplorer::findEndPoints() {
	MatrixXd endPoints(2*int(zMax),3);
	int index = 0;
	for(int i=0; i<zMax; i++) {
		bool nextLayer = false;
		for(int j=0; j<map3D[i].rows(); j++) {
			for(int k=0; k<map3D[i].cols(); k++) {
				if(map3D[i](j,k)==BUFFER_SAFE) {
					Vector3d point(k,j,i);
					Matrix3d neighbors = checkNeighbors(point,4,OCCLUSION);
					if (neighbors.sum()==1) {
						endPoints.row(index) = point;
						map3D[i](j,k) = END_POINT;
						index++;
						if(index%2==0) {
							nextLayer = true;
							break;
						}
					}
				}
			}
			if(nextLayer==true) break;
		}
	}
	return endPoints;
}

/* orderEndPoints(MatrixXd ep)
 * ----------------------------------------------------------------------------
 * For every pair of end points on each level, it should order them so the 
 * up down motion of the quad goes back and forth like an S:
 * ------------
 * |
 * ------------
 *            |
 * ------------
 */
MatrixXd SurfaceExplorer::orderEndPoints(MatrixXd ep) { // < -- FIX THIS
	MatrixXd endPoints = ep;
	Vector3d lastEndPoint = ep.block(1,0,1,3).transpose(); // 2nd point in first layer
	for(int i=zStep+2; i<2*zMax; i+=2*zStep) {
		double d1 = (ep.block(i,0,1,3).transpose()-lastEndPoint).squaredNorm();
		double d2 = (ep.block(i+1,0,1,3).transpose()-lastEndPoint).squaredNorm();
		if(d1>d2) {
			endPoints.block(i,0,1,3) = ep.block(i+1,0,1,3);
			endPoints.block(i+1,0,1,3) = ep.block(i,0,1,3);
		}
		lastEndPoint = ep.block(i+1,0,1,3).transpose();
	}

	return endPoints;
}

/* spatiallyVariantFilter(Vector3d location, int filterRadius, int value)
 * ----------------------------------------------------------------------------
 * The goal here is to look over a region of size filterRadius around a cell
 * and change all the cells to the value given if they are empty
 */
bool SurfaceExplorer::spatiallyVariantFilter(Vector3d location, int filterRadius, int value) {
	int x = location(0); int y = location(1); int z = location(2);
	for(int i=x-filterRadius; i<=x+filterRadius; i++) {
		for(int j=y-filterRadius; j<=y+filterRadius; j++) {
			if (map3D[z](j,i) == 0) {
				if (i>=0 && i<map3D[z].cols() && j>=0 && j<map3D[z].rows()) {
					map3D[z](j,i) = value;
				}
			}
		}
	}
	return true;
}

/* expandSurface()
 * ----------------------------------------------------------------------------
 * Creates a buffer around the obstacle at half the distance of the sensor
 * capability.  This will create a space for the quad to actively explore
 * without risk of hitting the obstacle.  
 */
bool SurfaceExplorer::expandSurface() {
	// buffer surfaces by half the distance the quad can view it
	int filterRadius = int(distance2Surface/2/resolution);
	for(int i=0; i<zMax; i++) {
		for(int j=0; j<map3D[i].rows(); j++) {
			for(int k=0; k<map3D[i].cols(); k++) {
				if (map3D[i](j,k) == OBSTACLE) {
					Vector3d loc(k,j,i);
					spatiallyVariantFilter(loc,filterRadius-1,BUFFER_OBS);
				}
			}
		}
	}
	// Add safe region at exactly half the distance
	for(int i=0; i<zMax; i++) {
		for(int j=0; j<map3D[i].rows(); j++) {
			for(int k=0; k<map3D[i].cols(); k++) {
				if (map3D[i](j,k) == BUFFER_OBS) {
					Vector3d loc(k,j,i);
					spatiallyVariantFilter(loc,1,BUFFER_SAFE);
				}
			}
		}
	}
	return true;
}

/* bresenhamLine(Vector3d a, Vector3d b)
 * ----------------------------------------------------------------------------
 * This is the 3D version of this algorithm adapted for use with Eigen
 * Finds the discrete occupied cells of a line in 3 dimensions, based on
 * this implementation: https://gist.github.com/yamamushi/5823518
 */
MatrixXd SurfaceExplorer::bresenhamLine(Vector3d a, Vector3d b) {
	Vector3d point = a;
	Vector3d dPoint = b-a;
	Vector3d magDPoint = dPoint.cwiseAbs();
	Vector3d dirDPoint((dPoint(0)<0) ? -1: 1, (dPoint(1)<0) ? -1: 1, (dPoint(2)<0) ? -1: 1);
	Vector3d dPoint2(int(magDPoint(0))<<1,int(magDPoint(1))<<1,int(magDPoint(2))<<1);
	Vector2d error;
	MatrixXd output(int(magDPoint.sum()+1),3);

	// cout << point << endl << dPoint << endl << magDPoint << endl << dirDPoint << endl << dPoint2 << endl;

	int j=0;
	if (magDPoint.maxCoeff()==magDPoint(0)) {
		error << dPoint2(1)-magDPoint(0), dPoint2(2)-magDPoint(0);
		for (int i=0; i<magDPoint(0); i++) {
			output.row(j) = point; j++;
			if(error(0)>0) {
				point(1) += dirDPoint(1);
				error(0) -= dPoint2(0);
			}
			if(error(1)>0) {
				point(2) += dirDPoint(2);
				error(1) -= dPoint2(0);
			}
			error(0) += dPoint2(1);
			error(1) += dPoint2(2);
			point(0) += dirDPoint(0);
		}
	}
	else if (magDPoint.maxCoeff()==magDPoint(1)) {
		error << dPoint2(0)-magDPoint(1), dPoint2(2)-magDPoint(0);
		for (int i=0; i<magDPoint(1); i++) {
			output.row(j) = point; j++;
			if(error(0)>0) {
				point(0) += dirDPoint(0);
				error(0) -= dPoint2(1);
			}
			if(error(1)>0) {
				point(2) += dirDPoint(2);
				error(1) -= dPoint2(1);
			}
			error(0) += dPoint2(0);
			error(1) += dPoint2(2);
			point(1) += dirDPoint(1);
		}
	}
	else {
		error << dPoint2(1)-magDPoint(2), dPoint2(0)-magDPoint(2);
		for(int i=0; i<magDPoint(2); i++) {
			output.row(j) = point; j++;
			if(error(0)>0) {
				point(1) += dirDPoint(1);
				error(0) -= dPoint2(2);
			}
			if(error(1)>0) {
				point(0) += dirDPoint(0);
				error(1) -= dPoint2(2);
			}
			error(0) += dPoint2(1);
			error(1) += dPoint2(0);
			point(2) += dirDPoint(2);
		}
	}
	output.row(j) = point;
	MatrixXd temp = output.topRows(j+1);
	output = temp;

	return output;
}

/* fillOccludedSpace()
 * ----------------------------------------------------------------------------
 * This algorithm iterates through the map and finds cells occluded from the 
 * sensor.  It uses rays generated by the bresenhamLine algorithm.  
 */
bool SurfaceExplorer::fillOccludedSpace() {
	for(int i=0; i<zMax; i++) {
		for(int j=0; j<map3D[i].rows(); j++) {
			for(int k=0; k<map3D[i].cols(); k++) {
				// Cell is traversable
				if(map3D[i](j,k)!=OBSTACLE) {
					Vector3d point(k,j,i);
					MatrixXd ray = bresenhamLine(observer,point);
					// Right now I trace a line for every cell
					// Could be faster if I occlude all bad cells
					// on detected obstacle
					for(int l=0; l<ray.rows(); l++) {
						// ray point is in map
						if(ray(l,0)<map3D[i].cols() && ray(l,1)<map3D[i].rows() && ray(l,2)<zMax
							&& ray(l,0)>=0 && ray(l,1)>=0 && ray(l,2)>=0) {
							// ray hits obstacle
							if(map3D[ray(l,2)](ray(l,1),ray(l,0))==OBSTACLE) {
								// Occlude the cell
								map3D[i](j,k)=OCCLUSION;
								break;
							}
						}
					}
				}
			}
		}
	}
	return true;
}

bool SurfaceExplorer::orientTheta() {
	// For each point in the plan, look at each layer and find all of the obstacles
	// Place a distance on each obstacle.  Then find the min distance and take that
	// obstacle and compute the angle between it and the path point.  Update the 
	// plan matrix
	plan.conservativeResize(plan.rows(),4);
	vector<pair<Vector2d, double> > obstaclesAndDistances;
	for(int i=0; i<plan.rows(); i++) {
		int z = plan(i,2);
		Vector2d planPoint = plan.block(i,0,1,2).transpose();
		for(int j=0; j<map3D[z].rows(); j++) {
			for(int k=0; k<map3d[z].cols(); k++) {
				if(map3D[z](j,k) == OBSTACLE) {
					Vector2d obstaclePoint(k,j);
					double distance = (planPoint-obstaclePoint).squaredNorm();
					obstaclesAndDistances.push_back(obstaclePoint,distance);
				}
			}
		}
		pair<Vector2d, double> min = obstaclesAndDistances[0];
		for(int j=1; j<obstaclesAndDistances.size(); j++) {
			if(obstaclesAndDistances[j].second<min.second) {
				min = obstaclesAndDistances[j];
			}
		}
		plan(i,3) = acos(min.dot(planPoint)/planPoint.norm()/min.norm());
	}
	return true;
}

MatrixXd SurfaceExplorer::plan() {
	fillOccludedSpace();
	expandSurface();
	MatrixXd endPoints = findEndPoints();
	MatrixXd ordered = orderEndPoints(endPoints);
	Dstar *dstar;
	dstar = new Dstar();
	Vector2d g1 = endPoints.block(0,0,1,2).transpose();
	Vector2d s1 = start.block(0,0,2,1);
	dstar->init(s1,g1);
	dstar->setMap(map3D[0]);
	dstar->replan();
	plan = dstar->getPath(0);
	int numRows = plan.rows()-1;
	int numCols = plan.cols();
	plan.conservativeResize(numRows,numCols);
	for(int z=0; z<zMax; z+=zStep) {
		Vector2d s = endPoints.block(z,0,1,2).transpose();
		Vector2d g = endPoints.block(z+1,0,1,2).transpose();
		dstar->init(s,g);
		dstar->setMap(map3D[z]);
		dstar->replan();
		MatrixXd morePlan = dstar->getPath(z);
		MatrixXd newPlan(plan.rows()+morePlan.rows(),plan.cols());
		newPlan << plan, morePlan;
		plan = newPlan;
	}
	orientTheta();

	return plan;
}

bool SurfaceExplorer::updateMap(MatrixXd locations, int value) {
	for(int i=0; i<locations.rows(); i++) {
		if(map3D[locations(i,2)](locations(i,1),locations(i,0)) != END_POINT)
			map3D[locations(i,2)](locations(i,1),locations(i,0)) = i;
	}
	return true;
}

/* getMap()
 * ----------------------------------------------------------------------------
 * Returns the 3D map of the environment.  Internally, this map is modified
 * based on what stage the planner is at.  Printing layers of this map is 
 * usefufl for debugging.  The format is a std::vector<MatrixXd> where each
 * matrix is the same size.
 */
vector<MatrixXd> SurfaceExplorer::getMap() {
	return map3D;
}

int main(int argc, char **argv) {
	vector<MatrixXd> map(10,MatrixXd::Zero(30,30));
	map[0].block(10,10,10,10) = MatrixXd::Ones(10,10);
	map[1].block(12,12,8,8) = MatrixXd::Ones(8,8);
	map[2].block(12,12,7,7) = MatrixXd::Ones(7,7);
	map[3].block(13,12,6,7) = MatrixXd::Ones(6,7);
	map[4].block(15,15,2,2) = MatrixXd::Ones(2,2);
	map[5].block(15,15,2,2) = MatrixXd::Ones(2,2);
	map[6].block(15,15,2,2) = MatrixXd::Ones(2,2);
	map[7].block(15,15,2,2) = MatrixXd::Ones(2,2);
	map[8].block(15,15,2,2) = MatrixXd::Ones(2,2);
	map[9].block(16,16,1,1) = MatrixXd::Ones(1,1);

	double resolution = 1; double zMax = 10; double zStep = 2; double d2s = 4;
	SurfaceExplorer se = SurfaceExplorer(resolution,d2s,zStep,zMax);
	Vector3d start(1,0,0);
	Vector3d observer(15,0,0);
	se.init(start,observer,map);
	vector<MatrixXd> initial = se.getMap();

	// cout << initial[0] << endl << initial[1] << endl << initial[2] << endl << endl;

	// se.fillOccludedSpace();
	// vector<MatrixXd> result = se.getMap();

	// cout << result[0] << endl << result[1] << endl << result[2] << endl << endl;

	// se.expandSurface();
	// vector<MatrixXd> expanded = se.getMap();

	// cout << expanded[0] << endl << expanded[1] << endl << expanded[2] << endl << endl;

	// MatrixXd n = se.findEndPoints();
	// vector<MatrixXd> endPoints = se.getMap();

	// cout << endPoints[0] << endl << endPoints[1] << endl << endPoints[2] << endl << endl;

	// MatrixXd o = se.orderEndPoints(n);

	// cout << o << endl;

	MatrixXd path = se.plan();
	cout << path << endl;

	se.updateMap(path.block(0,0,path.rows(),3),PATH_POINT);
	se.updateMap(start.transpose(),END_POINT);
	vector<MatrixXd> result = se.getMap();
	for(int i=0; i< result.size(); i++) {
		cout << result[i] << endl;
	}

	return 1;
}