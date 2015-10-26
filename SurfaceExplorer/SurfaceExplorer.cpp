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

	// Make a bunch of dstar planners for each level of the map
	// The can be accessed individually to alter stuff
	dstar =  new Dstar*[int(zMax+1)];
	for(int z=0; z<zMax+1; z++) {
		dstar[z] = new Dstar();
	}

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
	for(int z=0; z<zMax; z++) {
		bool nextLayer = false;
		for(int y=0; y<map3D[z].rows(); y++) {
			for(int x=0; x<map3D[z].cols(); x++) {
				if(map3D[z](y,x)==BUFFER_SAFE) {
					Vector3d point(x,y,z);
					Matrix3d neighbors = checkNeighbors(point,4,OCCLUSION);
					if (neighbors.sum()>0) {
						endPoints.row(index) = point.transpose();
						map3D[z](y,x) = END_POINT;
						index++;
						if(index%2==0) {
							nextLayer = true;
							break;
						}
					}
				}
			}
			if(nextLayer) break;
		}
		// Handles semi-occluded layer cases
		if(!nextLayer && z>0) {
			endPoints.row(index) = endPoints.row(index-2*zStep);
			endPoints(index,2)+=zStep;
			index++;
			endPoints.row(index) = endPoints.row(index-2*zStep);
			endPoints(index,2)+=zStep;
			index++;
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
MatrixXd SurfaceExplorer::orderEndPoints(MatrixXd ep) {
	MatrixXd endPoints = ep;
	Vector3d lastEndPoint = ep.block(1,0,1,3).transpose(); // 2nd point in first layer
	for(int z=2*zStep; z<2*zMax; z+=2*zStep) {
		double d1 = (ep.block(z,0,1,3).transpose()-lastEndPoint).squaredNorm();
		double d2 = (ep.block(z+1,0,1,3).transpose()-lastEndPoint).squaredNorm();
		// cout << d1 << " " << d2 << ep.block(z,0,1,3) << endl;
		if(d1>d2) {
			endPoints.block(z,0,1,3) = ep.block(z+1,0,1,3);
			endPoints.block(z+1,0,1,3) = ep.block(z,0,1,3);
		}
		lastEndPoint = endPoints.block(z+1,0,1,3).transpose();
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
	for(int z=0; z<zMax; z++) {
		for(int y=0; y<map3D[z].rows(); y++) {
			for(int x=0; x<map3D[z].cols(); x++) {
				if (map3D[z](y,x) == OBSTACLE) {
					Vector3d loc(x,y,z);
					spatiallyVariantFilter(loc,filterRadius-1,BUFFER_OBS);
				}
			}
		}
	}
	// Add safe region at exactly half the distance
	for(int z=0; z<zMax; z++) {
		for(int y=0; y<map3D[z].rows(); y++) {
			for(int x=0; x<map3D[z].cols(); x++) {
				if (map3D[z](y,x) == BUFFER_OBS) {
					Vector3d loc(x,y,z);	
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
	for(int z=0; z<zMax; z++) {
		for(int y=0; y<map3D[z].rows(); y++) {
			for(int x=0; x<map3D[z].cols(); x++) {
				// Cell is traversable
				if(map3D[z](y,x)!=OBSTACLE) {
					Vector3d point(x,y,z);
					MatrixXd ray = bresenhamLine(observer,point);
					// Right now I trace a line for every cell
					// Could be faster if I occlude all bad cells
					// on detected obstacle
					for(int l=0; l<ray.rows(); l++) {
						// ray point is in map
						if(ray(l,0)<map3D[z].cols() && ray(l,1)<map3D[z].rows() && ray(l,2)<zMax
							&& ray(l,0)>=0 && ray(l,1)>=0 && ray(l,2)>=0) {
							// ray hits obstacle
							if(map3D[ray(l,2)](ray(l,1),ray(l,0))==OBSTACLE) {
								// Occlude the cell
								map3D[z](y,x)=OCCLUSION;
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
	MatrixXd mins(2,plan.rows());
	for(int i=0; i<plan.rows(); i++) {
		int z = plan(i,2);
		vector<pair<Vector2d, double> > obstaclesAndDistances;
		Vector2d planPoint = plan.block(i,0,1,2).transpose();
		for(int y=0; y<map3D[z].rows(); y++) {
			for(int x=0; x<map3D[z].cols(); x++) {
				if(map3D[z](y,x) == OBSTACLE) {
					Vector2d obstaclePoint(y,x);
					double distance = (planPoint-obstaclePoint).squaredNorm();
					pair<Vector2d,double> newObs(obstaclePoint,distance);
					obstaclesAndDistances.push_back(newObs);
				}
			}
		}
		pair<Vector2d, double> min = obstaclesAndDistances[0];
		for(int j=1; j<obstaclesAndDistances.size(); j++) {
			if(obstaclesAndDistances[j].second<min.second) {
				min = obstaclesAndDistances[j];
			}
		}
		mins.block(0,i,2,1) = Vector2d(min.first);
		Vector2d xAxis(1,0);
		double signAngle = (min.first(1)<planPoint(1))? 1:-1;
		plan(i,3) = signAngle*acos(xAxis.dot(min.first-planPoint)/(min.first-planPoint).norm());
	}
	// cout << endl << mins << endl;
	return true;
}

MatrixXd SurfaceExplorer::planPath() {
	fillOccludedSpace();
	expandSurface();
	MatrixXd endPoints = findEndPoints();
	// cout << endPoints << endl;
	MatrixXd ordered = orderEndPoints(endPoints);
	// cout << ordered << endl;
	// dstar[0] = new Dstar();
	Vector2d g1 = ordered.block(0,0,1,2).transpose();
	Vector2d s1 = start.block(0,0,2,1);
	dstar[0]->init(s1,g1);
	dstar[0]->setMap(map3D[0]);
	dstar[0]->replan();
	plan = dstar[0]->getPath(0);
	int numRows = plan.rows()-1;
	int numCols = plan.cols();
	plan.conservativeResize(numRows,numCols);
	for(int z=0; z<zMax; z+=zStep) {
		Vector2d s = ordered.block(2*z,0,1,2).transpose();
		Vector2d g = ordered.block(2*z+1,0,1,2).transpose();
		// if (z>0) dstar[z] = new Dstar();
		dstar[z]->init(s,g);
		dstar[z]->setMap(map3D[z]);
		if(!dstar[z]->replan()) {
			cout << "Cannot generate plan for layer " << z << endl;
		}
		MatrixXd morePlan = dstar[z]->getPath(z);
		MatrixXd endPlan(0,0);
		// Connect the layers with a straight line (assume no obstacles)
		if(z+zStep<zMax) {
			Vector3d lastGoal = ordered.block(2*z+1,0,1,3).transpose();
			Vector3d nextStart = ordered.block(2*(z+zStep),0,1,3).transpose();
			endPlan = bresenhamLine(lastGoal,nextStart);
			// Remove duplicate start and ends
			endPlan = endPlan.block(1,0,endPlan.rows()-2,3);
		}
		// Return to home
		else {
			dstar[z+1]->init(g,start.block(0,0,2,1));
			dstar[z+1]->setMap(map3D[z]);
			if(!dstar[z+1]->replan()) {
				cout << "Dstar failed to plan to home" << endl;
			}
			MatrixXd lastDStar = dstar[z+1]->getPath(z);
			Vector3d aboveStart3d = start;
			aboveStart3d(2) = z;
			// Lower the quad
			MatrixXd moreEndPlan = bresenhamLine(aboveStart3d,start);
			MatrixXd temp = moreEndPlan.block(1,0,moreEndPlan.rows()-1,3);
			moreEndPlan = temp;
			endPlan.conservativeResize(lastDStar.rows()+moreEndPlan.rows(),lastDStar.cols());
			endPlan << lastDStar, moreEndPlan;
			temp = endPlan.block(1,0,endPlan.rows()-1,3);
			endPlan = temp;
		}
		MatrixXd newPlan(plan.rows()+morePlan.rows()+endPlan.rows(),plan.cols());
		newPlan << plan, morePlan, endPlan;
		plan = newPlan;
	}
	// Focus the sensor to the nearest obstacle
	orientTheta(); 	

	return plan;
}

/* updateMapByLocation(MatrixXd locations, int value)
 * ----------------------------------------------------------------------------
 * Allows specific locations on the map to be labeled differently
 */
bool SurfaceExplorer::updateMapByLocation(MatrixXd locations, int value) {
	for(int i=0; i<locations.rows(); i++) {
		// if(map3D[locations(i,2)](locations(i,1),locations(i,0)) != END_POINT)
			map3D[locations(i,2)](locations(i,1),locations(i,0)) = i;
	}
	return true;
}

/* updateMap(vector<MatrixXd> map)
 * ----------------------------------------------------------------------------
 * This allows the user to update the entire map.  It also updates all of the 
 * layers in the dStars.  This means the dStars can be asked to replan on the 
 * new map --> see replan().
 */
bool SurfaceExplorer::updateMap(vector<MatrixXd> map) {
	map3D = map;
	for(int z=0; z<zMax; z++) {
		dstar[z]->setMap(map3D[z]);
	}
	return true;
}

/* updateMapLayer(int layer, MatrixXd mapLayer)
 * ----------------------------------------------------------------------------
 * This lets the user update just one layer of the map
 */
bool SurfaceExplorer::updateMapLayer(int layer, MatrixXd mapLayer) {
	map3D[layer] = mapLayer;
	dstar[layer]->setMap(mapLayer);
	return true;
}

/* updateMapCell(Vector3d cell, int value)
 * ----------------------------------------------------------------------------
 * Update individual cells on the map. A 1 represents an obstacle. 0 is free.
 * Useful when called before a replanLayer().  
 */
bool SurfaceExplorer::updateMapCell(Vector3d cell, int value) {
	map3D[int(cell(2))](cell(1),cell(0)) = value;
	if(value==1)
		value=-1;
	dstar[int(cell(2))]->updateCell(cell.block(0,0,2,1),value);
	return true;
}

/* replanAll()
 * ----------------------------------------------------------------------------
 * Lets the user replan all layers of the explorer.  Usually applied after
 * significan map updates.  
 */
bool SurfaceExplorer::replanAll() {
	for(int z=0; z<zMax+1; z++) {
		if(!dstar[z]->replan()) {
			cout << "Replanning failed at " << z << "th layer" << endl;
			return false;
		}
	}
	return true;
}
/* replanLayer(int Layer)
 * ----------------------------------------------------------------------------
 * Lets the user replan a specific layer of the explorer.  Usually applied
 * after layer updates.  
 */
 bool SurfaceExplorer::replanLayer(int layer) {
 	if(!dstar[layer]->replan()) {
 		cout << "Replanning " << layer << "th layer failed" << endl;
 		return false;
 	}
 	return true;
 }
/* getDStar(int layerHeight)
 * ----------------------------------------------------------------------------
 * Allows the user to get access to planners at specific levels.  See Dstar.h
 * for public calls on these objects.  
 */
Dstar* SurfaceExplorer::getDStar(int layerHeight) {
	return dstar[layerHeight];
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

	double resolution = 1; double zMax = 12; double zStep = 2; double d2s = 4;
	SurfaceExplorer se = SurfaceExplorer(resolution,d2s,zStep,zMax);
	Vector3d start(1,0,0);
	Vector3d observer(0,10,0);
	se.init(start,observer,map);
	vector<MatrixXd> initial = se.getMap();

	MatrixXd path = se.planPath();
	for(int i=0; i<path.rows(); i++) {
		path(i,3) = path(i,3)*180/M_PI;
	}
	cout << path << endl;

	se.updateMapByLocation(path.block(0,0,path.rows(),3),PATH_POINT);
	se.updateMapByLocation(start.transpose(),END_POINT);
	vector<MatrixXd> result = se.getMap();
	for(int z=0; z<result.size(); z++) {
		cout << result[z] << endl;
	}

	return 1;
}