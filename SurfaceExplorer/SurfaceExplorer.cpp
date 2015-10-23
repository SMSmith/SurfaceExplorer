// #include "Dstar.h"
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
				if(map3D[i](j,k)==BUFFER_OBS) {
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

MatrixXd SurfaceExplorer::findPathPoints() {
	MatrixXd pathPoints(0,0);
	return pathPoints;
}

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

bool SurfaceExplorer::expandSurface() {
	// buffer surfaces by half the distance the quad can view it
	int filterRadius = int(distance2Surface/2/resolution);
	for(int i=0; i<zMax; i++) {
		for(int j=0; j<map3D[i].rows(); j++) {
			for(int k=0; k<map3D[i].cols(); k++) {
				if (map3D[i](j,k) == OBSTACLE) {
					Vector3d loc(k,j,i);
					spatiallyVariantFilter(loc,filterRadius,BUFFER_OBS);
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
	vector<MatrixXd> map(4,MatrixXd::Zero(30,30));
	// map[0] << 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0;
	// map[1] << 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0;
	// map[2] << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0;
	// map[3] << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
	map[0].block(4,4,5,5) = MatrixXd::Ones(5,5);
	map[1].block(4,5,4,4) = MatrixXd::Ones(4,4);
	map[2].block(5,6,2,2) = MatrixXd::Ones(2,2);
	map[3].block(7,7,1,1) = MatrixXd::Ones(1,1);

	double resolution = 1; double zMax = 3; double zStep = 1; double d2s = 2;
	SurfaceExplorer se = SurfaceExplorer(resolution,d2s,zStep,zMax);
	Vector3d start(1,0,0);
	Vector3d observer(0,0,0);
	se.init(start,observer,map);
	vector<MatrixXd> initial = se.getMap();

	cout << initial[0] << endl << initial[1] << endl << initial[2] << endl << endl;

	se.fillOccludedSpace();
	vector<MatrixXd> result = se.getMap();

	cout << result[0] << endl << result[1] << endl << result[2] << endl << endl;

	se.expandSurface();
	vector<MatrixXd> expanded = se.getMap();

	cout << expanded[0] << endl << expanded[1] << endl << expanded[2] << endl << endl;

	MatrixXd n = se.findEndPoints();
	vector<MatrixXd> endPoints = se.getMap();

	cout << endPoints[0] << endl << endPoints[1] << endl << endPoints[2] << endl << endl;

	return 1;s
}