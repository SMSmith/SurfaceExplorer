#ifndef SURFACE_EXPLORER_H
#define SURFACE_EXPLORER_H

#include <Eigen>
#include <vector>
#include <iostream>
#include <math.h>

using namespace Eigen;
using namespace std;

#define OCCLUSION 3
#define OBSTACLE 1
#define BUFFER_OBS 2
#define BUFFER_SAFE 4
#define END_POINT 5
#define PATH_POINT 7

class SurfaceExplorer {
	public:
		SurfaceExplorer(double r, double d2s, double zs, double zm);
		bool 				init(Vector3d s, Vector3d obs, vector<MatrixXd>& m);
		vector<MatrixXd> 	getMap();
		bool				fillOccludedSpace();
		bool				expandSurface();
		MatrixXd 			findEndPoints();
		MatrixXd			orderEndPoints(MatrixXd ep);
		MatrixXd			plan();
		bool				updateMap(MatrixXd locations, int value);

		MatrixXd plan;

	private:
		vector<MatrixXd> map3D;
		vector<MatrixXd> expandedMap3D;
		Vector3d start;
		Vector3d observer;
		double resolution, zMax, zStep, distance2Surface;

		MatrixXd			bresenhamLine(Vector3d a, Vector3d b);
		bool				spatiallyVariantFilter(Vector3d loc, int filterRadius, int value);
		Matrix3d 			checkNeighbors(Vector3d location, int connectivity, int value);
};

class strangeMapException: public exception {
	virtual const char* what() const throw() {
		return "Your map has inconsistant matrix sizes on each layer";
	}
} mapException;

#endif