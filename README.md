# SurfaceExplorer
A planning library for exploring surfaces in 3D while in view of a sensor

## Usage
Initiate a SurfaceExplorer by providing the following items:
* resolution (r) - the physical size of a cell (distance/cell) in whatever units you use (keep them consistant)
* distance2Surface (d2s) - the maximum physical distance the quad's sensor can see the surface
* zStep (zs) - the physical distance the quad should advance between 2d scans
* zMax (zm) - the physical max height the quad should reach

You can initialize the explorer by passing the following items to init:
* start (s) - Vector3d which is the 3 dimensional start point (assumes it starts at the 0th layer)
* observer (o) - Vector3d which is the 3 dimensional point corresponding to the location of a sensor that the quad must stay in line-of-sight with
* map (m) - vector<MatrixXd> the map which is a vector af matrices corresponding to the locations of the obstacles (0 is free, 1 is obstacle)

Iniating the SurfaceExplorer and calling its init function allows you to call planPath() which will generate the 3D plan from start point, through imaging the structure and back to the start point, while staying in sight of the sensor. planPath will return a MatrixXd, size Nx4 where N is the number of points in the path and each row corresponds to a point in (x,y,z,theta). Theta is in radians.

An example workflow is seen in TestSurfaceExplorer.cpp

Call it like this:
./surfaceExplorer -s 10 3 2 -o 0 0 1 -m SurfaceExplorer/map3D.txt