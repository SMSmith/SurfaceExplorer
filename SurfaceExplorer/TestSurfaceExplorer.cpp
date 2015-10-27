#include "SurfaceExplorer.h"
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <fstream>

namespace po = boost::program_options;

bool loadMap3D(string fileName, vector<MatrixXd>& map);

/* main(int argc, char ** argv)
 * ----------------------------------------------------------------------------
 * Example usage on some fake map data.  The data has 10 levels with one 
 * obstacle that diminishes in size as you go up the layers (kind of like a 
 * cell phone tower).  
 */
int main(int argc, char **argv) {
	po::options_description desc("This serves to test SurfaceExplorer by giving it different initial conditions");
	desc.add_options()
		("help,h", "See the options below")
		("map,m", po::value<string>(), "Provide a map file, example in repo")
		("start,s", po::value<vector<int> >()->multitoken(), "Provide the start point")
		("observer,o", po::value<vector<int> >()->multitoken(), "Provide the sensor's observation point")
		("z_max,z", po::value<double>(), "Provide the maximum physical height")
		("z_step,i", po::value<double>(), "Provide the physical step height")
		("resolution,r", po::value<double>(), "Provide the distance/cell")
		("sensor_max,d", po::value<double>(), "Provide the maximum distance to surface")
		;
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

	if(!vm["help"].empty()) {
		cout << desc;
		return 0;
	}

	// Start point as x y z Eigen vector
	Vector3d start(0,0,0);
	vector<int> startPoint;
	if(!vm["start"].empty() && (startPoint=vm["start"].as<vector<int> >()).size()==3) {
		start(0) = startPoint[0];
		start(1) = startPoint[1];
		start(2) = startPoint[2];
	}
	// Observer point as x y z Eigen vector
	Vector3d observer(0,0,0);
	vector<int> observerPoint;
	if(!vm["start"].empty() && (observerPoint=vm["start"].as<vector<int> >()).size()==3) {
		observer(0) = observerPoint[0];
		observer(1) = observerPoint[1];
		observer(2) = observerPoint[2];
	}
	// zMax value from arguments
	double zMax = 12;
	if(!vm["z_max"].empty()) {
		zMax = vm["z_max"].as<double>();
	}
	// zStep value from arguments
	double zStep = 2;
	if(!vm["z_step"].empty()) {
		zStep = vm["z_step"].as<double>();
	}
	// distance to surface value from arguments
	double d2s = 4;
	if(!vm["sensor_max"].empty()) {
		d2s = vm["sensor_max"].as<double>();
	}
	// resolution value from arguments
	double resolution = 1;
	if(!vm["resolution"].empty()) {
		resolution = vm["resolution"].as<double>();
	}
	// Load the map - example in repo
	vector<MatrixXd> map(10,MatrixXd::Zero(30,30));
	if(!vm["map"].empty()) {
		loadMap3D(vm["map"].as<string>(),map);
	}
	else {
		// Create a map
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
	}

	// The parameters for the map - note zMax can be higher than the map's
	// highest layer.  
	SurfaceExplorer se = SurfaceExplorer(resolution,d2s,zStep,zMax);
	if(!se.init(start,observer,map)) {
		return 0;
	}

	// Plan the path and convert orientation to degrees for readability
	MatrixXd path = se.planPath();
	for(int i=0; i<path.rows(); i++) {
		path(i,3) = path(i,3)*180/M_PI;
	}
	cout << path << endl;

	// Draw the path points on the map for visualization
	se.updateMapByLocation(path.block(0,0,path.rows(),3),PATH_POINT);
	// se.updateMapByLocation(start.transpose(),END_POINT);

	// Draw the map layers to show the decision path
	vector<MatrixXd> result = se.getMap();
	for(int z=0; z<result.size(); z++) {
		cout << result[z] << endl;
	}

	return 1;
}

/* loadmap(string fileName, vector<MatrixXd> map)
 * ----------------------------------------------------------------------------
 * loads a map from file - see repo for example
 */
bool loadMap3D(string fileName, vector<MatrixXd>& map) {
	vector<vector<string> > lines;
	ifstream mapFile(fileName.c_str());
	vector<string> layer;
	while(!mapFile.eof()) {
		string line;
		getline(mapFile,line);
		if(line.at(0)==',') {
			lines.push_back(layer);
			layer.clear();
		}
		else layer.push_back(line);
	}

	for(int z=0; z<lines.size(); z++) {
		for(int y=0; y<lines[z].size(); y++) {
			for(int x=0; x<lines[z][y].length(); x++) {
				map[z](y,x) = ((int) (lines[z][y].at(x)-'0'));
			}
		}
	}
	return true;
}