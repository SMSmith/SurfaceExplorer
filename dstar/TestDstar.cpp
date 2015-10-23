#include "Dstar.h"
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <fstream>

namespace po = boost::program_options; 

Dstar *dstar;
bool loadMap(string fileName, MatrixXd& map);

/* main(int argc, char **argv)
 * -----------------------------------------------------------
 * The main point of this program is to test the capabilities
 * of the D* lite implementation. Generate a map file of 
 * ASCII 1's and 0's with 1's being obstacles and you can 
 * test whether it produces optimal solutions from start
 * to goal.
 */
int main(int argc, char **argv) {
    // Handle Arguments
    po::options_description desc("This serves to test the dstar implementation");
    desc.add_options()
        ("help,h", "See the options below")
        ("map,m", po::value<string>(), "Pass in the binary file for the map")
        ("start,s", po::value<vector<int> >()->multitoken(), "give it a start position '-s 3 13'")
        ("end,e", po::value<vector<int> >()->multitoken(), "give it an end position '-e 23 0'");
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    // Start point as x y Eigen vector
    vector<int> startPoint;
    Vector2d start(0,0);
    if (!vm["start"].empty() && (startPoint = vm["start"].as<vector<int> >()).size()==2) {
        start(0) = startPoint[0];
        start(1) = startPoint[1];
    }
    // End point as xy Eigen vector
    vector<int> endPoint;
    Vector2d goal(1,1);
    if (!vm["end"].empty() && (endPoint = vm["end"].as<vector<int > >()).size()==2) {
        goal(0) = endPoint[0];
        goal(1) = endPoint[1];
    }
    // The map as an ASCII file - Could be binary, but then harder to play with in testing
    MatrixXd map(0,0);
    if (!vm["map"].empty()) {
        loadMap(vm["map"].as<string>(),map);
    }

    // Call the dstar algorithm and run
    dstar = new Dstar();
    dstar->init(start,goal);
    dstar->setMap(map);
    dstar->replan();
    cout << map << endl;
    dstar->printPath();

    return 1;
}

/* bool loadMap(string fileName, MatrixXd& map)
 * ----------------------------------------------
 * Loads a map into an eigen matrix from file
 * The map is ascii format, 0 is free, 1 is 
 * obstacle.  Binary would be better, but hard
 * to hand code.
 */
bool loadMap(string fileName, MatrixXd& map) {
    vector<string> lines;
    ifstream mapFile(fileName.c_str());
    while(!mapFile.eof()) {
        string line;
        getline(mapFile, line);
        lines.push_back(line);
    }

    map.resize(lines.size(),lines[0].length());
    for(int i=0;i<lines.size(); i++) {
        for(int j=0; j<lines[0].length(); j++) {
            map(i,j) = ((int) (lines[i].at(j)-'0'));
        }
    }

    return true;
}