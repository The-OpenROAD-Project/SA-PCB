#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctype.h>
#include <map>
#include <vector>
#include <string>
#include <numeric>
#include <cmath>
#include <unistd.h>
#include <ctime>
#include <ratio>
#include <chrono>

// boost version 1.69.0
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/geometry.hpp>
#include <boost/assign.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/io/io.hpp>
#include <boost/geometry/algorithms/area.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/foreach.hpp>

//#include "taskflow/taskflow.hpp"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>


#include "readFiles.h"
#include "time.h"
#include "readScl.h"


void printEverything();
void initialPlacement();
void printCellMap();
void printRowMap();
void updateCenter();
void CalcBoundaries();
void SetInitParameters();
int macroPlacement();
void validateMove(Node* node, double rx, double ry);
double cost(int temp_debug = 0);
double cost_partial(int temp_debug, map < string, Node > &nodes);
double cellOverlap();
double wireLength();
double cellOverlap_partial(map < string, Node > &nodes);
double wireLength_partial(map < string, Node > &nodes);
double rudy();
float timberWolfAlgorithm(int outer_loop_iter, int inner_loop_iter, double eps, double t_0,bool var);
float multistart();
double varanelli_cohoon();
void update_Temperature();
double initiateMove();
bool checkMove(double prevCost);
void project_soln();
void randomPlacement(int xmin, int xmax, int ymin, int ymax, Node n);
void gen_report(map<string, vector<double> > &report);
map < string, Node > ::iterator random_node();

struct boundaries {
  double minX, maxX, minY, maxY = 0.0;
};
