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
#include <boost/geometry/index/rtree.hpp>

//#include "taskflow/taskflow.hpp"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>


#include "readFiles.h"
#include "time.h"
//#include "readScl.h"


void printEverything();
void initialPlacement();
void printCellMap();
void printRowMap();
void updateCenter();
void CalcBoundaries();
void SetInitParameters(std::pair <double,double> *wl_normalization,
                       std::pair <double,double> *area_normalization,
                       std::pair <double,double> *routability_normalization,
                       map<int, vector<Pin> > &netToCell);
int macroPlacement();
void validateMove(Node* node, double rx, double ry);
double cost(
            std::pair <double,double> &wl_normalization,
            std::pair <double,double> &area_normalization,
            std::pair <double,double> &routability_normalization,
            map<int, vector<Pin> > &netToCell,
            int temp_debug = 0);
double cost_partial(int temp_debug,
                    vector < Node > &nodes,
                    std::pair <double,double> &wl_normalization,
                    std::pair <double,double> &area_normalization,
                    std::pair <double,double> &routability_normalization,
                    map<int, vector<Pin> > &netToCell);
double cellOverlap();
double wireLength(map<int, vector<Pin> > &netToCell);
double cellOverlap_partial(vector < Node > &nodes);
double wireLength_partial(vector < Node > &nodes, map<int, vector<Pin> > &netToCell);
double rudy(map<int, vector<Pin> > &netToCell);
float timberWolfAlgorithm(int outer_loop_iter, int inner_loop_iter, double eps, double t_0,bool var, map<int, vector<Pin> > &netToCell);
float multistart();
double varanelli_cohoon(vector< int > &accept_history,
                        double & Temperature,
                        std::pair <double,double> &wl_normalization,
                        std::pair <double,double> &area_normalization,
                        std::pair <double,double> &routability_normalization,
                        map<int, vector<Pin> > &netToCell);
void update_Temperature(double* Temperature);
double initiateMove(vector< int > *accept_history,
                    double & Temperature,
                    std::pair <double,double> &wl_normalization,
                    std::pair <double,double> &area_normalization,
                    std::pair <double,double> &routability_normalization,
                    map<int, vector<Pin> > &netToCell);
bool checkMove(double prevCost,
               vector< int > *accept_history,
               double & Temperature,
               std::pair <double,double> &wl_normalization,
               std::pair <double,double> &area_normalization,
               std::pair <double,double> &routability_normalization,
               map<int, vector<Pin> > &netToCell);
void project_soln();
void randomPlacement(int xmin, int xmax, int ymin, int ymax, Node n);
void gen_report(map<string, vector<double> > &report,
                vector< double > &accept_ratio_history,
                std::pair <double,double> &wl_normalization,
                std::pair <double,double> &area_normalization,
                std::pair <double,double> &routability_normalization,
                map<int, vector<Pin> > &netToCell);
void update_accept_history(vector< int > &accept_history, vector< double > *accept_ratio_history, float *accept_ratio);
void update_rtree(int idx);
vector < Node > ::iterator random_node();

struct boundaries {
  double minX, maxX, minY, maxY = 0.0;
};
