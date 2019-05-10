#include <iostream>
#include <vector>
#include <fstream>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/algorithms/area.hpp>
#include <string>
#include <iostream>
#include <sstream>
#include <regex>
#include <map>
#include <cmath>
#include <stdio.h>
#include <unistd.h>
#include<ctype.h>
#include<stdlib.h>
#include "readFiles.h"

void printEverything();
void initialPlacement();
void printCellMap();
void printRowMap();
void updateCenter();
void CalcBoundaries();
int macroPlacement();
double cost(int temp_debug = 0);
double cellOverlap();
double wireLength();
float timberWolfAlgorithm();
void update_Temperature();
void initiateMove();
bool checkMove(long int prevCost);
void crandomPlacement(int xmin, int xmax, int ymin, int ymax, Node n);
map < string, Node > ::iterator random_node();


struct boundaries {
  int minX, maxX, minY, maxY;
};
