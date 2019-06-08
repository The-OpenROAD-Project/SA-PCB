#include <map>
#include "readFiles.h"

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
double cost_partial(int temp_debug, map < string, Node > nodes);
double cellOverlap();
double wireLength();
double cellOverlap_partial(map < string, Node > nodes);
double wireLength_partial(map < string, Node > nodes);
double rudy();
float timberWolfAlgorithm(int outer_loop_iter, int inner_loop_iter, double eps, double t_0,bool var);
float multistart();
double varanelli_cohoon();
void update_Temperature();
double initiateMove();
bool checkMove(double prevCost);
void project_soln();
void randomPlacement(int xmin, int xmax, int ymin, int ymax, Node n);
void gen_report(map<string, vector<double> > report);
map < string, Node > ::iterator random_node();

struct boundaries {
  double minX, maxX, minY, maxY = 0.0;
};
