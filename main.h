#include <map>
#include "readFiles.h"

void printEverything();
void initialPlacement();
void printCellMap();
void printRowMap();
void updateCenter();
void CalcBoundaries();
int macroPlacement();
void validateMove(Node* node, double rx, double ry);
double cost(int temp_debug = 0);
double cellOverlap();
double wireLength();
float timberWolfAlgorithm();
float multistart();
void update_Temperature();
void initiateMove();
bool checkMove(long int prevCost);
void randomPlacement(int xmin, int xmax, int ymin, int ymax, Node n);
map < string, Node > ::iterator random_node();


struct boundaries {
  int minX, maxX, minY, maxY;
};
