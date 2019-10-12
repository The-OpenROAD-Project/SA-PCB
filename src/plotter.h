
#ifndef PCBROUTER_PLOTTER_H
#define PCBROUTER_PLOTTER_H

#include <iostream>
#include <stdio.h>
#include <fstream>
#include <string>
#include <math.h>
#include <iomanip>
#include <list>
//#include "point.h"
//#include "rectangle.h"
//#include "segment.h"
//#include "path.h"

using namespace std;

// TODO
// Plot the via/wire obstacle grids
class Plotter
{
public:
  Plotter(char *outputFileName);
  ~Plotter();

  bool isInit();

  void outputString(const string &s);
  void outputDouble(double d, size_t digit = 2);
  //void outputString(char* s);

  //plot
  void plotLine(double x1, double y1, double x2, double y2);
  //void plotLine(Point<double> p1, Point<double> p2);
  //void plotLine(Segment s);
  void plotRec(double &x1, double &y1, double &x2, double &y2);
  void plotRec(double x1, double y1, double x2, double y2);
  //void plotRec(Point<double> p1, Point<double> p2);
  //void plotRec(Rectangle r);
  //void plotSquare(Point<double> cp, double length);
  //void plotOct(Point<double> cp, double &radius);
  //void plotPath(Path &p);
  //void plotRectilinear(list<Point<double> > &rectilinear);

  //plot Obj
  //void plotCircleObj(Point<double> cp, double radius, const string &rgb = "\"#000000\"");

private:
  //void outputPoint(Point<double> p);

private:
  ofstream _outFile;
  //size_t _objCounter;
};

//Congestion color map
//Color	    	Congestion  Value
//Blue	    	+1		3
//Green	    	0		2
//Light Brown   -1		8
//Red	    	-2		1
//Magenta	    -3		4
//Light Blue    -4 or less	5

//	string colorTable[20] = {"7FFFD4", "FFE4C4", "8A2BE2", "7FFF00", "D2691E",
//                                 "6495ED", "B8860B", "006400", "556B2F", "CD5C5C",
//                                 "F0E68C", "90EE90"};

//	string colorTable2[20] = {"7FFFB4", "FFE4A4", "8A2BC2", "7FDF00", "D2690E",
//                                 "6495CD", "B8660B", "004400", "556B0F", "CD5C3C",
//                                 "F0E66C", "90EE70"};

//'-' w l lt rgb \"#006000\"

#endif
