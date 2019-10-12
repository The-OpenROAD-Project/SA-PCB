#ifndef PCBROUTER_GLOBALPARAM_H
#define PCBROUTER_GLOBALPARAM_H

#include <string>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include "util.h"

using namespace std;

#define PI 3.14159265358979323846264338

#define BOUND_CHECKS

// TODO
// See router.h -> class net
struct NetClass
{

};

// TODO
// See clearance matrix in the EAGLE/KiCad format
class GlobalParam
{
public:
  static int gLayerNum;
  static double gEpsilon;
  static float HPWLWeight;
  static float DensityWeight;
  static float RoutabilityWeight;
  static bool g90DegreeMode;

  //Outputfile
  static int gOutputPrecision;
  static string gOutputFolder;

  //Log
  static string gLogFolder;

  const static double gSqrt2;
  const static double gTan22_5;

  static int gSeed;

  static util::TimeUsage runTime;

  static void setFolders();
  static void setLayerNum(int l) { gLayerNum = l; }

  static void showParam();
  static void showCurrentUsage(const string comment);
  static void showFinalUsage(const string comment);
  static void setUsageStart();
};

#endif
