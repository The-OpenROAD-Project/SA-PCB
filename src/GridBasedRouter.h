// GridBasedRouter.h
#ifndef PCBROUTER_GRID_BASED_ROUTER_H
#define PCBROUTER_GRID_BASED_ROUTER_H

#include <fstream>
#include <vector>
#include <string>
#include <cstdio>
#include <iostream>
#include "BoardGrid.h"
#include "kicadPcbDataBase.h"
#include "globalParam.h"
#include "util.h"

class GridBasedRouter
{
public:
  //ctor
  GridBasedRouter(kicadPcbDataBase &db) : mDb(db) {}
  //dtor
  ~GridBasedRouter(){}

  void test_router();
  bool outputResults2KiCadFile(std::vector<MultipinRoute> &nets);

private:
  bool writeNets(std::vector<MultipinRoute> &multipinNets, std::ofstream &ofs);

  // Utility
  int dbLengthToGridLength(const double dbLength) { return (int)ceil(dbLength * inputScale); }

  bool dbPointToGridPoint(const point_2d &dbPt, point_2d &gridPt);
  bool gridPointToDbPoint(const point_2d &gridPt, point_2d &dbPt);
  void addPinCost(const pin &, const float);
  void addPinCost(const padstack &, const instance &, const float);
  void add_pin_cost_to_via_cost(const pin &, const float);
  void add_pin_cost_to_via_cost(const padstack &, const instance &, const float);

private:
  BoardGrid mBg;
  kicadPcbDataBase &mDb;

  std::vector<std::string> mGridLayerToName;
  std::unordered_map<std::string, int> mLayerNameToGrid;

  // TODO
  // Temporary value
  const float pinCost = 100000.0;
  // Put below stuff to globalParam:: ??
  double mMinX = std::numeric_limits<double>::max();
  double mMaxX = std::numeric_limits<double>::min();
  double mMinY = std::numeric_limits<double>::max();
  double mMaxY = std::numeric_limits<double>::min();
  // Take const to below?
  const unsigned int inputScale = 10;
  const unsigned int enlargeBoundary = 10;
  const float grid_factor = 0.1; //For outputing
};

#endif