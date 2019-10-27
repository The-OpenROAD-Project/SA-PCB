///////////////////////////////////////////////////////////////////////////////
// Authors: Chester Holtz, Devon Merrill, James (Ting-Chou) Lin, Connie (Yen-Yi) Wu
//          (respective Ph.D. advisors: Chung-Kuan Cheng, Andrew B. Kahng, Steven Swanson).
//
// BSD 3-Clause License
//
// Copyright (c) 2018, The Regents of the University of California
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
//   list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holder nor the names of its
//   contributors may be used to endorse or promote products derived from
//   this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
///////////////////////////////////////////////////////////////////////////////

#define BOOST_NO_AUTO_PTR

#include <cstdio>
#include "BoardGrid.h"
#include "kicadPcbDataBase.h"
#include "globalParam.h"
#include "util.h"

#include <functional>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctype.h>
#include <map>
#include <vector>
#include <unordered_set>
#include <algorithm>
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

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "readFiles.hpp"
//#include "readScl.h"
#include "time.h"
//#include "taskflow/taskflow.hpp"

class GridBasedPlacer
{
public:
  //ctor
  GridBasedPlacer(kicadPcbDataBase &db) : mDb(db) {}
  GridBasedPlacer() {}
  //dtor
  ~GridBasedPlacer(){}

  kicadPcbDataBase &test_placer_flow();

  struct boundaries {
    double minX, maxX, minY, maxY = 0.0;
  };
  kicadPcbDataBase &getDb() { return mDb; }
private:
  // Utility
  int dbLengthToGridLength(const double dbLength) { return (int)ceil(dbLength * inputScale); }
 /*
  bool dbPointToGridPoint(const point_2d &dbPt, point_2d &gridPt);
  bool gridPointToDbPoint(const point_2d &gridPt, point_2d &dbPt);
  void addPinCost(const pin &, const float);
  void addPinCost(const padstack &, const instance &, const float);
  void add_pin_cost_to_via_cost(const pin &, const float);
  void add_pin_cost_to_via_cost(const padstack &, const instance &, const float);
*/
  kicadPcbDataBase &mDb;
  void random_initial_placement(int rotate_flag);
  void set_boundaries();
  void initialize_params(std::pair <double,double> &wl_normalization,
                         std::pair <double,double> &area_normalization,
                         std::pair <double,double> &routability_normalization,
                         map<int, vector<Pin> > &netToCell,
                         int rotate_flag);
  void validate_move(Node &node, double rx, double ry);
  double cost(
              std::pair <double,double> &wl_normalization,
              std::pair <double,double> &area_normalization,
              std::pair <double,double> &routability_normalization,
              map<int, vector<Pin> > &netToCell,
              int temp_debug = 0);
  double cost_partial(vector < Node *> &nodes,
                      std::pair <double,double> &wl_normalization,
                      std::pair <double,double> &area_normalization,
                      std::pair <double,double> &routability_normalization,
                      map<int, vector<Pin> > &netToCell);
  double cell_overlap();
  double wirelength(map<int, vector<Pin> > &netToCell);
  double cell_overlap_partial(vector < Node* > &nodes);
  double wirelength_partial(vector < Node* > &nodes, map<int, vector<Pin> > &netToCell);
  double rudy(map<int, vector<Pin> > &netToCell);
  float annealer(int outer_loop_iter,
                            int inner_loop_iter,
                            double eps,
                            double t_0,bool var,
                            map<int, vector<Pin> > &netToCell,
                            string initial_pl,
                            int rotate_flag);
  float multistart();
  double initialize_temperature(vector< int > &accept_history,
                                double &Temperature,
                                std::pair <double,double> &wl_normalization,
                                std::pair <double,double> &area_normalization,
                                std::pair <double,double> &routability_normalization,
                                map<int, vector<Pin> > &netToCell,
                                int rotate_flag);
  void update_temperature(double& Temperature);
  double initiate_move(double current_cost,
                       vector< int > &accept_history,
                       double & Temperature,
                       std::pair <double,double> &wl_normalization,
                       std::pair <double,double> &area_normalization,
                       std::pair <double,double> &routability_normalization,
                       map<int, vector<Pin> > &netToCell,
                       int rotate_flag);
  bool check_move(double prevCost,
                  double newCost,
                  vector< int > &accept_history,
                  double &Temperature);
  void project_soln();
  void random_placement(int xmin, int xmax, int ymin, int ymax, Node &n, int rotate_flag);
  void gen_report(map<string, vector<double> > &report,
                  vector< double > &accept_ratio_history,
                  std::pair <double,double> &wl_normalization,
                  std::pair <double,double> &area_normalization,
                  std::pair <double,double> &routability_normalization,
                  map<int, vector<Pin> > &netToCell);
  void update_accept_history(vector< int > &accept_history, vector< double > &accept_ratio_history, float &accept_ratio);
  void update_rtree(int idx);
  vector < Node > ::iterator random_node();

  int fact(int n);

private:
  BoardGrid mBg;

  std::vector<std::string> mGridLayerToName;
  std::unordered_map<std::string, int> mLayerNameToGrid;

  //bnu::matrix<double> D (static_cast<int>(abs(b.maxY)+abs(b.minY) + 1), static_cast<int>(abs(b.maxX)+abs(b.minX) + 1), 0.0);

  std::pair <double,double> *wl_normalization = nullptr;
  std::pair <double,double> *area_normalization = nullptr;
  std::pair <double,double> *routability_normalization = nullptr;
  map<int, vector<Pin> > *netToCell = nullptr;
  vector < vector < Pin > > *netToCellVec = nullptr;
  vector< int > *accept_history = nullptr;
  double *Temperature = nullptr;
  bool rotate_flag = false;
  const int debug = 0;

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
