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
//#include "kicadPcbDataBase.h"
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
  GridBasedPlacer() {}
  //dtor
  ~GridBasedPlacer(){}

  void test_placer_flow();
  void HPlace();


private:
  // Utility
  int dbLengthToGridLength(const double dbLength) { return (int)ceil(dbLength * inputScale); }

  void random_initial_placement();
  void set_boundaries();
  void initialize_params(map<int, vector<Pin> > &netToCell);
  void validate_move(Node &node, double rx, double ry);
  double cost(map<int, vector<Pin> > &netToCell,
              int temp_debug = 0);
  double cost_partial(vector < Node *> &nodes, map<int, vector<Pin> > &netToCell);
  double cell_overlap();
  double wirelength(map<int, vector<Pin> > &netToCell);
  double cell_overlap_partial(vector < Node* > &nodes);
  double wirelength_partial(vector < Node* > &nodes, map<int, vector<Pin> > &netToCell);
  double rudy(map<int, vector<Pin> > &netToCell);
  float annealer(map<int, vector<Pin> > &netToCell,
                string initial_pl);
  float multistart();
  double initialize_temperature(double &Temperature, map<int, vector<Pin> > &netToCell);
  void update_temperature(double& Temperature);
  double initiate_move(double current_cost,
                       double & Temperature,
                       map<int, vector<Pin> > &netToCell);
  bool check_move(double prevCost,
                  double newCost,
                  double &Temperature);
  void project_soln();
  void random_placement(int xmin, int xmax, int ymin, int ymax, Node &n);
  void gen_report(map<string, vector<double> > &report,
                  vector< double > &accept_ratio_history,
                  map<int, vector<Pin> > &netToCell);
  void update_accept_history(vector< double > &accept_ratio_history, float &accept_ratio);
  void update_rtree(int idx);
  vector < Node > ::iterator random_node();

  int fact(int n);

private:
  BoardGrid mBg;
  //bgi::rtree<std::pair<box, int>, bgi::quadratic<16>> rtree;

  std::vector<std::string> mGridLayerToName;
  std::unordered_map<std::string, int> mLayerNameToGrid;

  std::pair <double,double> wl_normalization;
  std::pair <double,double> area_normalization;
  std::pair <double,double> routability_normalization;
  map<int, vector<Pin> > *netToCell = nullptr;
  vector < vector < Pin > > *netToCellVec = nullptr;
  vector< int > accept_history;
  double Temperature = 0.0;
  int outer_loop_iter = 101;
  int inner_loop_iter = 20;
  double eps = -1.0;
  double t_0 = 1.0;
  bool var = false;

  bool rotate_flag = 0;
  boost::mt19937 rng;

  double mMinX = std::numeric_limits<double>::max();
  double mMaxX = std::numeric_limits<double>::min();
  double mMinY = std::numeric_limits<double>::max();
  double mMaxY = std::numeric_limits<double>::min();
  // Take const to below?
  const unsigned int inputScale = 10;
  const unsigned int enlargeBoundary = 10;
  const float grid_factor = 0.1; 
  const int debug = 0;
};
