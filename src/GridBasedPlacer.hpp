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
#include "HPlacerUtils.hpp"
//#include "readScl.h"
#include "time.h"
//#include "taskflow/taskflow.hpp"
typedef boost::geometry::model::d2::point_xy<double> Point;

class GridBasedPlacer {
public:
  //ctor
  GridBasedPlacer(kicadPcbDataBase &db) : mDb(db) {}
  GridBasedPlacer() {}
  //dtor
  ~GridBasedPlacer(){}

  /* SWIG Methods Start */
  void test_hplacer_flow();
  kicadPcbDataBase &test_placer_flow();
  kicadPcbDataBase &getDb() { return mDb; }
  kicadPcbDataBase &mDb;

  double get_total_cost() { return cost_hist.back(); }
  double get_wirelength_cost() { return cost_hist.back(); }
  double get_overlap_cost() { return cost_hist.back(); }
  double get_temperature() { return Temperature; }

  void set_overlap_weight(double _cst) {l1 = _cst;}
  void set_wirelength_weight(double _cst) {l1 = 1-_cst;}
  void set_two_sided(bool _2sided) {two_sided = _2sided; shift_proba += layer_change_proba; layer_change_proba -= layer_change_proba; }
  void set_initial_move_radius(double _eps) { vector < Node > ::iterator nodeit = nodeId.begin(); 
                                              for (nodeit = nodeId.begin(); nodeit != nodeId.end(); ++nodeit) { 
                                                nodeit->sigma = _eps; } }
  void set_rtree(bool _rt) { rt = _rt; }
  void set_lam(bool _lam) { lam = _lam; }
  void set_lamtemp_update(double _coef) { lamtemp_update = _coef; }

  void set_num_iterations(int iter) {outer_loop_iter = iter;}
  void set_iterations_moves(int iter) {inner_loop_iter = iter;}
  void set_initial_temperature(double tmp) {t_0 = tmp;}

  /* SWIG Methods End */


  double h_cellDensity();
  void h_initialize_params(map<int, vector<Module *> > &netToCell);
  void h_validate_move(Module *node, double rx, double ry);
  double h_cost(map<int, vector <Module *> > &netToCell, int temp_debug = 0);
  double h_cost_partial(vector < Module *> &nodes, map<int, vector<Module *> > &netToCell);
  double h_cell_overlap();
  double h_wirelength(map<int, vector<Module *> > &netToCell);
  double h_cell_overlap_partial(vector < Module * > &nodes);
  double h_wirelength_partial(vector < Module * > &nodes, map<int, vector<Module *> > &netToCell);
  double h_rudy(map<int, vector<Module *> > &netToCell);
  float h_annealer(map<int, vector<Module *> > &netToCell, string initial_pl,int level=0);
  double h_initialize_temperature(double &Temperature,map<int, vector<Module *> > &netToCell);
  double h_initiate_move(double current_cost, map<int, vector<Module *> > &netToCell);
  void h_random_placement(int xmin, int xmax, int ymin, int ymax, Module &n);
  void h_random_initial_placement();
  bool h_check_move(double prevCost, double newCost, double & Temperature);
  vector < Module * > ::iterator h_random_node();

private:
  void hplace(map<int, vector<pPin> > &netToCell, string initial_pl);

  int dbLengthToGridLength(const double dbLength) { return (int)ceil(dbLength * inputScale); }

  void random_initial_placement();
  void set_boundaries();
  void initialize_params(map<int, vector<pPin> > &netToCell);
  void validate_move(Node &node, double rx, double ry);
  double cost(map<int, vector<pPin> > &netToCell, int temp_debug = 0);
  double cost_partial(vector < Node *> &nodes, map<int, vector<pPin> > &netToCell);
  double cell_overlap();
  double wirelength(map<int, vector<pPin> > &netToCell);
  double cell_overlap_partial(vector < Node* > &nodes);
  double wirelength_partial(vector < Node* > &nodes, map<int, vector<pPin> > &netToCell);
  double rudy(map<int, vector<pPin> > &netToCell);
  float annealer(map<int, vector<pPin> > &netToCell, string initial_pl);
  float hannealer(map<int, vector<Module *> > &netToCell, string initial_pl, int level=0);
  float multistart();
  double initialize_temperature(map<int, vector<pPin> > &netToCell);
  void update_temperature();
  void modified_lam_update(int i);
  double initiate_move(double current_cost, map<int, vector<pPin> > &netToCell);
  bool check_move(double prevCost, double newCost);
  void project_soln();
  void random_placement(int xmin, int xmax, int ymin, int ymax, Node &n);
  void gen_report(map<string, vector<double> > &report, map<int, vector<pPin> > &netToCell);
  //void update_accept_history(vector< double > &accept_ratio_history, float &accept_ratio);
  void update_rtree(int idx);
  vector < Node > ::iterator random_node();

private:

  vector < double > l_hist;
  vector < double > density_hist;
  vector < double > var_hist;
  vector < double > temp_hist;

  bool do_hplace = false;

  std::pair <double,double> density_normalization;
  double densityAlpha = 0.0001;
  int densityFlag = 0;


  // best-so-far solution variables
  vector < Node > bestSol;

  bool hierachical = false;

  double best_wl = 0.0;
  double best_overlap = std::numeric_limits<double>::max();

  // rtree datastructure for fast overlap check
  bool rt = false;
  bgi::rtree<std::pair<boost::geometry::model::box< model::d2::point_xy<double> >, int>, bgi::quadratic<16> > rtree;

  std::vector<std::string> mGridLayerToName;
  std::unordered_map<std::string, int> mLayerNameToGrid;

  // cost histories
  vector < double > cost_hist;
  vector < double > wl_hist;
  vector < double > oa_hist;

  // cost normalization terms
  int initial_loop_iter = 100;
  std::pair <double,double> wl_normalization;
  std::pair <double,double> area_normalization;
  std::pair <double,double> routability_normalization;

  map<int, vector<pPin> > *netToCell = nullptr;
  vector < vector < pPin > > *netToCellVec = nullptr;
  vector< int > accept_history;

  // annealing parameters
  double t_0 = 0.5;
  double Temperature;
  int outer_loop_iter = 101;
  int inner_loop_iter = 20;
  double eps = -1.0;
  bool var = false;
  double l1 = 0.4;
  double shift_var = 1.0;
  double ssamp = 0.0;

  // annealing move parameters
  float rotate_proba = 0.15;
  float layer_change_proba = 0.1;
  float swap_proba = 0.25;
  float shift_proba = 0.5;
  bool rotate_flag = 0;

  // two-sided placement
  bool two_sided = false;

  // modified lam schedule params
  bool lam = true;
  double AcceptRate = 0.5;
  double LamRate = 0.5;
  double lamtemp_update = 0.85;

  // boost mt random number generator
  boost::mt19937 rng;

  // board boundaries
  double mMinX = std::numeric_limits<double>::max();
  double mMaxX = std::numeric_limits<double>::min();
  double mMinY = std::numeric_limits<double>::max();
  double mMaxY = std::numeric_limits<double>::min();

  const unsigned int inputScale = 10;
  const float grid_factor = 0.1; 
  const int debug = 0;
};
