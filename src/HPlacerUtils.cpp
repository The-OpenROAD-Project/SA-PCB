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

#include "HPlacerUtils.hpp"
#define PI 3.14159265

using namespace std::chrono;

BOOST_GEOMETRY_REGISTER_BOOST_TUPLE_CS(cs::cartesian)

namespace bg = boost::geometry;
namespace bnu = boost::numeric::ublas;
namespace bgi = boost::geometry::index;
typedef bg::model::box< bg::model::d2::point_xy<double> > box;
typedef boost::geometry::model::polygon<boost::geometry::model::d2::point_xy<double> > polygon;

/**
random_placement
Randomly place and orient a single component within a bounded region
*/
void GridBasedPlacer::h_random_placement(int xmin, int xmax, int ymin, int ymax, Module &n) {
  boost::uniform_int<> uni_distx(xmin,xmax);
  boost::variate_generator<boost::mt19937&, boost::uniform_int<> > unix(rng, uni_distx);
  int rx = unix();

  boost::uniform_int<> uni_disty(ymin,ymax);
  boost::variate_generator<boost::mt19937&, boost::uniform_int<> > uniy(rng, uni_disty);
  int ry = uniy();

  int ro = 0.0;

  n.setPos(rx, ry);
  this->h_validate_move(&n, rx, ry);
}

/**
initial_placement
Randomly place and orient all movable components in the board area
*/
void GridBasedPlacer::h_random_initial_placement() {
  vector < Module * > ::iterator itNode;
  for (itNode = moduleId.begin(); itNode != moduleId.end(); ++itNode) {
    if (!(*itNode) -> fixed) {
      this->h_random_placement(mMinX, mMaxX - max((*itNode) -> width, (*itNode) -> height), mMinY, mMaxY - max((*itNode) -> width, (*itNode) -> height), *(*itNode));
    }
  }
}

/**
wirelength
Computes HPWL for all nets
*/
double GridBasedPlacer::h_wirelength(map<int, vector<Module *> > &netToCell) {
  // compute HPWL
  //map<int, vector < Pin > > ::iterator itNet;
  //vector < Pin > ::iterator itCellList;
  map<int, vector<Module *> > ::iterator itNet;
  vector < Module * > ::iterator itCellList;
  double xVal, yVal, wireLength = 0;
  for (itNet = netToCell.begin(); itNet != netToCell.end(); ++itNet) {
    double minXW = mMaxX, minYW = mMaxY, maxXW = mMinX, maxYW = mMinY;
    for (itCellList = itNet -> second.begin(); itCellList != itNet -> second.end(); ++itCellList) {
      xVal = moduleId[(*itCellList)->idx-1]->xBy2;
      yVal = moduleId[(*itCellList)->idx-1]->yBy2;

      if (xVal < minXW)
        minXW = xVal;
      if (xVal > maxXW)
        maxXW = xVal;
      if (yVal < minYW)
        minYW = yVal;
      if (yVal > maxYW)
        maxYW = yVal;
    }
    wireLength += (abs((maxXW - minXW)) + abs((maxYW - minYW)));
  }
  return wireLength;
}

/**
cell_overlap
Compute sum squared overlap for all components
*/
double GridBasedPlacer::h_cell_overlap() {
  double overlap = 0.0;

  for(size_t i = 0; i < moduleId.size(); i++) {
    //if(nodeId[i].terminal) { continue; }
    for(size_t j = i; j < moduleId.size(); j++) {
      if (i == j) {continue;}
      //if(nodeId[j].terminal) { continue; }
      if(!intersects(moduleId[i]->poly, moduleId[j]->poly) || (moduleId[i]->fixed && moduleId[j]->fixed)) {
        continue;
      } else {
        double oa = 0.0;
        std::deque<polygon> intersect_poly;
        boost::geometry::intersection(moduleId[i]->poly, moduleId[j]->poly, intersect_poly);

        BOOST_FOREACH(polygon const& p, intersect_poly) {
            oa +=  bg::area(p);
        }
        if(oa < 1.0) {
            oa = 1.0;
        }
        overlap +=  pow(oa,2);
      }
    }
  }
  return overlap;
}

/**
wireLength_partial
Compute HPWL for select nets
*/
double GridBasedPlacer::h_wirelength_partial(vector < Module *> &nodes, map<int, vector<Module *> > &netToCell) {
  //vector < Pin > net;
  vector < Module* > net;
  vector < int > ::iterator itNet;
  //vector < Pin > ::iterator itCellList;
  vector < Module * > ::iterator itCellList;
  vector < Module *>::iterator itNode;
  unordered_set < int > net_history;

  double xVal, yVal, wireLength = 0;
  for (itNode = nodes.begin(); itNode != nodes.end(); ++itNode) {
    for (itNet = (*itNode)->Netlist.begin(); itNet != (*itNode)->Netlist.end(); ++itNet) {
      if (net_history.find(*itNet) == net_history.end()) {
        net_history.insert(*itNet);
      } else {
        continue;
      }
      net = netToCell.find(*itNet)->second;
      double minXW = mMaxX, minYW = mMaxY, maxXW = mMinX, maxYW = mMinY;
      for (itCellList = net.begin(); itCellList != net.end(); ++itCellList) {
        int orient = moduleId[(*itCellList)->idx]->orientation;
        xVal = moduleId[(*itCellList)->idx-1]->xBy2;
        yVal = moduleId[(*itCellList)->idx-1]->yBy2;

        if (xVal < minXW)
          minXW = xVal;
        if (xVal > maxXW)
          maxXW = xVal;
        if (yVal < minYW)
          minYW = yVal;
        if (yVal > maxYW)
          maxYW = yVal;
      }
      wireLength += (abs((maxXW - minXW)) + abs((maxYW - minYW)));
    }
  }
  return wireLength;
}

/**
rudy
Computes a routability score
*/
double GridBasedPlacer::h_cellDensity() {
  static bnu::matrix<double> D (static_cast<int>(abs(mMaxY)+abs(mMinY) + 1), static_cast<int>(abs(mMaxX)+abs(mMinX) + 1), 0.0);
  static bnu::matrix<double> D_sup (static_cast<int>(abs(mMaxY)+abs(mMinY) + 1), static_cast<int>(abs(mMaxX)+abs(mMinX) + 1), 1.0);

  D.clear();

  for(size_t i = 0; i < moduleId.size(); i++) {

    int minXW = floor(moduleId[i]->xCoordinate);
    int minYW = floor(moduleId[i]->yCoordinate);
    int maxXW = ceil(moduleId[i]->xCoordinate + moduleId[i]->width);
    int maxYW = ceil(moduleId[i]->yCoordinate + moduleId[i]->height);

    //double minOverflowX = abs(minXW - moduleId[i]->xCoordinate);
    //double maxOverflowX = abs();
    //double minOverflowY = abs();
    //double maxOverflowY = abs();

    for (unsigned i = max(minYW, 0); i < maxYW; ++ i) {
        for (unsigned j = max(minXW, 0); j < maxXW; ++ j) {
          //D (i,j) += 1;
          double xover = 1.0;
          double yover = 1.0;
          if (i < moduleId[i]->yCoordinate) {
            yover += moduleId[i]->yCoordinate - floor(moduleId[i]->yCoordinate);
          } else if (i < moduleId[i]->yCoordinate + moduleId[i]->height) {
            yover += ceil(moduleId[i]->yCoordinate + moduleId[i]->height) - moduleId[i]->yCoordinate + moduleId[i]->height;
          } else {
            yover += 1.0;
          }

          if (i < moduleId[j]->xCoordinate) {
            xover += moduleId[i]->xCoordinate - floor(moduleId[i]->xCoordinate);
          } else if (j < moduleId[i]->xCoordinate + moduleId[i]->width) {
            xover += ceil(moduleId[i]->xCoordinate + moduleId[i]->width) - moduleId[i]->xCoordinate + moduleId[i]->width;
          } else {
            xover += 1.0;
          }
          D(i,j) += xover * yover;
        }
      }
  }

  double r = 0.0;
  for (unsigned i = 0; i < D.size1 (); ++ i) {
      for (unsigned j = 0; j < D.size2 (); ++ j) {
        D(i,j) = exp(densityAlpha * (D(i,j) - D_sup(i,j))/D_sup(i,j));
        r += D(i,j);
      }
  }
  return r;
}

/**
cell_overlap_partial
Compute sum squared overlap for select components
*/
double GridBasedPlacer::h_cell_overlap_partial(vector < Module *> &nodes) {
  double overlap = 0.0;
  unordered_set < int > cell_history;
  for(size_t i = 0; i < nodes.size(); i++) {
    //if(nodes[i].terminal) { continue; }
    cell_history.insert(nodes[i]->idx);
    for(size_t j = 0; j < moduleId.size(); j++) {
      //if(nodeId[j].terminal) { continue; }
      /*
      for ( Rtree::const_query_iterator it = rtree.qbegin(index::intersects(nodeId[i].envelope)) ;
        it != rtree.qend() ; ++it ) {
      }*/
      if (cell_history.find(moduleId[j]->idx) != cell_history.end()) {
        continue;
      }
      if(!intersects(nodes[i]->poly, moduleId[j]->poly) || (nodes[i]->fixed && moduleId[j]->fixed)) {
        continue;
      } else {
        double oa = 0.0;
        std::deque<polygon> intersect_poly;
        boost::geometry::intersection(nodes[i]->poly, moduleId[j]->poly, intersect_poly);

        BOOST_FOREACH(polygon const& p, intersect_poly) {
            oa +=  bg::area(p);
        }
        if(oa < 1.0) {
            oa = 1.0;
        }
        overlap +=  pow(oa,2);
      }
    }
  }
  return overlap;
}

/**
rudy
Computes a routability score
*/
double GridBasedPlacer::h_rudy(map<int, vector<Module *> > &netToCell) {
  static bnu::matrix<double> D (static_cast<int>(abs(mMaxY)+abs(mMinY) + 1), static_cast<int>(abs(mMaxX)+abs(mMinX) + 1), 0.0);
  static bnu::matrix<double> D_route_sup (static_cast<int>(abs(mMaxY)+abs(mMinY) + 1), static_cast<int>(abs(mMaxX)+abs(mMinX) + 1), 1.0);

  D.clear();

  //map<int, vector < Pin > > ::iterator itNet;
  map<int, vector < Module* > > ::iterator itNet;
  //vector < Pin > ::iterator itCellList;
  vector < Module* > ::iterator itCellList;

  for (itNet = netToCell.begin(); itNet != netToCell.end(); ++itNet) {
    double xVal, yVal, hpwl = 0.0;
    double minXW = mMaxX, minYW = mMaxY, maxXW = mMinX, maxYW = mMinY;
    double rudy = 0.0;
    for (itCellList = itNet -> second.begin(); itCellList != itNet -> second.end(); ++itCellList) {
      int orient = nodeId[(*itCellList)->idx].orientation;
 
      xVal = nodeId[(*itCellList)->idx-1].xBy2;
      yVal = nodeId[(*itCellList)->idx-1].yBy2;

      // compute pin position from orientation & offsets
      if(orient == 0) { // 0
        xVal = xVal + (*itCellList)->x_offset;
        yVal = yVal + (*itCellList)->y_offset;
      } else if(orient == 2) { // 90
        xVal = xVal + (*itCellList)->y_offset;
        yVal = yVal - (*itCellList)->x_offset;
      } else if(orient == 4) { // 180
        xVal = xVal - (*itCellList)->x_offset;
        yVal = yVal - (*itCellList)->y_offset;
      } else if(orient == 6) { // 270
        xVal = xVal - (*itCellList)->y_offset;
        yVal = yVal + (*itCellList)->x_offset;
      } else {
        double rad = (orient*45.0*PI/180.0);
        xVal = (*itCellList)->y_offset*sin(rad) + (*itCellList)->x_offset*cos(rad) + xVal;
        yVal = (*itCellList)->y_offset*cos(rad) - (*itCellList)->x_offset*sin(rad) + yVal;
      }

      if (xVal < minXW)
        minXW = xVal;
      if (xVal > maxXW)
        maxXW = xVal;
      if (yVal < minYW)
        minYW = yVal;
      if (yVal > maxYW)
        maxYW = yVal;
      // set read_net
      
      //for (unsigned i = xVal-5; i < xVal+5; ++ i) {
      //    for (unsigned j = yVal-5; j < yVal+5; ++ j) {
      //      D (i,j) += 5;
      //    }
      //}
    }
    minXW = max(minXW, mMinX);
    minYW = max(minYW, mMinY);
    maxXW = min(maxXW, mMaxX);
    maxYW = min(maxYW, mMaxY);
    // now have boundary of net
    hpwl = (abs((maxXW - minXW)) + abs((maxYW - minYW)));
    rudy = hpwl / (max((maxXW - minXW)*(maxYW - minYW), 1.0)); // rudy density
    // set read_net

    for (unsigned i = max(minYW,0.0); i < maxYW; ++ i) {
        for (unsigned j = max(minXW,0.0); j < maxXW; ++ j) {
          D (i,j) += rudy;
        }
      }
  }
  double r = 0.0;
  for (unsigned i = 0; i < D.size1 (); ++ i) {
      for (unsigned j = 0; j < D.size2 (); ++ j) {
        D(i,j) = exp((D(i,j) - D_route_sup(i,j))/D_route_sup(i,j));
        r += D(i,j);
      }
  }


  return r;
}

double GridBasedPlacer::h_cost(map<int, vector<Module *> > &netToCell, int temp_debug) {
  double l2 = 1 - l1;
  double wirelength_cost = l1*0.9*(this->h_wirelength(netToCell) - wl_normalization.first)/(wl_normalization.second - wl_normalization.first);
  double overlap_cost = 0.0;
  if (densityFlag == 0) {
    overlap_cost = l2  * (this->h_cell_overlap() - area_normalization.first)/(area_normalization.second - area_normalization.first);
  } else {
    overlap_cost = l2 * this->cellDensity();
  }
  //double routability_cost = l1 * 0.1 * (this->rudy(netToCell) - routability_normalization.first)/(routability_normalization.second - routability_normalization.first);
  double total_cost = wirelength_cost + overlap_cost;// + routability_cost;
  cout << "wirelength: " << this->h_wirelength(netToCell) << endl;
  cout << "overlap: " << this->h_cell_overlap() << endl;
  cout << "l1: " << l1 << " l2: " << l2 << endl;
  cout << "wirelength_cost: " << wirelength_cost  << endl;
  cout << "overlap_cost: " << overlap_cost  << endl;
  //cout << "routability_cost: " << routability_cost  << endl;
  cout << "cost: " << total_cost << endl;
  return total_cost;
}

double GridBasedPlacer::h_cost_partial(vector < Module *> &nodes, map<int, vector<Module *> > &netToCell) {
  double l2 = 1-l1;
  if(densityFlag == 0) {
    return l1 * 0.9 * (this->h_wirelength_partial(nodes, netToCell) - wl_normalization.first)/(wl_normalization.second - wl_normalization.first) +
           l2 * (this->h_cell_overlap_partial(nodes) - area_normalization.first)/(area_normalization.second - area_normalization.first); //+
           //l1 * 0.1 * (this->rudy(netToCell) - routability_normalization.first)/(routability_normalization.second - routability_normalization.first);
  } else {
    return l1 * 0.9 * (this->h_wirelength_partial(nodes, netToCell) - wl_normalization.first)/(wl_normalization.second - wl_normalization.first) +
           l2 * (this->cellDensity()); //+
           //l1 * 0.1 * (this->rudy(netToCell) - routability_normalization.first)/(routability_normalization.second - routability_normalization.first);
  }
}

/**
random_node
Select random node from set of nodes
*/
vector < Module * >::iterator GridBasedPlacer::h_random_node() {
  vector < Module * > ::iterator itNode = moduleId.begin();
  int size = moduleId.size();
  boost::uniform_int<> uni_dist(0,size-1);
  boost::variate_generator<boost::mt19937&, boost::uniform_int<> > uni(rng, uni_dist);
  int randint = uni();
  std::advance(itNode, randint);
  return itNode;
}

/**
validate_move
validates a shift, project within board boundary
*/
void GridBasedPlacer::h_validate_move(Module *node, double rx, double ry) {
  rx = max(rx, mMinX);
  ry = max(ry, mMinY);

  rx = min(rx, mMaxX - node->width);
  ry = min(ry, mMaxY - node->height);

  node->setPos(rx,ry);
}

double c = 0.0;
double GridBasedPlacer::h_initiate_move(double current_cost, double & Temperature, map<int, vector<Module *> > &netToCell) {
  // Initate a transition
  int state = -1;
  double prevCost = 0.0;
  if(debug > 1) {
    cout << "current_cost: " << current_cost << endl;
  }

  vector < Module * > ::iterator rand_node1;
  vector < Module * > ::iterator rand_node2;

  vector < Module * > perturbed_nodes;

  int r = 0;
  rand_node1 = this->h_random_node();
  while((*rand_node1)->terminal || (*rand_node1)->fixed) {
    rand_node1 = this->h_random_node();
  }
  perturbed_nodes.push_back(*rand_node1);

  //rtree.remove(std::make_pair(rand_node1->envelope, rand_node1->idx));
  double rand_node1_orig_x = (*rand_node1)->xCoordinate;
  double rand_node1_orig_y = (*rand_node1)->yCoordinate;
  double rand_node2_orig_x = 0.0;
  double rand_node2_orig_y = 0.0;

  if(debug > 1) {
    cout << "=======" << endl;
    cout <<  "name: " << (*rand_node1)->name << endl;
  }

  boost::uniform_real<> uni_dist(0,1);
  boost::variate_generator<boost::mt19937&, boost::uniform_real<> > uni(rng, uni_dist);
  double i = uni();

  if (i< p_swap) { // swap
    state = 0;
    if(debug > 1) {
      cout << "swap" << endl;
    }
    rand_node2 = this->h_random_node();
    while((*rand_node2)->terminal || (*rand_node2)->idx == (*rand_node1)->idx || (*rand_node2)->fixed) {    
      rand_node2 = this->h_random_node();
    }
    perturbed_nodes.push_back(*rand_node2);
    prevCost = this->h_cost_partial(perturbed_nodes,netToCell);

    //rtree.remove(std::make_pair(rand_node2->envelope, rand_node2->idx));
    rand_node2_orig_x = (*rand_node2)->xCoordinate;
    rand_node2_orig_y = (*rand_node2)->yCoordinate;

    h_validate_move(*rand_node1, rand_node2_orig_x, rand_node2_orig_y);
    h_validate_move(*rand_node2, rand_node1_orig_x, rand_node1_orig_y);
    //rtree.insert(std::make_pair(rand_node1->envelope, rand_node1->idx));
    //rtree.insert(std::make_pair(rand_node2->envelope, rand_node2->idx));
  } else if (i < p_shift) { // shift
    state = 1;
    if(debug > 1) {
      cout << "shift" << endl;
    }
    prevCost = this->h_cost_partial(perturbed_nodes, netToCell);

    double sigma = (*rand_node1)->sigma;

    boost::normal_distribution<> nd(0.0, sigma);
    boost::variate_generator<boost::mt19937&,
                             boost::normal_distribution<> > var_nor(rng, nd);

    double dx = var_nor();
    double dy = var_nor();

    double rx = rand_node1_orig_x + dx;
    double ry = rand_node1_orig_y + dy;

    h_validate_move(*rand_node1, rx, ry);
    //rtree.insert(std::make_pair(rand_node1->envelope, rand_node1->idx));
  } else if (i < p_rotate) { // rotate
    state = 2;
    if(debug > 1) {
      cout << "rotate" << endl;
    }
    prevCost = this->h_cost_partial(perturbed_nodes, netToCell);

    if (rotate_flag == 0) {
      boost::uniform_int<> uni_dist(0,3);
      boost::variate_generator<boost::mt19937&, boost::uniform_int<> > uni(rng, uni_dist);
      r = uni();
      r = r*2;
    } else {
      boost::uniform_int<> uni_dist(0,7);
      boost::variate_generator<boost::mt19937&, boost::uniform_int<> > uni(rng, uni_dist);
      r = uni();
    }
    //(*rand_node1)->setRotation(r);
    h_validate_move(*rand_node1, rand_node1_orig_x, rand_node1_orig_y);

    //rtree.insert(std::make_pair(rand_node1->envelope, rand_node1->idx));
  }
  double transition_cost = this->h_cost_partial(perturbed_nodes,netToCell);
  double updated_cost = current_cost - prevCost + transition_cost;

  bool accept = this->check_move(current_cost, updated_cost, Temperature);
  if (!accept) {
    if(debug > 1) {
      cout << "reject" << endl;
    }
    // revert state
    //rtree.remove(std::make_pair(rand_node1->envelope, rand_node1->idx));
    if (state == 0) {
      //rtree.remove(std::make_pair(rand_node2->envelope, rand_node2->idx));
      (*rand_node1)->setPos(rand_node1_orig_x,rand_node1_orig_y);
      (*rand_node2)->setPos(rand_node2_orig_x,rand_node2_orig_y);
      //rtree.insert(std::make_pair(rand_node2->envelope, rand_node2->idx));
    } else if (state == 1) {
      (*rand_node1)->setPos(rand_node1_orig_x,rand_node1_orig_y);
    } else if (state == 2) {
      //(*rand_node1)->setRotation(8-r);
      (*rand_node1)->setPos(rand_node1_orig_x,rand_node1_orig_y);
    }
    //rtree.insert(std::make_pair(rand_node1->envelope, rand_node1->idx));
    return current_cost;
  }
  else {
    if(debug > 1) {
      cout << "accept" << endl;
    }
    return updated_cost;
  }
}

/**
check_move
either accept or reject the move based on current & previous temperature & cost
*/
bool GridBasedPlacer::h_check_move(double prevCost, double newCost, double & Temperature) {
  double delCost = 0;
  boost::uniform_real<> uni_dist(0,1);
  boost::variate_generator<boost::mt19937&, boost::uniform_real<> > uni(rng, uni_dist);
  double prob = uni();
  if(debug > 1) {
    cout << "new cost: " << newCost << endl;
  }
  delCost = newCost - prevCost;
  if (delCost <= 0 || prob <= (exp(-delCost/Temperature))) {
    prevCost = newCost;
    accept_history.push_back(1);
    return true;
  } else {
    accept_history.push_back(0);
    return false;
  }
}

double GridBasedPlacer::h_initialize_temperature(double &Temperature, map<int, vector<Module *> > &netToCell) {
  double t = 0.0;
  double emax = 0.0;
  double emin = 0.0;
  double xt = 1.0;
  double x0 = 0.84;
  double p = 2.0;
  this->random_initial_placement();
  for(int i=1; i<=10; i++){
    for(int j=1; j<=10; j++){
      this->random_initial_placement();
      emax += exp(h_cost(netToCell)/t);
      this->h_initiate_move(0.0, Temperature, netToCell);
      emin += exp(h_cost(netToCell)/t);
    }
    xt = emax/emin;
    t = t * pow(log(xt),1/p)/log(x0);
  }
  return t;
}

/**
annealer
main loop for sa algorithm
*/
float GridBasedPlacer::h_annealer(map<int, vector<Module *> > &netToCell, string initial_pl, int level) {
  double Temperature = t_0;
  int num_components = 0;

  vector < Module * > ::iterator itNode;
  map < string, vector < double > > report;


  vector< int > accept_history;
  float accept_ratio = 0.0;
  vector< double > accept_ratio_history;

  cout << "calculating initial params..." << endl;
  this->h_initialize_params(netToCell);

  cout << "calculating initial cost..." << endl;  
  double cst = this->h_cost(netToCell,-1);
  
  if(var) {
    Temperature = this->h_initialize_temperature(Temperature, netToCell);
  }
  int idx = 0;

  for (itNode = moduleId.begin(); itNode != moduleId.end(); ++itNode) {
    //rtree.insert(std::make_pair(itNode -> envelope, idx));
    idx+=1;
    if(!(*itNode) -> terminal && !(*itNode) -> fixed) {
      num_components += 1;
    }
  }

  high_resolution_clock::time_point t1 = high_resolution_clock::now();
  long long int ii = 0; // outer loop iterator
  int i = 0; // inner loop iterator
  cout << "beginning iterations..." << endl;
  while (ii < outer_loop_iter) {
    i = inner_loop_iter*num_components; 
    if(ii > 1 && ii % 5 == 0) {
      high_resolution_clock::time_point t2 = high_resolution_clock::now();
      duration<double> time_span = duration_cast< duration<double> >(t2 - t1);

      cout << "******" << ii << "******" << endl;
      cout << "iteration: " << ii << endl;
      cout << "time: " <<  time_span.count() << " (s)" << endl;
      cout << "move/time: " <<  i*ii/time_span.count() << endl;
      cout << "time remaining: " <<  time_span.count()/ii * (outer_loop_iter-ii) << " (s)" << endl;
      cout << "temperature: " << Temperature << endl;
      cout << "acceptance ratio: " << accept_ratio << endl;
      wl_hist.push_back((h_wirelength(netToCell) - wl_normalization.first)/(wl_normalization.second - wl_normalization.first));
      if (densityFlag == 0) {
        oa_hist.push_back((h_cell_overlap() - area_normalization.first)/(area_normalization.second - area_normalization.first));
      } else {
        oa_hist.push_back(cellDensity());
      }
      l_hist.push_back(l1);
      temp_hist.push_back(Temperature);
      this->h_cost(netToCell,-1);

      nodeId = H.update_cell_positions_at_level(nodeId, level);
      writePlFile("./cache/"+std::to_string( level )+"_"+std::to_string( ii )+".pl");

    }
    while (i > 0) {
      cst = this->h_initiate_move(cst, Temperature, netToCell);
      this->update_accept_history(accept_ratio_history, accept_ratio);
      report["cost_hist"].push_back(cst);
      i -= 1;
    }

    this->update_temperature(Temperature);
    ii += 1;
  }
  return this->h_cost(netToCell);
}