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

#include "GridBasedPlacer.hpp"
#define PI 3.14159265

using namespace std::chrono;

BOOST_GEOMETRY_REGISTER_BOOST_TUPLE_CS(cs::cartesian)

namespace bg = boost::geometry;
namespace bnu = boost::numeric::ublas;
namespace bgi = boost::geometry::index;
typedef bg::model::box< bg::model::d2::point_xy<double> > box2d;
typedef boost::geometry::model::polygon<boost::geometry::model::d2::point_xy<double> > polygon;
typedef bgi::rtree<std::pair<boost::geometry::model::box< model::d2::point_xy<double> >, int>, bgi::quadratic<16> > rtree_t;

vector < Node > nodeId;
vector < Node > bestSol;
double best_wl = 0.0;
map < string, int > name2id;

kicadPcbDataBase &GridBasedPlacer::test_placer_flow() {
    srand(time(NULL));
    int opt;
    string out_file;

    long long int idx = -1;

    string parg = "";
    string initial_pl = "";
    int bound_inc = 0;

    std::cout << "=================test_placer==================" << std::endl;

    map<int, vector<pPin> > netToCell;

    std::vector<instance> &instances =  mDb.getInstances();
    std::vector<net> &nets = mDb.getNets();
    for (int i = 0; i < instances.size(); ++i) {
      Node n;
      nodeId.push_back(n);
    }

    int fixed = 0;
    for (auto &inst : instances) {
        point_2d bbox;
        mDb.getCompBBox(inst.getComponentId(), &bbox); 

        Node n;
        n.setParameterNodes(inst.getName(), bbox.m_x, bbox.m_y, 0, inst.getId());
        nodeId[inst.getId()] = n;
        name2id.insert(pair < string, int > (inst.getName(), inst.getId()));

        double angle = inst.getAngle();
        string ang = "";
        if (angle == 0) {
            ang = "N";
        } else if (angle == 90) {
            ang = "E";
        } else if (angle == 180) {
            ang = "S";
        } else if (angle == 270) {
            ang = "W";
        }
       if(inst.isLocked()) {
           fixed = 1;
       } 
      
        nodeId[name2id[inst.getName()]].setParameterPl(inst.getX() - bbox.m_x/2, inst.getY() - bbox.m_y/2, ang, fixed);
        nodeId[name2id[inst.getName()]].printParameter();
    }

    cout << "calculating boundaries..." << endl;
    points_2d b = mDb.getBoardBoundary();
    mMinX = b[0].m_x;
    mMaxX = b[1].m_x;
    mMinY = b[0].m_y;
    mMaxY = b[1].m_y;

    for (auto &net : nets){

      vector < pPin > pinTemp;
      for (auto &pin : net.getPins())
      {
          pPin p;
          auto &inst = mDb.getInstance(pin.getInstId());
          nodeId[name2id[inst.getName()]].setNetList(net.getId());

          auto &comp = mDb.getComponent(pin.getCompId());
          point_2d pos;
          mDb.getPinPosition(pin, &pos);

          auto &pad = comp.getPadstack(pin.getPadstackId());
          p.set_params(inst.getName(), pos.m_x - nodeId[name2id[inst.getName()]].xBy2, pos.m_y - nodeId[name2id[inst.getName()]].yBy2, nodeId[name2id[inst.getName()]].idx);
          pinTemp.push_back(p);
      }
      netToCell.insert(pair < int, vector< pPin > > (net.getId(), pinTemp));
    }

  cout << "annealing" << endl;
  float cost = this->annealer(netToCell, initial_pl);
  nodeId = bestSol; 
  // write back to db
  cout << "writing back to db..." << endl;
  for (auto &inst : instances) {
      point_2d bbox;
      mDb.getCompBBox(inst.getComponentId(), &bbox);

      int angle = nodeId[inst.getId()].orientation;
      double ang = 0;
      if (angle == 0) {
          ang = 0;
      } else if (angle == 2) {
          ang = 90;
      } else if (angle == 4) {
          ang = 180;
      } else if (angle == 6) {
          ang = 270;
      }
      inst.setAngle(ang);
      inst.setX(nodeId[inst.getId()].xBy2);
      inst.setY(nodeId[inst.getId()].yBy2);
  }
  return mDb;
}


/*
initialize_params
Empirically finds normalization parameters for scaling cost terms.
Currently, we scale by 1/(f) where f is the cost of an expected placement.
*/
void GridBasedPlacer::initialize_params(map<int, vector<pPin> > &netToCell) {

  vector < std::pair <double,double> > normalization_terms;

  int num_components = 0;
  vector < Node > ::iterator itNode;
  for (itNode = nodeId.begin(); itNode != nodeId.end(); ++itNode) {
    if(!itNode -> terminal) {
      num_components += 1;
    }
  }
  std::pair < double, double > wl;
  std::pair < double, double > area;
  std::pair < double, double > rn;

  double sum_wl = 0.0;
  double sum_oa = 0.0;
  double sum_rn = 0.0;
  for (int i = 0; i<100; i++) {
    this->random_initial_placement();
    sum_wl += this->wirelength(netToCell);
    sum_oa += this->cell_overlap();
    sum_rn += this->rudy(netToCell);
  }

  wl.first = 0.0;
  wl.second = sum_wl / 100.0;

  area.first = 0.0;
  area.second = sum_oa / 100.0;
  area.second = max(area.second, 1.0); 

  rn.first = 0.0;
  rn.second = sum_rn / 100.0;

  normalization_terms.push_back(wl);
  normalization_terms.push_back(area);
  normalization_terms.push_back(rn);

  wl_normalization = normalization_terms[0];
  area_normalization = normalization_terms[1];
  routability_normalization = normalization_terms[2];
}

/*
set_boundaries
Calculate board boundaries by finding the maximum-area
rectangular envelope encompasing terminals & fixed modules.
*/
void GridBasedPlacer::set_boundaries() {
  double xval, yval, width, height;
  vector < Node > ::iterator itNode;
  double total_area = 0.0;
  for (itNode = nodeId.begin(); itNode != nodeId.end(); ++itNode) {
    width = itNode -> width;
    height =itNode -> height;
    total_area += width * height;
    if (itNode -> terminal || itNode->fixed) {
      xval = itNode -> xCoordinate;
      yval = itNode -> yCoordinate;
      if (xval < mMinX) {
        mMinX = xval - 2;
      }
      if (xval > mMaxX) {
        mMaxX = xval + width + 2;
      }
      if (yval < mMinY) {
        mMinY = yval - 2;
      }
      if (yval > mMaxY) {
        mMaxY = yval + height + 2;
      }
    }
  }

  if (mMinX == 0 && mMaxX == 0 && mMinY == 0 && mMaxY == 0) {
    cout << "no boundary found" << endl;
    mMaxX = sqrt(1.8*total_area);
    mMaxY = sqrt(1.8*total_area);
  }
  cout << "min: " << mMinX << " " << mMinY << endl;
  cout << "max: " << mMaxX << " " << mMaxY << endl;
}

/*
random_placement
Randomly place and orient a single component within a bounded region
*/
void GridBasedPlacer::random_placement(int xmin, int xmax, int ymin, int ymax, Node &n) {
  boost::uniform_int<> uni_distx(xmin,xmax);
  boost::variate_generator<boost::mt19937&, boost::uniform_int<> > unix(rng, uni_distx);
  int rx = unix();

  boost::uniform_int<> uni_disty(ymin,ymax);
  boost::variate_generator<boost::mt19937&, boost::uniform_int<> > uniy(rng, uni_disty);
  int ry = uniy();

  int ro = 0.0;
  if (rotate_flag == 0) {
    boost::uniform_int<> uni_disto(0,3);
    boost::variate_generator<boost::mt19937&, boost::uniform_int<> > unio(rng, uni_disto);
    ro = unio();
    ro = ro * 2;
  } else {
    boost::uniform_int<> uni_disto(0,7);
    boost::variate_generator<boost::mt19937&, boost::uniform_int<> > unio(rng, uni_disto);
    ro = unio();
  }

  string ostr = n.orient2str(ro);
  n.setParameterPl(rx, ry, ostr, n.fixed);
  this->validate_move(n, rx, ry);
}

/*
initial_placement
Randomly place and orient all movable components in the board area
*/
void GridBasedPlacer::random_initial_placement() {
  vector < Node > ::iterator itNode;
  for (itNode = nodeId.begin(); itNode != nodeId.end(); ++itNode) {
    if (!itNode -> fixed) {
      this->random_placement(mMinX, mMaxX - max(itNode -> width, itNode -> height), mMinY, mMaxY - max(itNode -> width, itNode -> height), *itNode);
    }
  }
}

/*
project_soln
projects the final solution to enforce 0 overlap by iteratively
shifting overlapped components until cumulative overlap reduces below an eps.
*/
void GridBasedPlacer::project_soln() {
  double oa = cell_overlap();
  double eps = 10.0e-5;

  vector < Node > ::iterator nodeit1;
  vector < Node > ::iterator nodeit2;

  box2d env1;
  box2d env2;

  while (oa > eps) {
    for (nodeit1 = nodeId.begin(); nodeit1 != nodeId.end(); ++nodeit1) {
      if (!nodeit1->terminal) {
        env1 = nodeit1->envelope;
        for (nodeit2 = nodeit1++; nodeit2 != nodeId.end(); ++nodeit2) {
          if(nodeit2->idx == nodeit1->idx){continue;}
          if(!intersects(nodeit1->poly, nodeit2->poly) || (nodeit1->fixed && nodeit2->fixed)) {
            continue;
          } else {
            // Project out of envelope
            env2 = nodeit2->envelope;

            double dx = 0.0;
            double dy = 0.0;
            double rx = nodeit1->xCoordinate;
            double ry = nodeit1->yCoordinate;
            double env1minx = env1.min_corner().x();
            double env1miny = env1.min_corner().y();
            double env1maxx = env1.max_corner().x();
            double env1maxy = env1.max_corner().y();

            double env2minx = env2.min_corner().x();
            double env2miny = env2.min_corner().y();
            double env2maxx = env2.max_corner().x();
            double env2maxy = env2.max_corner().y();

            // will get stuck in boundary edge case
            if (abs(env1maxx - env2minx) < abs(env1minx - env2maxx)) {
              // shift left
              dx = env2minx - env1maxx;
            } else {
              // shift right
              dx = env2maxx - env1minx;
            }
            if (abs(env1maxy - env2miny) < abs(env1miny - env2maxy)) {
              // shift down
              dy = env2miny - env1maxy;
            } else {
              // shift up
              dy = env2maxy - env1miny;
            }
            if (dx < dy) {
              // project in x direction
              this->validate_move(*nodeit1,rx + dx,ry);
            } else {
              // project in y direction
              this->validate_move(*nodeit1,rx,ry + dy);
            }
          }
        }
      }
    }
  }
}

/*
wirelength
Computes HPWL for all nets
*/
double GridBasedPlacer::wirelength(map<int, vector<pPin> > &netToCell) {
  // compute HPWL
  map<int, vector < pPin > > ::iterator itNet;
  vector < pPin > ::iterator itCellList;
  double xVal, yVal, wireLength = 0;
  for (itNet = netToCell.begin(); itNet != netToCell.end(); ++itNet) {
    double minXW = mMaxX, minYW = mMaxY, maxXW = mMinX, maxYW = mMinY;
    for (itCellList = itNet -> second.begin(); itCellList != itNet -> second.end(); ++itCellList) {
      if(itCellList->name == "") {
        continue;
      }
      int orient = nodeId[itCellList->idx].orientation;
      xVal = nodeId[itCellList->idx].xBy2;
      yVal = nodeId[itCellList->idx].yBy2;

      if(orient == 0) { // 0
        xVal = xVal + itCellList->x_offset;
        yVal = yVal + itCellList->y_offset;
      } else if(orient == 2) { // 90
        xVal = xVal + itCellList->y_offset;
        yVal = yVal - itCellList->x_offset;
      } else if(orient == 4) { // 180
        xVal = xVal - itCellList->x_offset;
        yVal = yVal - itCellList->y_offset;
      } else if(orient == 6) { // 270
        xVal = xVal - itCellList->y_offset;
        yVal = yVal + itCellList->x_offset;
      } else {
        double rad = (orient*45.0*PI/180.0);
        xVal = itCellList->y_offset*sin(rad) + itCellList->x_offset*cos(rad) + xVal;
        yVal = itCellList->y_offset*cos(rad) - itCellList->x_offset*sin(rad) + yVal;
      }

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

/*
cell_overlap
Compute sum squared overlap for all components
*/
double GridBasedPlacer::cell_overlap() {
  double overlap = 0.0;

  for(size_t i = 0; i < nodeId.size(); i++) {
    //nodeId[i].printParameter();
    //if(nodeId[i].terminal) { continue; }
    for(size_t j = i; j < nodeId.size(); j++) {
      if (i == j) {continue;}
      //if(nodeId[j].terminal) { continue; }
      if(!intersects(nodeId[i].poly, nodeId[j].poly) || (nodeId[i].fixed && nodeId[j].fixed)) {
        continue;
      } else {
        double oa = 0.0;
        std::deque<polygon> intersect_poly;
        boost::geometry::intersection(nodeId[i].poly, nodeId[j].poly, intersect_poly);

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

/*
wireLength_partial
Compute HPWL for select nets
*/
double GridBasedPlacer::wirelength_partial(vector < Node *> &nodes, map<int, vector<pPin> > &netToCell) {
  vector < pPin > net;
  vector < int > ::iterator itNet;
  vector < pPin > ::iterator itCellList;
  vector < Node *>::iterator itNode;
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
        if(itCellList->name == "") {
          continue;
        }
        int orient = nodeId[itCellList->idx].orientation;
        xVal = nodeId[itCellList->idx].xBy2;
        yVal = nodeId[itCellList->idx].yBy2;

        if(orient == 0) { // 0
          xVal = xVal + itCellList->x_offset;
          yVal = yVal + itCellList->y_offset;
        } else if(orient == 2) { // 90
          xVal = xVal + itCellList->y_offset;
          yVal = yVal - itCellList->x_offset;
        } else if(orient == 4) { // 180
          xVal = xVal - itCellList->x_offset;
          yVal = yVal - itCellList->y_offset;
        } else if(orient == 6) { // 270
          xVal = xVal - itCellList->y_offset;
          yVal = yVal + itCellList->x_offset;
        } else {
          double rad = (orient*45.0*PI/180.0);
          xVal = itCellList->y_offset*sin(rad) + itCellList->x_offset*cos(rad) + xVal;
          yVal = itCellList->y_offset*cos(rad) - itCellList->x_offset*sin(rad) + yVal;
        }

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

/*
cell_overlap_partial
Compute sum squared overlap for select components
*/
double GridBasedPlacer::cell_overlap_partial(vector < Node *> &nodes) {
  double overlap = 0.0;
  unordered_set < int > cell_history;
  for(size_t i = 0; i < nodes.size(); i++) {
    cell_history.insert(nodes[i]->idx);
    if(rt) {
      //for(size_t j = 0; j < nodeId.size(); j++) {
      for ( rtree_t::const_query_iterator it = rtree.qbegin(index::intersects(nodeId[i].envelope)) ;
        it != rtree.qend() ; ++it ) {
        size_t j = it->second;
        if (cell_history.find(nodeId[j].idx) != cell_history.end()) {
          continue;
        }
        if(!intersects(nodes[i]->poly, nodeId[j].poly) || (nodes[i]->fixed && nodeId[j].fixed)) {
          continue;
        } else {
          double oa = 0.0;
          std::deque<polygon> intersect_poly;
          boost::geometry::intersection(nodes[i]->poly, nodeId[j].poly, intersect_poly);

          BOOST_FOREACH(polygon const& p, intersect_poly) {
              oa +=  bg::area(p);
          }
          if(oa < 1.0) {
              oa = 1.0;
          }
          overlap +=  pow(oa,2);
        }
      }
    } else {
      for(size_t j = 0; j < nodeId.size(); j++) {
        if (cell_history.find(nodeId[j].idx) != cell_history.end()) {
          continue;
        }
        if(!intersects(nodes[i]->poly, nodeId[j].poly) || (nodes[i]->fixed && nodeId[j].fixed)) {
          continue;
        } else {
          double oa = 0.0;
          std::deque<polygon> intersect_poly;
          boost::geometry::intersection(nodes[i]->poly, nodeId[j].poly, intersect_poly);

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
  }
  return overlap;
}

/*
rudy
Computes a routability score
*/
double GridBasedPlacer::rudy(map<int, vector<pPin> > &netToCell) {
  static bnu::matrix<double> D (static_cast<int>(abs(mMaxY)+abs(mMinY) + 1), static_cast<int>(abs(mMaxX)+abs(mMinX) + 1), 0.0);
  static bnu::matrix<double> D_route_sup (static_cast<int>(abs(mMaxY)+abs(mMinY) + 1), static_cast<int>(abs(mMaxX)+abs(mMinX) + 1), 1.0);

  D.clear();

  map<int, vector < pPin > > ::iterator itNet;
  vector < pPin > ::iterator itCellList;

  for (itNet = netToCell.begin(); itNet != netToCell.end(); ++itNet) {
    double xVal, yVal, hpwl = 0.0;
    double minXW = mMaxX, minYW = mMaxY, maxXW = mMinX, maxYW = mMinY;
    double rudy = 0.0;
    for (itCellList = itNet -> second.begin(); itCellList != itNet -> second.end(); ++itCellList) {
      int orient = nodeId[itCellList->idx].orientation;
 
      xVal = nodeId[itCellList->idx].xBy2;
      yVal = nodeId[itCellList->idx].yBy2;

      // compute pin position from orientation & offsets
      if(orient == 0) { // 0
        xVal = xVal + itCellList->x_offset;
        yVal = yVal + itCellList->y_offset;
      } else if(orient == 2) { // 90
        xVal = xVal + itCellList->y_offset;
        yVal = yVal - itCellList->x_offset;
      } else if(orient == 4) { // 180
        xVal = xVal - itCellList->x_offset;
        yVal = yVal - itCellList->y_offset;
      } else if(orient == 6) { // 270
        xVal = xVal - itCellList->y_offset;
        yVal = yVal + itCellList->x_offset;
      } else {
        double rad = (orient*45.0*PI/180.0);
        xVal = itCellList->y_offset*sin(rad) + itCellList->x_offset*cos(rad) + xVal;
        yVal = itCellList->y_offset*cos(rad) - itCellList->x_offset*sin(rad) + yVal;
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
      /*
      for (unsigned i = xVal-5; i < xVal+5; ++ i) {
          for (unsigned j = yVal-5; j < yVal+5; ++ j) {
            D (i,j) += 5;
          }
      }*/
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

  // write tab separated matrix to file
  /*
  if(iii % 10 == 0) {
      ofstream dat("cache/route/"+std::to_string(iii) + ".txt");
      for (unsigned i = 0; i < D.size1() ; i++) {
          for (unsigned j = 0; j < D.size2(); j++) {
              dat << D(i, j) << "\t";
          }
          dat << endl;
      }
  }*/

  return r;
}

double GridBasedPlacer::cost(map<int, vector<pPin> > &netToCell, int temp_debug) {
  double l2 = 1 - l1;
  double wirelength_cost = l1*0.9*(this->wirelength(netToCell) - wl_normalization.first)/(wl_normalization.second - wl_normalization.first);
  double overlap_cost = l2  * (this->cell_overlap() - area_normalization.first)/(area_normalization.second - area_normalization.first);
  double routability_cost = l1 * 0.1 * (this->rudy(netToCell) - routability_normalization.first)/(routability_normalization.second - routability_normalization.first);
  double total_cost = wirelength_cost + overlap_cost + routability_cost;
  if(debug > 1 || temp_debug == -1) {
    double wl = this->wirelength(netToCell);
    double cell_oa = this->cell_overlap();
    cout << "wirelength: " << this->wirelength(netToCell) << endl;
    cout << "overlap: " << this->cell_overlap() << endl;
    cout << "l1: " << l1 << " l2: " << l2 << endl;
    cout << "wirelength_cost: " << wirelength_cost  << endl;
    cout << "overlap_cost: " << overlap_cost  << endl;
    cout << "routability_cost: " << routability_cost  << endl;
    cout << "cost: " << total_cost << endl;

    if (cell_oa < 1 && wl < best_wl) {
        bestSol = nodeId;
        best_wl = wl;
    }
  }
  return total_cost;
}

double GridBasedPlacer::cost_partial(vector < Node *> &nodes, map<int, vector<pPin> > &netToCell) {
  double l2 = 1-l1;
  return l1 * 0.9 * (this->wirelength_partial(nodes, netToCell) - wl_normalization.first)/(wl_normalization.second - wl_normalization.first) +
         l2 * (this->cell_overlap_partial(nodes) - area_normalization.first)/(area_normalization.second - area_normalization.first) +
         l1 * 0.1 * (this->rudy(netToCell) - routability_normalization.first)/(routability_normalization.second - routability_normalization.first);

}

/*
random_node
Select random node from set of nodes
*/
vector < Node >::iterator GridBasedPlacer::random_node() {
  vector < Node > ::iterator itNode = nodeId.begin();
  int size = nodeId.size();
  boost::uniform_int<> uni_dist(0,size-1);
  boost::variate_generator<boost::mt19937&, boost::uniform_int<> > uni(rng, uni_dist);
  int randint = uni();
  std::advance(itNode, randint);
  return itNode;
}

/*
validate_move
validates a shift, project within board boundary
*/
void GridBasedPlacer::validate_move(Node &node, double rx, double ry) {
  rx = max(rx, mMinX);
  ry = max(ry, mMinY);

  rx = min(rx, mMaxX - node.width);
  ry = min(ry, mMaxY - node.height);

  node.setPos(rx,ry);
}

double c = 0.0;
double GridBasedPlacer::initiate_move(double current_cost, double Temperature, map<int, vector<pPin> > &netToCell) {
  // Initate a transition
  int state = -1;
  double prevCost = 0.0;
  if(debug > 1) {
    cout << "current_cost: " << current_cost << endl;
  }
  vector < Node > ::iterator rand_node1;
  vector < Node > ::iterator rand_node2;

  vector < Node* > perturbed_nodes;

  int r = 0;
  rand_node1 = this->random_node();
  while(rand_node1->terminal || rand_node1->name == "" || rand_node1->fixed) {
    rand_node1 = this->random_node();
  }
  perturbed_nodes.push_back(&(*rand_node1));

  if(rt) {
    rtree.remove(std::make_pair(rand_node1->envelope, rand_node1->idx));
  }
  double rand_node1_orig_x = rand_node1->xCoordinate;
  double rand_node1_orig_y = rand_node1->yCoordinate;
  double rand_node2_orig_x = 0.0;
  double rand_node2_orig_y = 0.0;

  if(debug > 1) {
    cout << "=======" << endl;
    cout <<  "name: " << rand_node1->name << endl;
  }

  boost::uniform_real<> uni_dist(0,1);
  boost::variate_generator<boost::mt19937&, boost::uniform_real<> > uni(rng, uni_dist);
  double i = uni();

  if (i< 0.2) { // swap
    state = 0;
    if(debug > 1) {
      cout << "swap" << endl;
    }
    rand_node2 = this->random_node();
    while(rand_node2->terminal || rand_node2->idx == rand_node1->idx || rand_node2->name == "" || rand_node2->fixed) {
      rand_node2 = this->random_node();
    }
    perturbed_nodes.push_back(&(*rand_node2));
    prevCost = this->cost_partial(perturbed_nodes,netToCell);

    if(rt) {
      rtree.remove(std::make_pair(rand_node2->envelope, rand_node2->idx));
    }
    rand_node2_orig_x = rand_node2->xCoordinate;
    rand_node2_orig_y = rand_node2->yCoordinate;

    validate_move(*rand_node1, rand_node2_orig_x, rand_node2_orig_y);
    validate_move(*rand_node2, rand_node1_orig_x, rand_node1_orig_y);
    if(rt) {
      rtree.insert(std::make_pair(rand_node1->envelope, rand_node1->idx));
      rtree.insert(std::make_pair(rand_node2->envelope, rand_node2->idx));
    }
  } else if (i < 0.85) { // shift
    state = 1;
    if(debug > 1) {
      cout << "shift" << endl;
    }
    prevCost = this->cost_partial(perturbed_nodes, netToCell);

    double sigma = rand_node1->sigma;

    boost::normal_distribution<> nd(0.0, sigma);
    boost::variate_generator<boost::mt19937&,
                             boost::normal_distribution<> > var_nor(rng, nd);

    double dx = var_nor();
    double dy = var_nor();

    double rx = rand_node1_orig_x + dx;
    double ry = rand_node1_orig_y + dy;

    this->validate_move(*rand_node1, rx, ry);
    if(rt) {
      rtree.insert(std::make_pair(rand_node1->envelope, rand_node1->idx));
    }
  } else { // rotate
    state = 2;
    if(debug > 1) {
      cout << "rotate" << endl;
    }
    prevCost = this->cost_partial(perturbed_nodes, netToCell);

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
    rand_node1->setRotation(r);
    this->validate_move(*rand_node1, rand_node1_orig_x, rand_node1_orig_y);

    if(rt) {
      rtree.insert(std::make_pair(rand_node1->envelope, rand_node1->idx));
    }
  }
  double transition_cost = this->cost_partial(perturbed_nodes,netToCell);
  double updated_cost = current_cost - prevCost + transition_cost;

  bool accept = this->check_move(current_cost, updated_cost, Temperature);

  if (!accept) {
    AcceptRate = 1/500 *(499*AcceptRate);
    if(debug > 1) {
      cout << "reject" << endl;
    }
    // revert state
    if(rt) {
      rtree.remove(std::make_pair(rand_node1->envelope, rand_node1->idx));
    }
    if (state == 0) {
      if(rt) {
        rtree.remove(std::make_pair(rand_node2->envelope, rand_node2->idx));
      }
      rand_node1->setPos(rand_node1_orig_x,rand_node1_orig_y);
      rand_node2->setPos(rand_node2_orig_x,rand_node2_orig_y);
      if(rt) {
        rtree.insert(std::make_pair(rand_node2->envelope, rand_node2->idx));
      }
    } else if (state == 1) {
      rand_node1->setPos(rand_node1_orig_x,rand_node1_orig_y);
    } else if (state == 2) {
      rand_node1->setRotation(8-r);
      rand_node1->setPos(rand_node1_orig_x,rand_node1_orig_y);
    }
    if(rt) {
      rtree.insert(std::make_pair(rand_node1->envelope, rand_node1->idx));
    }
    return current_cost;
  }
  else {
    if(debug > 1) {
      cout << "accept" << endl;
    }
    AcceptRate = 1/500 *(499*AcceptRate + 1);
    return updated_cost;
  }
}

/*
update_Temperature
Update the SA parameters according to annealing schedule
*/
void GridBasedPlacer::update_temperature(double Temperature) {
  vector < Node > ::iterator nodeit = nodeId.begin();
  l1 = 0.98*l1;
  if (Temperature > 50e-3) {
    Temperature = (0.985) * Temperature;
//    if (l1 > 92e-2) {
//      l1 -= 2e-4;
//    }
    for (nodeit = nodeId.begin(); nodeit != nodeId.end(); ++nodeit) {
      nodeit->sigma =  max(0.985*nodeit->sigma,3/4 * nodeit->sigma);
    }
  } else if (Temperature > 10e-3) {
    Temperature = (0.9992) * Temperature;
//    if (l1 > 88.5e-2) {
//      l1 -= 4e-4;
//    }
    for (nodeit = nodeId.begin(); nodeit != nodeId.end(); ++nodeit) {
      nodeit->sigma =  max(0.9992*nodeit->sigma,2/4 * nodeit->sigma);
    }
  } else if (Temperature > 50e-4) {
    Temperature = (0.9955) * Temperature;
//    if (l1 > 88e-2) {
//      l1 -= 8e-4;
//    }
    for (nodeit = nodeId.begin(); nodeit != nodeId.end(); ++nodeit) {
      nodeit->sigma =  max(0.9955*nodeit->sigma,1/4 * nodeit->sigma);
    }
  } else if (Temperature > 10e-4) {
    Temperature = (0.9965) * Temperature;
//    if (l1 > 85e-2) {
//      l1 -= 16e-4;
//    }
    for (nodeit = nodeId.begin(); nodeit != nodeId.end(); ++nodeit) {
      nodeit->sigma =  max(0.9965*nodeit->sigma,2.0);
    }
  } else {
    if (Temperature > 10e-8) {
      Temperature = (0.885) * Temperature;
//      if (l1 > 80e-2) {
//        l1 -= 10e-4;
//      }
    } else {
//      if (l1 > 10e-2) {
//        l1 -= 10e-4;
//      }
      Temperature = 0.0000000001;
    }
//    if (l1 > 10e-4) {
//      l1 -= 10e-4;
//    }
    for (nodeit = nodeId.begin(); nodeit != nodeId.end(); ++nodeit) {
      nodeit->sigma =  max(0.855*nodeit->sigma,1.0);
    }
  }
}

void GridBasedPlacer::update_accept_history(vector< double > &accept_ratio_history, float &accept_ratio) {
  if (accept_history.size() > 100) {
    accept_ratio = ((accept_ratio*100) - (accept_history[accept_history.size()-101]) + accept_history[accept_history.size()-1])/100;
  } else {
    accept_ratio = accumulate(accept_history.begin(), accept_history.end(),0.0) / accept_history.size();
  }
  accept_ratio_history.push_back(accept_ratio);
}

/*
check_move
either accept or reject the move based on current & previous temperature & cost
*/
bool GridBasedPlacer::check_move(double prevCost, double newCost, double Temperature) {
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

double GridBasedPlacer::initialize_temperature(double Temperature, map<int, vector<pPin> > &netToCell) {
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
      emax += exp(cost(netToCell)/t);
      this->initiate_move(0.0, Temperature, netToCell);
      emin += exp(cost(netToCell)/t);
    }
    xt = emax/emin;
    t = t * pow(log(xt),1/p)/log(x0);
  }
  return t;
}

/*
gen_report
generates a report and outputs files to ./reports/ directory
*/
void GridBasedPlacer::gen_report(map<string, vector<double> > &report,
                vector< double > &accept_ratio_history,
                map<int, vector<pPin> > &netToCell) {
    vector < double > cost_hist = report["cost_hist"];
    vector < double > wl_hist   = report["wl_hist"];
    vector < double > oa_hist   = report["oa_hist"];

    time_t t = time(0);
    struct tm * now = localtime( & t );

    char buf [80];
    strftime (buf,80,"%Y-%m-%d-%H-%M-%S",now);
    string buffer = string(buf);

    double wl = this->wirelength(netToCell);
    double oa = this->cell_overlap();
    double routability = this->rudy(netToCell);
    double normalized_wl = (wl - wl_normalization.first)/(wl_normalization.second - wl_normalization.first);
    double normalized_oa = (oa - area_normalization.first)/(area_normalization.second - area_normalization.first);
    double normalized_routability = routability;
    double l2 = 1 - l1;
    double cost = l1 * normalized_wl +
                  l2 * 0.9 * normalized_oa +
                  l2 * 0.1 * normalized_routability;

    std::ofstream f0("./reports/"+buffer+"_summary.txt");
    f0 << "wirelength: " << wl << '\n';
    f0 << "overlap: " << oa << '\n';
    f0 << "normalized wirelength: " << normalized_wl << '\n';
    f0 << "normalized overlap: " << normalized_oa << '\n';
    f0 << "routability: " << routability << '\n';
    f0 << "cost: " << cost << '\n';
    f0.close();

    std::ofstream f1("./reports/"+buffer+"_cost.txt");
    for(vector<double>::const_iterator i = cost_hist.begin(); i != cost_hist.end(); ++i) {
        f1 << *i << '\n';
    }
    f1.close();

    std::ofstream f2("./reports/"+buffer+"_wl.txt");
    for(vector<double>::const_iterator i = wl_hist.begin(); i != wl_hist.end(); ++i) {
        f2 << *i << '\n';
    }
    f2.close();

    std::ofstream f3("./reports/"+buffer+"_oa.txt");
    for(vector<double>::const_iterator i = oa_hist.begin(); i != oa_hist.end(); ++i) {
        f3 << *i << '\n';
    }
    f3.close();

    std::ofstream f4("./reports/"+buffer+"_accept_ratio.txt");
    for(vector<double>::const_iterator i = accept_ratio_history.begin(); i != accept_ratio_history.end(); ++i) {
        f4 << *i << '\n';
    }
    f4.close();
}

/*
annealer
main loop for sa algorithm
*/
float GridBasedPlacer::annealer(map<int, vector<pPin> > &netToCell, string initial_pl) {
  double Temperature = t_0;
  int num_components = 0;

  vector < Node > ::iterator itNode;
  map < string, vector < double > > report;
  report["cost_hist"] = cost_hist;
  report["wl_hist"] = wl_hist;
  report["oa_hist"] = oa_hist;

  vector< int > accept_history;
  float accept_ratio = 0.0;
  vector< double > accept_ratio_history;

  cout << "calculating initial params..." << endl;
  this->initialize_params(netToCell);
  if (initial_pl != "") {
    readPlFile(initial_pl);
  } else {
    this->random_initial_placement();
  }

  double cst = this->cost(netToCell,-1);
  if(var) {
    Temperature = this->initialize_temperature(Temperature, netToCell);
  }
  int idx = 0;

  for (itNode = nodeId.begin(); itNode != nodeId.end(); ++itNode) {
    rtree.insert(std::make_pair(itNode -> envelope, idx));
    idx+=1;
    if(!itNode -> terminal && !itNode -> fixed) {
      num_components += 1;
    }
  }

  high_resolution_clock::time_point t1 = high_resolution_clock::now();
  long long int ii = 0; // outer loop iterator
  int i = 0; // inner loop iterator
  while (ii < outer_loop_iter) {
    i = inner_loop_iter*num_components; 
    if(ii > 1) {
      high_resolution_clock::time_point t2 = high_resolution_clock::now();
      duration<double> time_span = duration_cast< duration<double> >(t2 - t1);

      cout << "******" << ii << "******" << endl;
      cout << "iteration: " << ii << endl;
      cout << "time: " <<  time_span.count() << " (s)" << endl;
      cout << "move/time: " <<  i*ii/time_span.count() << endl;
      cout << "time remaining: " <<  time_span.count()/ii * (outer_loop_iter-ii) << " (s)" << endl;
      cout << "temperature: " << Temperature << endl;
      cout << "acceptance ratio: " << accept_ratio << endl;
      if (ii % 25 == 0) {
          this->cost(netToCell,-1);
      }

      //this->gen_report(report,
      //           accept_ratio_history,
      //           netToCell);
    }
    while (i > 0) {

      cst = this->initiate_move(cst, Temperature, netToCell);
      this->update_accept_history(accept_ratio_history, accept_ratio);
      cost_hist.push_back(cst);
      wl_hist.push_back((wirelength(netToCell) - wl_normalization.first)/(wl_normalization.second - wl_normalization.first));
      oa_hist.push_back((cell_overlap() - area_normalization.first)/(area_normalization.second - area_normalization.first));

      i -= 1;
    }

    // convergence criterion
    if(eps > 0 && abs(cost_hist.end()[-1] - cost_hist.end()[-2]) < eps) {
      break;
    }
    this->update_temperature(Temperature);
    ii += 1;
  }
  return this->cost(netToCell);
}

int GridBasedPlacer::fact(int n) {
    if (n <= 1) return 1;
    else return n*fact(n-1);
}
