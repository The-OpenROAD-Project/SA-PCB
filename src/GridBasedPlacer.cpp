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
typedef bg::model::box< bg::model::d2::point_xy<double> > box;
typedef boost::geometry::model::polygon<boost::geometry::model::d2::point_xy<double> > polygon;

int debug = 1;
GridBasedPlacer::boundaries b;
boost::mt19937 rng;

vector < Node > nodeId;
map < string, int > name2id;
//bgi::rtree<std::pair<box, int>, bgi::quadratic<16>> rtree;

float l1 = 0.4;

//int main(int argc, char *argv[]) {
void GridBasedPlacer::test_placer() {
  srand(time(NULL));
  int opt;
  int outer_loop_iter = 105;
  int inner_loop_iter = 2;
  double eps = -1.0;
  double t_0 = 1.0;
  string out_file;
  int rotate_flag = 0;

  long long int idx = -1;

  string parg = "";
  string initial_pl = "";
  int bound_inc = 0;
  bool var = false;

  std::cout << "=================test_placer==================" << std::endl;
  /*
  while ((opt = getopt(argc,argv,"i:x:j:t:f:p:b:r::h")) != EOF) {
      switch(opt) {
          case 'x': idx = atoi(optarg); break;
          case 'p': parg=string(optarg); break;
          case 'l': initial_pl=string(optarg); break;
          case 'd': debug=atoi(optarg); break;
          case 'i': outer_loop_iter = atoi(optarg); break;
          case 'j': inner_loop_iter = atoi(optarg); break;
          case 't': t_0 = atof(optarg); break;
          case 'e': eps = atof(optarg); break;
          case 'v': var = true; break;
          case 'f': out_file=string(optarg); break;
          case 'r': rotate_flag=atoi(optarg); break;
          case 'b':
            bound_inc += 1;
            switch(bound_inc) { // TODO: have 4 different argument flags for each point
              case 1: b.minX = atoi(optarg) - 1;
              case 2: b.minY = atoi(optarg) - 1;
              case 3: b.maxX = atoi(optarg) + 1;
              case 4: b.maxY = atoi(optarg) + 1;
            }
          case 'h': fprintf(stderr, "./sa usuage is \n \
                                    -i <optional, value> : for denoting # outer iterations PER SA INSTANCE \n \
                                    -j <optional, value> : for denoting 'j'*#nodes inner iterations \n \
                                    -t <optional, value> : for denoting initial temperature \n \
                                    -f <optional, str>   : for output filename \n \
                                    -e <optional, float> : convergence epsilon \n \
                                    -v <optional>        : ben-amur flag \n \
                                    -x <optional, int>   : simulated annealing instance index \n \
                                    -p <required, string>: prefix of board to place \n \
                                    -l <optional, string>: prefix of initial placement \n \
                                    -d <optional, {0-3}> : debug verbosity \n \
                                    -r <optional, {0-3}> : rotation \n \
                                    EXAMPLE: ./sa -i 20000 -j 20 -t 1 -p input -f output");
          default: cout<<endl; abort();
      }
  }
  if(parg == "")   {
    cout << "./main: option requires an argument -- p\n";
    exit(1);
  }
  if(parg == "")   {
    out_file = std::to_string( idx );
  } 
  */

  
  parg = "../designs/bm1";
  cout << "circuit: " << parg << endl;

  string nodesfname  = parg + ".nodes";
  string shapesfname = parg + ".shapes";
  string netsfname   = parg + ".nets";
  string plfname     = "";

  if (initial_pl != "") {
    plfname          = initial_pl + ".pl";
  } else {
    plfname          = parg + ".pl";
  }

  string wtsfname    = parg + ".wts";

  cout << nodesfname << endl;
  cout << shapesfname << endl;
  cout << netsfname << endl;
  cout << plfname << endl;
  
  map<int, vector<Pin> > netToCell;

  if(debug) { cout << "reading nodes..." << endl; }
  //readNodesFile(nodesfname);
  //readShapesFile(shapesfname);
  //readWtsFile(wtsfname);
  //if(debug) { cout << "reading pl..." << endl; }
  //readPlFile(plfname);

  //mDb.printInst();
  std::vector<instance> &instances =  mDb.getInstances();
  std::vector<net> &nets = mDb.getNets();
  
  for (int i = 0; i < instances.size(); ++i) {
    Node n;
    nodeId.push_back(n);
  }
  for (auto &inst : instances) {
      Node n;
      n.setParameterNodes(inst.getName(), inst.getWidth(), inst.getHeight(), 1, inst.getId());
      //nodeId.push_back(n);
      nodeId[inst.getId()] = n;
      name2id.insert(pair < string, int > (inst.getName(), inst.getId()));
      nodeId[name2id[inst.getName()]].setParameterPl(inst.getX() - inst.getWidth()/2, inst.getY() - inst.getHeight()/2, "N", 0);
  }
  if(debug) { cout << "reading nets..." << endl; }
  netToCell = readNetsFile(netsfname);
  if(debug) { cout << "calculating boundaries..." << endl; }
  set_boundaries();
  if(debug) { cout << "annealing" << endl; }

    for (auto &net : nets){

      vector < Pin > pinTemp;
      for (auto &pin : net.getPins())
      {
          Pin p;
          auto &inst = mDb.getInstance(pin.m_inst_id);
          nodeId[name2id[inst.getName()]].setNetList(net.getId());

          auto &comp = mDb.getComponent(pin.m_comp_id);
          point_2d pos;
          mDb.getPinPosition(pin, &pos);

          auto &pad = comp.getPadstack(pin.m_padstack_id);
          //cout << pos.m_x << " " << inst.getX() <<  " " << nodeId[name2id[inst.getName()]].xCoordinate << " " <<  pos.m_x - nodeId[name2id[inst.getName()]].xCoordinate << endl;
          p.set_params(inst.getName(), pos.m_x - nodeId[name2id[inst.getName()]].xBy2, pos.m_y - nodeId[name2id[inst.getName()]].yBy2, nodeId[name2id[inst.getName()]].idx);
          pinTemp.push_back(p);

      }
      netToCell.insert(pair < int, vector< Pin > > (net.getId(), pinTemp));
    }

  cout << "annealing" << endl;
  //float cost = this->annealer(outer_loop_iter, inner_loop_iter, eps, t_0, var, netToCell, initial_pl, rotate_flag);
  writePlFile("./final_placement.pl");
  //cout << " " << idx << " " << cost << endl;
}

/*
initialize_params
Empirically finds normalization parameters for scaling cost terms.
Currently, we scale by 1/(f) where f is the cost of an expected placement.
*/
void GridBasedPlacer::initialize_params(std::pair <double,double> &wl_normalization,
                       std::pair <double,double> &area_normalization,
                       std::pair <double,double> &routability_normalization,
                       map<int, vector<Pin> > &netToCell,
                       int rotate_flag) {

  vector < std::pair <double,double> > normalization_terms;

  int num_components = 0;
  vector < Node > ::iterator itNode;
  for (itNode = nodeId.begin(); itNode != nodeId.end(); ++itNode) {
    if(!itNode -> terminal) {
      num_components += 1;
    }
  }
  // min/max wl
  std::pair < double, double > wl;
//  double min_wl = std::numeric_limits<double>::max();
//  double max_wl = 0.0;

  // min/max overlap area
  std::pair < double, double > area;
//  double min_area = 0.0;
//  double max_area = 0.0;

  std::pair < double, double > rn;
//  double min_rudy = 0.0;
//  double max_rudy = 0.0;

/*
  vector < Node > ::iterator nodeit1;
  vector < Node > ::iterator nodeit2;

  for (nodeit1 = nodeId.begin(); nodeit1 != nodeId.end(); ++nodeit1) {
    if(nodeit1->terminal) {
      continue;
    } else {
        double cell_wl = nodeit1->width + nodeit1->height;
        if (cell_wl < min_wl) {
          min_wl = cell_wl;
        }
        for (nodeit2 = nodeit1++; nodeit2 != nodeId.end(); ++nodeit2) {
          if(nodeit2->idx == nodeit1->idx){continue;}
          if(nodeit1->fixed && nodeit2->fixed) {
              if(intersects(nodeit1->poly, nodeit2->poly)) {
                double oa = 0.0;
                std::deque<polygon> intersect_poly;
                boost::geometry::intersection(nodeit1->poly, nodeit2->poly, intersect_poly);

                BOOST_FOREACH(polygon const& p, intersect_poly) {
                    oa += bg::area(p);
                }
                min_area += pow(oa,2);
                max_area += pow(oa,2);
              }
            } else {
              max_area += pow(min(nodeit1->width, nodeit2->width) * min(nodeit1->height, (nodeit2->width)), 2);
            }
        }
      }
  }

  min_wl = min_wl*num_components;
  max_wl = ((b.maxX - b.minX) + (b.maxY - b.minY))*netToCell.size();

  wl.first = 0.0; // min_area
  wl.second = max(wirelength(netToCell), 1.0); // max_area

  area.first = 0.0; // min_area;
  area.second = max(cell_overlap(), 1.0); // max_area;

  routability_normalization.first = 0.0;
  routability_normalization.second = max(rudy(netToCell),1.0);
*/


  double sum_wl = 0.0;
  double sum_oa = 0.0;
  double sum_rn = 0.0;
  for (int i = 0; i<1000; i++) {
    this->random_initial_placement(rotate_flag);
    sum_wl += this->wirelength(netToCell);
    sum_oa += this->cell_overlap();
    sum_rn += this->rudy(netToCell);
  }

  wl.first = 0.0;
  wl.second = sum_wl / 1000.0;

  area.first = 0.0;
  area.second = sum_oa / 1000.0;

  rn.first = 0.0;
  rn.second = sum_rn / 1000.0;

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
      if (xval < b.minX) {
        b.minX = xval - 2;
      }
      if (xval > b.maxX) {
        b.maxX = xval + width + 2;
      }
      if (yval < b.minY) {
        b.minY = yval - 2;
      }
      if (yval > b.maxY) {
        b.maxY = yval + height + 2;
      }
    }
  }

  if (b.minX == 0 && b.maxX == 0 && b.minY == 0 && b.maxY == 0) {
    cout << "no boundary found" << endl;
    b.maxX = sqrt(1.8*total_area);
    b.maxY = sqrt(1.8*total_area);
  }

  if(debug) {
    cout << "min: " << b.minX << " " << b.minY << endl;
    cout << "max: " << b.maxX << " " << b.maxY << endl;
  }
}

/*
random_placement
Randomly place and orient a single component within a bounded region
*/
void GridBasedPlacer::random_placement(int xmin, int xmax, int ymin, int ymax, Node &n, int rotate_flag) {
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
void GridBasedPlacer::random_initial_placement(int rotate_flag) {
  vector < Node > ::iterator itNode;
  for (itNode = nodeId.begin(); itNode != nodeId.end(); ++itNode) {
    if (!itNode -> fixed) {
      this->random_placement(b.minX, b.maxX - max(itNode -> width, itNode -> height), b.minY, b.maxY - max(itNode -> width, itNode -> height), *itNode, rotate_flag);
    }
  }
}

/*
project_soln
projects the final solution to enforce 0 overlap by iteratively
shifting overlapped components until cumulative overlap reduces below an eps.
*/
void GridBasedPlacer::project_soln() {
  // ensure zero overlap
  double oa = cell_overlap();
  double eps = 10.0e-5;

  vector < Node > ::iterator nodeit1;
  vector < Node > ::iterator nodeit2;

  box env1;
  box env2;

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
double GridBasedPlacer::wirelength(map<int, vector<Pin> > &netToCell) {
  // compute HPWL
  map<int, vector < Pin > > ::iterator itNet;
  vector < Pin > ::iterator itCellList;
  double xVal, yVal, wireLength = 0;
  for (itNet = netToCell.begin(); itNet != netToCell.end(); ++itNet) {
    double minXW = b.maxX, minYW = b.maxY, maxXW = b.minX, maxYW = b.minY;
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
double GridBasedPlacer::wirelength_partial(vector < Node *> &nodes, map<int, vector<Pin> > &netToCell) {
  //map < int, vector < Pin > > ::iterator itNet;
  vector < Pin > net;
  vector < int > ::iterator itNet;
  vector < Pin > ::iterator itCellList;
  vector < Node *>::iterator itNode;
  unordered_set < int > net_history;

  double xVal, yVal, wireLength = 0;
  //for (itNet = netToCell.begin(); itNet != netToCell.end(); ++itNet) {
  for (itNode = nodes.begin(); itNode != nodes.end(); ++itNode) {
    for (itNet = (*itNode)->Netlist.begin(); itNet != (*itNode)->Netlist.end(); ++itNet) {
      if (net_history.find(*itNet) == net_history.end()) {
        net_history.insert(*itNet);
      } else {
        continue;
      }
      net = netToCell.find(*itNet)->second;
      double minXW = b.maxX, minYW = b.maxY, maxXW = b.minX, maxYW = b.minY;
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
    //if(nodes[i].terminal) { continue; }
    cell_history.insert(nodes[i]->idx);
    for(size_t j = 0; j < nodeId.size(); j++) {
      //if(nodeId[j].terminal) { continue; }
      /*
      for ( Rtree::const_query_iterator it = rtree.qbegin(index::intersects(nodeId[i].envelope)) ;
        it != rtree.qend() ; ++it ) {
          // yeet
      }*/
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
  return overlap;
}

/*
rudy
Computes a routability score
*/
double GridBasedPlacer::rudy(map<int, vector<Pin> > &netToCell) {
  static bnu::matrix<double> D (static_cast<int>(abs(b.maxY)+abs(b.minY) + 1), static_cast<int>(abs(b.maxX)+abs(b.minX) + 1), 0.0);
  static bnu::matrix<double> D_route_sup (static_cast<int>(abs(b.maxY)+abs(b.minY) + 1), static_cast<int>(abs(b.maxX)+abs(b.minX) + 1), 1.0);

  D.clear();

  map<int, vector < Pin > > ::iterator itNet;
  vector < Pin > ::iterator itCellList;

  for (itNet = netToCell.begin(); itNet != netToCell.end(); ++itNet) {
    double xVal, yVal, hpwl = 0.0;
    double minXW = b.maxX, minYW = b.maxY, maxXW = b.minX, maxYW = b.minY;
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

double GridBasedPlacer::cost(
            std::pair <double,double> &wl_normalization,
            std::pair <double,double> &area_normalization,
            std::pair <double,double> &routability_normalization,
            map<int, vector<Pin> > &netToCell,
            int temp_debug) {
  double l2 = 1 - l1;
  double wirelength_cost = l1*(this->wirelength(netToCell) - wl_normalization.first)/(wl_normalization.second - wl_normalization.first);
  double overlap_cost = l2  * 0.9 * (this->cell_overlap() - area_normalization.first)/(area_normalization.second - area_normalization.first);
  double routability_cost = l1 * 0.1 * (this->rudy(netToCell) - routability_normalization.first)/(routability_normalization.second - routability_normalization.first);
  double total_cost = wirelength_cost + overlap_cost + routability_cost;
  if(debug > 1 || temp_debug == -1) {
    cout << "wirelength: " << this->wirelength(netToCell) << endl;
    cout << "overlap: " << this->cell_overlap() << endl;
    cout << "l1: " << l1 << " l2: " << l2 << endl;
    cout << "wirelength_cost: " << wirelength_cost  << endl;
    cout << "overlap_cost: " << overlap_cost  << endl;
    cout << "routability_cost: " << routability_cost  << endl;
    cout << "cost: " << total_cost << endl;
  }
  return total_cost;
}

double GridBasedPlacer::cost_partial(vector < Node *> &nodes,
                    std::pair <double,double> &wl_normalization,
                    std::pair <double,double> &area_normalization,
                    std::pair <double,double> &routability_normalization,
                    map<int, vector<Pin> > &netToCell) {
  double l2 = 1-l1;
  return l1 * 0.9 * (this->wirelength_partial(nodes, netToCell) - wl_normalization.first)/(wl_normalization.second - wl_normalization.first) +
         l2 * this->cell_overlap() + //(cell_overlap_partial(nodes) - area_normalization.first)/(area_normalization.second - area_normalization.first) +
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
  rx = max(rx, b.minX);
  ry = max(ry, b.minY);

  rx = min(rx, b.maxX - node.width);
  ry = min(ry, b.maxY - node.height);

  node.setPos(rx,ry);
}

double c = 0.0;
double GridBasedPlacer::initiate_move(double current_cost,
                     vector< int > &accept_history,
                     double & Temperature,
                     std::pair <double,double> &wl_normalization,
                     std::pair <double,double> &area_normalization,
                     std::pair <double,double> &routability_normalization,
                     map<int, vector<Pin> > &netToCell,
                     int rotate_flag) {
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

  //rtree.remove(std::make_pair(rand_node1->envelope, rand_node1->idx));
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
    prevCost = this->cost_partial(perturbed_nodes, wl_normalization, area_normalization, routability_normalization, netToCell);

    //rtree.remove(std::make_pair(rand_node2->envelope, rand_node2->idx));
    rand_node2_orig_x = rand_node2->xCoordinate;
    rand_node2_orig_y = rand_node2->yCoordinate;

    validate_move(*rand_node1, rand_node2_orig_x, rand_node2_orig_y);
    validate_move(*rand_node2, rand_node1_orig_x, rand_node1_orig_y);
    //rtree.insert(std::make_pair(rand_node1->envelope, rand_node1->idx));
    //rtree.insert(std::make_pair(rand_node2->envelope, rand_node2->idx));
  } else if (i < 0.85) { // shift
    state = 1;
    if(debug > 1) {
      cout << "shift" << endl;
    }
    prevCost = this->cost_partial(perturbed_nodes, wl_normalization, area_normalization, routability_normalization, netToCell);

    double sigma = rand_node1->sigma;

    boost::normal_distribution<> nd(0.0, sigma);
    boost::variate_generator<boost::mt19937&,
                             boost::normal_distribution<> > var_nor(rng, nd);

    double dx = var_nor();
    double dy = var_nor();

    double rx = rand_node1_orig_x + dx;
    double ry = rand_node1_orig_y + dy;

    this->validate_move(*rand_node1, rx, ry);
    //rtree.insert(std::make_pair(rand_node1->envelope, rand_node1->idx));
  } else { // rotate
    state = 2;
    if(debug > 1) {
      cout << "rotate" << endl;
    }
    prevCost = this->cost_partial(perturbed_nodes, wl_normalization, area_normalization, routability_normalization, netToCell);

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

    //rtree.insert(std::make_pair(rand_node1->envelope, rand_node1->idx));
  }
  double transition_cost = this->cost_partial(perturbed_nodes, wl_normalization, area_normalization, routability_normalization, netToCell);
  double updated_cost = current_cost - prevCost + transition_cost;

  bool accept = this->check_move(current_cost,
                           updated_cost,
                           accept_history,
                           Temperature);
  if (!accept) {
    if(debug > 1) {
      cout << "reject" << endl;
    }
    // revert state
    //rtree.remove(std::make_pair(rand_node1->envelope, rand_node1->idx));
    if (state == 0) {
      //rtree.remove(std::make_pair(rand_node2->envelope, rand_node2->idx));
      rand_node1->setPos(rand_node1_orig_x,rand_node1_orig_y);
      rand_node2->setPos(rand_node2_orig_x,rand_node2_orig_y);
      //rtree.insert(std::make_pair(rand_node2->envelope, rand_node2->idx));
    } else if (state == 1) {
      rand_node1->setPos(rand_node1_orig_x,rand_node1_orig_y);
    } else if (state == 2) {
      rand_node1->setRotation(8-r);
      rand_node1->setPos(rand_node1_orig_x,rand_node1_orig_y);
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

/*
update_Temperature
Update the SA parameters according to annealing schedule
*/
void GridBasedPlacer::update_temperature(double& Temperature) {
  vector < Node > ::iterator nodeit = nodeId.begin();
  if (Temperature > 50e-3) {
    Temperature = (0.985) * Temperature;
    if (l1 > 92e-2) {
      l1 -= 2e-4;
    }
    for (nodeit = nodeId.begin(); nodeit != nodeId.end(); ++nodeit) {
      nodeit->sigma =  max(0.985*nodeit->sigma,3/4 * nodeit->sigma);
    }
  } else if (Temperature > 10e-3) {
    Temperature = (0.9992) * Temperature;
    if (l1 > 88.5e-2) {
      l1 -= 4e-4;
    }
    for (nodeit = nodeId.begin(); nodeit != nodeId.end(); ++nodeit) {
      nodeit->sigma =  max(0.9992*nodeit->sigma,2/4 * nodeit->sigma);
    }
  } else if (Temperature > 50e-4) {
    Temperature = (0.9955) * Temperature;
    if (l1 > 88e-2) {
      l1 -= 8e-4;
    }
    for (nodeit = nodeId.begin(); nodeit != nodeId.end(); ++nodeit) {
      nodeit->sigma =  max(0.9955*nodeit->sigma,1/4 * nodeit->sigma);
    }
  } else if (Temperature > 10e-4) {
    Temperature = (0.9965) * Temperature;
    if (l1 > 85e-2) {
      l1 -= 16e-4;
    }
    for (nodeit = nodeId.begin(); nodeit != nodeId.end(); ++nodeit) {
      nodeit->sigma =  max(0.9965*nodeit->sigma,2.0);
    }
  } else {
    if (Temperature > 10e-8) {
      Temperature = (0.885) * Temperature;
      if (l1 > 80e-2) {
        l1 -= 10e-4;
      }
    } else {
      if (l1 > 10e-2) {
        l1 -= 10e-4;
      }
      Temperature = 0.0000000001;
    }
    if (l1 > 10e-4) {
      l1 -= 10e-4;
    }
    for (nodeit = nodeId.begin(); nodeit != nodeId.end(); ++nodeit) {
      nodeit->sigma =  max(0.855*nodeit->sigma,1.0);
    }
  }
}

void GridBasedPlacer::update_accept_history(vector< int > &accept_history, vector< double > &accept_ratio_history, float &accept_ratio) {
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
bool GridBasedPlacer::check_move(double prevCost,
                double newCost,
                vector< int > &accept_history,
                double & Temperature) {
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

double GridBasedPlacer::initialize_temperature(vector< int > &accept_history,
                              double &Temperature,
                              std::pair <double,double> &wl_normalization,
                              std::pair <double,double> &area_normalization,
                              std::pair <double,double> &routability_normalization,
                              map<int, vector<Pin> > &netToCell,
                              int rotate_flag) {
  double t = 0.0;
  double emax = 0.0;
  double emin = 0.0;
  double xt = 1.0;
  double x0 = 0.84;
  double p = 2.0;
  this->random_initial_placement(rotate_flag);
  for(int i=1; i<=10; i++){
    for(int j=1; j<=10; j++){
      this->random_initial_placement(rotate_flag);
      emax += exp(cost(wl_normalization, area_normalization, routability_normalization, netToCell)/t);
      this->initiate_move(0.0, accept_history, Temperature, wl_normalization, area_normalization, routability_normalization, netToCell, rotate_flag);
      emin += exp(cost(wl_normalization, area_normalization, routability_normalization, netToCell)/t);
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
                std::pair <double,double> &wl_normalization,
                std::pair <double,double> &area_normalization,
                std::pair <double,double> &routability_normalization,
                map<int, vector<Pin> > &netToCell) {
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

    //writePlFile("./reports/"+buffer+"_pl.pl");
}

/*
annealer
main loop for sa algorithm
*/
float GridBasedPlacer::annealer(int outer_loop_iter,
               int inner_loop_iter,
               double eps,
               double t_0,
               bool var,
               map<int, vector<Pin> > &netToCell,
               string initial_pl,
               int rotate_flag) {
  double Temperature = t_0;
  int num_components = 0;

  vector < Node > ::iterator itNode;
  map < string, vector < double > > report;
  vector < double > cost_hist;
  vector < double > wl_hist;
  vector < double > oa_hist;
  report["cost_hist"] = cost_hist;
  report["wl_hist"] = wl_hist;
  report["oa_hist"] = oa_hist;

  std::pair <double,double> wl_normalization;
  std::pair <double,double> area_normalization;
  std::pair <double,double> routability_normalization;

  vector< int > accept_history;
  float accept_ratio = 0.0;
  vector< double > accept_ratio_history;

  if(debug) { cout << "calculating initial params..." << endl; }
  cout << "calculating initial params..." << endl;
  this->initialize_params(wl_normalization, area_normalization, routability_normalization, netToCell, rotate_flag);
  cout << "initial palcement..." << endl;
  if (initial_pl != "") {
    readPlFile(initial_pl);
  } else {
    this->random_initial_placement(rotate_flag);
  }

  double cst = this->cost(wl_normalization, area_normalization, routability_normalization, netToCell);
  if(var) {
    Temperature = this->initialize_temperature(accept_history,
                                         Temperature,
                                         wl_normalization,
                                         area_normalization,
                                         routability_normalization,
                                         netToCell,
                                         rotate_flag);
  }
  int idx = 0;

  cout << "counting components..." << endl;
  for (itNode = nodeId.begin(); itNode != nodeId.end(); ++itNode) {
    //rtree.insert(std::make_pair(itNode -> envelope, idx));
    idx+=1;
    if(!itNode -> terminal && !itNode -> fixed) {
      num_components += 1;
    }
  }

  high_resolution_clock::time_point t1 = high_resolution_clock::now();
  long long int ii = 0; // outer loop iterator
  int i = 0; // inner loop iterator
  if(debug) { cout << "annealing..." << endl; }
  while (ii < outer_loop_iter) {
    i = inner_loop_iter*num_components;
    if(ii > 1 && ii % 10 == 0 && debug) {
      high_resolution_clock::time_point t2 = high_resolution_clock::now();
      duration<double> time_span = duration_cast< duration<double> >(t2 - t1);

      cout << "******" << ii << "******" << endl;
      cout << "iteration: " << ii << endl;
      cout << "time: " <<  time_span.count() << " (s)" << endl;
      cout << "move/time: " <<  i*ii/time_span.count() << endl;
      cout << "time remaining: " <<  time_span.count()/ii * (outer_loop_iter-ii) << " (s)" << endl;
      cout << "temperature: " << Temperature << endl;
      cout << "acceptance ratio: " << accept_ratio << endl;
      this->cost(wl_normalization, area_normalization, routability_normalization,netToCell,-1);

      this->gen_report(report,
                 accept_ratio_history,
                 wl_normalization,
                 area_normalization,
                 routability_normalization,
                 netToCell);
    }
    while (i > 0) {

      cst = this->initiate_move(cst, accept_history, Temperature, wl_normalization, area_normalization, routability_normalization, netToCell, rotate_flag);
      this->update_accept_history(accept_history, accept_ratio_history, accept_ratio);
      report["cost_hist"].push_back(cst);

      //wl_hist.push_back((wirelength(netToCell) - wl_normalization.first)/(wl_normalization.second - wl_normalization.first));
      //oa_hist.push_back((cell_overlap() - area_normalization.first)/(area_normalization.second - area_normalization.first));

      i -= 1;
    }

    // convergence criterion
    if(eps > 0 && abs(cost_hist.end()[-1] - cost_hist.end()[-2]) < eps) {
      break;
    }
    this->update_temperature(Temperature);
    ii += 1;
    if (ii % 10 == 0) {
      writePlFile("./cache/"+std::to_string( ii )+".pl");
    }
  }
  return this->cost(wl_normalization, area_normalization, routability_normalization, netToCell);
}

int GridBasedPlacer::fact(int n) {
    if (n <= 1) return 1;
    else return n*fact(n-1);
}
