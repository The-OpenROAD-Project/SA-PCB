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
vector < Module * > moduleId;
Hierarchy H;
map < string, int > name2id;

void GridBasedPlacer::test_hplacer_flow() {
  srand(time(NULL));
  int opt;
  string out_file;
  do_hplace = true;

  long long int idx = -1;

  string parg = "";
  string initial_pl = "../designs/ispd2005/adaptec2/adaptec2_sol";

  std::cout << "=================test_hplacer==================" << std::endl;
  
  parg = "../designs/ispd2005/adaptec2/adaptec2";
  cout << "circuit: " << parg << endl;

  string nodesfname  = parg + ".nodes";
  string netsfname   = parg + ".nets";
  string clstname = parg + ".clst";
  string plfname     = "";

  if (initial_pl != "") {
    plfname          = initial_pl + ".pl";
  } else {
    plfname          = parg + ".pl";
  }

  string wtsfname    = parg + ".wts";

  cout << nodesfname << endl;
  cout << netsfname << endl;
  cout << plfname << endl;
  cout << clstname << endl;
 
  map<int, vector<pPin> > netToCell;
  cout << "reading nodes..." << endl;
  readNodesFile(nodesfname);
  cout << "# Nodes: " << nodeId.size() << endl;
  cout << "reading pl..." << endl;
  readPlFile(plfname);
  cout << "reading nets..." << endl;
  netToCell = readNetsFile(netsfname); // mapping from net ids to vectors of pins - need to add weights
  cout << "# Nets: " << netToCell.size() << endl;
  cout << "reading cluster heirarchy..." << endl;
  readClstFile(clstname);
  cout << "computing netlist heirarchy..." << endl;
  H.set_netlist_hierarchy(netToCell);
  cout << "composing module geometries..." << endl;
  H.set_module_geometries(nodeId);

  cout << "calculating boundaries..." << endl;
  set_boundaries();

  cout << "annealing" << endl;
  hplace(netToCell, initial_pl);
}

void GridBasedPlacer::hplace(map<int, vector<pPin> > &netToCell, string initial_pl) {
  // top down placement
  map<int, vector<pPin> > tmp_netToCell = netToCell;
  //H.print_param(1); 
  int l = 0;
  for (auto &lvl : H.levels) {
    cout << "===== level: " << l << " =====" << endl;
    if (l == 0){ 
      l ++;
      continue;
    }
    moduleId = lvl.modules;
    map<int, vector<Module *> > netToCell = lvl.netToModule;

    float cost = h_annealer(netToCell, initial_pl, l);

    nodeId = H.update_cell_positions_at_level(nodeId, l);

    // update coords

    cout << "writting result" << endl;
    writePlFile("./final_placement_"+std::to_string(l)+".pl");

    std::ofstream f1("./oa_"+std::to_string(l)+".txt");
    for(vector<double>::const_iterator i = oa_hist.begin(); i != oa_hist.end(); ++i) {
        f1 << *i << '\n';
    }
    f1.close();

    std::ofstream f2("./wl_"+std::to_string(l)+".txt");
    for(vector<double>::const_iterator i = wl_hist.begin(); i != wl_hist.end(); ++i) {
        f2 << *i << '\n';
    }
    f2.close();

    std::ofstream f3("./l_"+std::to_string(l)+".txt");
    for(vector<double>::const_iterator i = l_hist.begin(); i != l_hist.end(); ++i) {
        f3 << *i << '\n';
    }
    f3.close();

    std::ofstream f4("./var_"+std::to_string(l)+".txt");
    for(vector<double>::const_iterator i = var_hist.begin(); i != var_hist.end(); ++i) {
        f4 << *i << '\n';
    }
    f4.close();

    std::ofstream f5("./temp_"+std::to_string(l)+".txt");
    for(vector<double>::const_iterator i = temp_hist.begin(); i != temp_hist.end(); ++i) {
        f5 << *i << '\n';
    }
    f5.close();

    wl_hist.clear();
    oa_hist.clear();
    l_hist.clear();
    density_hist.clear();
    var_hist.clear();
    temp_hist.clear();

    l++;
    initial_loop_iter = ceil(initial_loop_iter / 2);
    t_0 = 0.01*t_0;

    cout << endl;
  } 
/*
  cout << "===== greedy placement =====" << endl;
  t_0 = 0.0;
  initial_loop_iter = 1;
  nodeId = H.update_cell_positions(nodeId);
  vector < Module * > moduleId_tmp;
  map<int, vector<Module *> > netToCell_tmp;
  for (auto &node : nodeId) {
    Module *tmp = new Module;
    tmp->init_module(node.idx, -1, false);
    tmp->setParameterNodes(node.width, node.height);
    tmp->setParameterPl(node.xCoordinate, node.yCoordinate);
    tmp->terminal = node.terminal;
    tmp->fixed = node.fixed;
    moduleId_tmp.push_back(tmp);
  }
  moduleId = moduleId_tmp;
  int netidx = 1;
  for (auto &net : netToCell) { 
    int netid = net.first;
    vector<Pin> pvec = net.second; //pvec_i.idx -> cell id
    vector<Module *> ms; 
    for (auto &pin : pvec) {
      int cell_id = pin.idx;
      Module *m = moduleId[cell_id-1];
      m->setNetList(netidx);
      ms.push_back(m);
    }
    netToCell_tmp.insert(pair < int, vector < Module * > > (netidx, ms));
    netidx ++;
  }
  float cost = this->annealer(netToCell_tmp, initial_pl);*/
}

kicadPcbDataBase &GridBasedPlacer::test_placer_flow() {
    srand(time(NULL));
    int opt;
    string out_file;

    long long int idx = -1;

    string parg = "";
    string initial_pl = "";

    std::cout << "=================test_placer==================" << std::endl;

    map<int, vector<pPin> > netToCell;

    std::vector<instance> &instances =  mDb.getInstances();
    std::vector<net> &nets = mDb.getNets();
    for (int i = 0; i < instances.size(); ++i) {
      Node n;
      nodeId.push_back(n);
    }

    bool fixed = false;
    int mirror = 1;
    for (auto &inst : instances) {
        point_2d bbox;
        mDb.getCompBBox(inst.getComponentId(), &bbox); 

        Node n;
        int layer = inst.getLayer();
        if(layer == 0) {
            mirror = 1;
        } else {
            mirror = -1;
        }
       
        n.setParameterNodes(inst.getName(), bbox.m_x + 0.2, bbox.m_y + 0.2, false, inst.getId(), mirror);
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
            fixed = true;
        } else {
            fixed = false;
        }

        nodeId[name2id[inst.getName()]].setParameterPl(inst.getX() - bbox.m_x/2, inst.getY() - bbox.m_y/2, ang, fixed);
        nodeId[name2id[inst.getName()]].printParameter();
    }

    cout << "calculating boundaries..." << endl;
    double pminx, pmaxx, pminy, pmaxy;
    mDb.getBoardBoundaryByPinLocation(pminx, pmaxx, pminy ,pmaxy);
    mMinX = pminx;
    mMaxX = pmaxx;
    mMinY = pminy;
    mMaxY = pmaxy;
    cout << mMinX << "," << mMinY << " " << mMaxX << "," << mMaxY << endl;

    for (auto &net : nets) {
      vector < pPin > pinTemp;
      for (auto &pin : net.getPins()) {
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
    float cost = annealer(netToCell, initial_pl);
    if (bestSol.size() > 0) {
        nodeId = bestSol; 
    }
    cout << "Solution: " << "overlap: " << cell_overlap() << " wirelength: " << wirelength(netToCell) << endl;
    // write back to db
    cout << "writing back to db..." << endl;
    int top = 0;
    int bot = 31;
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
        double cx = nodeId[inst.getId()].xBy2 + 0.1;
        double cy = nodeId[inst.getId()].yBy2 + 0.1;
        inst.setX(cx);
        inst.setY(cy);

        if (nodeId[inst.getId()].layer == 1) {
          inst.setLayer(top);
        } else {
          inst.setLayer(bot);
        }
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

  for (int i = 0; i<initial_loop_iter; i++) {
    random_initial_placement();
    sum_wl += wirelength(netToCell);
    sum_oa += cell_overlap();
    sum_rn += 0.0;//rudy(netToCell);
  }

  wl.first = 0.0;
  wl.second = sum_wl / (float)initial_loop_iter;

  area.first = 0.0;
  area.second = sum_oa / (float)initial_loop_iter;
  area.second = max(area.second, 1.0); 

  rn.first = 0.0;
  rn.second = sum_rn / (float)initial_loop_iter;

  normalization_terms.push_back(wl);
  normalization_terms.push_back(area);
  normalization_terms.push_back(rn);

  wl_normalization = normalization_terms[0];
  area_normalization = normalization_terms[1];
  routability_normalization = normalization_terms[2];
}

void GridBasedPlacer::h_initialize_params(map<int, vector<Module*> > &netToCell) {

  vector < std::pair <double,double> > normalization_terms;

  int num_components = 0;
  vector < Module* > ::iterator itNode;
  for (itNode = moduleId.begin(); itNode != moduleId.end(); ++itNode) {
    if(!(*itNode) -> terminal) {
      num_components += 1;
    }
  }
  std::pair < double, double > wl;
  std::pair < double, double > area;
  std::pair < double, double > rn;

  double sum_wl = 0.0;
  double sum_oa = 0.0;
  double sum_rn = 0.0;

  for (int i = 0; i<initial_loop_iter; i++) {
    h_random_initial_placement();
    sum_wl += h_wirelength(netToCell);
    sum_oa += h_cell_overlap();
    sum_rn += 0.0;//rudy(netToCell);
  }

  wl.first = 0.0;
  wl.second = sum_wl / (float)initial_loop_iter;

  area.first = 0.0;
  area.second = sum_oa / (float)initial_loop_iter;
  area.second = max(area.second, 1.0); 

  rn.first = 0.0;
  rn.second = sum_rn / (float)initial_loop_iter;

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

  int ro = 0;
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
  validate_move(n, rx, ry);
}

/*
initial_placement
Randomly place and orient all movable components in the board area
*/
void GridBasedPlacer::random_initial_placement() {
  vector < Node > ::iterator itNode;
  for (itNode = nodeId.begin(); itNode != nodeId.end(); ++itNode) {
    if (!itNode -> fixed && !itNode -> terminal) {
      //random_placement(mMinX, mMaxX - max(itNode -> width, itNode -> height), mMinY, mMaxY - max(itNode -> width, itNode -> height), *itNode);
        random_placement(mMinX, mMaxX, mMinY, mMaxY, *itNode);
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
              validate_move(*nodeit1,rx + dx,ry);
            } else {
              // project in y direction
              validate_move(*nodeit1,rx,ry + dy);
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
    for(size_t j = i; j < nodeId.size(); j++) {
      if (i == j) {continue;}
      if((nodeId[i].fixed && nodeId[j].fixed) || (nodeId[i].layer != nodeId[j].layer)) {
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
      for ( rtree_t::const_query_iterator it = rtree.qbegin(index::intersects(nodeId[i].envelope)) ;
        it != rtree.qend() ; ++it ) {
        size_t j = it->second;
        if (cell_history.find(nodeId[j].idx) != cell_history.end()) {
          continue;
        }
        if((nodes[i]->fixed && nodeId[j].fixed)  || (nodeId[i].layer != nodeId[j].layer)) {
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
        if((nodes[i]->fixed && nodeId[j].fixed) || (nodeId[i].layer != nodeId[j].layer)) {
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

/*
cost
Compute full cost over all nodes
*/
double GridBasedPlacer::cost(map<int, vector<pPin> > &netToCell, int temp_debug) {
  double l2 = 1 - l1;

  double wl = wirelength(netToCell);
  double oa = cell_overlap();
  double rd = 0.0;//rudy(netToCell);
  double normalized_wl = (wl - wl_normalization.first)/(wl_normalization.second - wl_normalization.first);
  double normalized_oa = (oa - area_normalization.first)/(area_normalization.second - area_normalization.first);
  double normalized_rd = 0.0;//(rd - routability_normalization.first)/(routability_normalization.second - routability_normalization.first);
  wl_hist.push_back(normalized_wl);
  oa_hist.push_back(normalized_oa);

  double wirelength_cost = l1 * normalized_wl;
  double overlap_cost = l2 * 0.85 * normalized_oa;
  double routability_cost = l2 * 0.15 * normalized_rd;
  double total_cost = wirelength_cost + overlap_cost + routability_cost;
  cout << "cost: " << total_cost << " wirelength: " << wl << " " << wirelength_cost << " overlap: " << oa << " " << overlap_cost << endl;

  if (oa < 1 && wl < best_wl) {
      bestSol = nodeId;
      best_wl = wl;
      best_overlap = oa;
  } else if (oa < best_overlap) {
     bestSol = nodeId;
     best_wl = wl;
     best_overlap = oa;
  }
  
  return total_cost;
}

/*
cost_partial
compute partial cost over subset of nodes
*/
double GridBasedPlacer::cost_partial(vector < Node *> &nodes, map<int, vector<pPin> > &netToCell) {
  double l2 = 1-l1;
  return l1 * (wirelength_partial(nodes, netToCell) - wl_normalization.first)/(wl_normalization.second - wl_normalization.first) +
         l2 * 0.85 * (cell_overlap_partial(nodes) - area_normalization.first)/(area_normalization.second - area_normalization.first) +
         l2 * 0.15 * 0.0;//(rudy(netToCell) - routability_normalization.first)/(routability_normalization.second - routability_normalization.first);

}

/*
random_node
Select random node from set of nodes - returns iterator
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
  double validated_rx = max(rx, mMinX);
  double validated_ry = max(ry, mMinY);

  validated_rx = min(validated_rx, mMaxX - node.width);
  validated_ry = min(validated_ry, mMaxY - node.height);

  node.setPos(validated_rx,validated_ry);
}

/*
initiate_move
initiate a single move for annealing
*/
//double c = 0.0;
double GridBasedPlacer::initiate_move(double current_cost, map<int, vector<pPin> > &netToCell) {
  int state = -1;
  double prevCost = 0.0;
  if(debug > 1) {
    cout << "current_cost: " << current_cost << endl;
  }
  vector < Node > ::iterator rand_node1;
  vector < Node > ::iterator rand_node2;

  vector < Node* > perturbed_nodes;

  int r = 0;
  rand_node1 = random_node();
  while(rand_node1->terminal || rand_node1->name == "" || rand_node1->fixed) {
    rand_node1 = random_node();
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
  if (i < swap_proba) { // swap
    state = 0;
    if(debug > 1) {
      cout << "swap" << endl;
    }
    rand_node2 = random_node();
    while(rand_node2->terminal || rand_node2->idx == rand_node1->idx || rand_node2->name == "" || rand_node2->fixed) {
      rand_node2 = random_node();
    }
    perturbed_nodes.push_back(&(*rand_node2));
    prevCost = cost_partial(perturbed_nodes,netToCell);

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
  } else if (swap_proba <= i && i < shift_proba + swap_proba) { // shift
    state = 1;
    if(debug > 1) {
      cout << "shift" << endl;
    }
    prevCost = cost_partial(perturbed_nodes, netToCell);

    double sigma = rand_node1->sigma;
    ssamp = sigma;
    boost::normal_distribution<> nd(sigma, shift_var);
    boost::variate_generator<boost::mt19937&,
                             boost::normal_distribution<> > var_nor(rng, nd);

    double dx = var_nor();
    double dy = var_nor();
    double rx = rand_node1_orig_x + dx;
    double ry = rand_node1_orig_y + dy;

    validate_move(*rand_node1, rx, ry);
    if(rt) {
      rtree.insert(std::make_pair(rand_node1->envelope, rand_node1->idx));
    }
  } else if(shift_proba + swap_proba <= i && i < shift_proba + swap_proba + rotate_proba) { // rotate
    state = 2;
    if(debug > 1) {
      cout << "rotate" << endl;
    }
    prevCost = cost_partial(perturbed_nodes, netToCell);

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
    validate_move(*rand_node1, rand_node1_orig_x, rand_node1_orig_y);

    if(rt) {
      rtree.insert(std::make_pair(rand_node1->envelope, rand_node1->idx));
    }
  } else { // layer change
    state = 3;
    if(debug > 1) {
      cout << "layer change" << endl;
    }
    prevCost = cost_partial(perturbed_nodes, netToCell);
    rand_node1->layerChange();
  }
  double transition_cost = cost_partial(perturbed_nodes,netToCell);
  double updated_cost = current_cost - prevCost + transition_cost;

  bool accept = check_move(current_cost, updated_cost);

  if (!accept) {
    AcceptRate = 1.0/500.0 *(499.0*AcceptRate);
    if(debug > 1) { 
      cout << "reject" << endl;
    }
    // revert state if reject
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
    } else if (state == 3) {
      rand_node1->layerChange();
    }
    if(rt) {
      rtree.insert(std::make_pair(rand_node1->envelope, rand_node1->idx));
    }
    return current_cost;
  } else {
    if(debug > 1) {
      cout << "accept" << endl;
    }
    AcceptRate = 1.0/500.0 *(499.0*AcceptRate + 1.0);
    return updated_cost;
  }
}

/*
update_Temperature
Update the SA parameters according to annealing schedule
*/
void GridBasedPlacer::update_temperature() {
  vector < Node > ::iterator nodeit = nodeId.begin();
  l1 = 0.98*l1;
  if (Temperature > 50e-3) {
    Temperature = (0.985) * Temperature;
    for (nodeit = nodeId.begin(); nodeit != nodeId.end(); ++nodeit) {
      nodeit->sigma =  max(0.985*nodeit->sigma,3/4 * nodeit->sigma);
    }
  } else if (Temperature > 10e-3) {
    Temperature = (0.9992) * Temperature;
    for (nodeit = nodeId.begin(); nodeit != nodeId.end(); ++nodeit) {
      nodeit->sigma =  max(0.9992*nodeit->sigma,2/4 * nodeit->sigma);
    }
  } else if (Temperature > 50e-4) {
    Temperature = (0.9955) * Temperature;
    for (nodeit = nodeId.begin(); nodeit != nodeId.end(); ++nodeit) {
      nodeit->sigma =  max(0.9955*nodeit->sigma,1/4 * nodeit->sigma);
    }
  } else if (Temperature > 10e-4) {
    Temperature = (0.9965) * Temperature;
    for (nodeit = nodeId.begin(); nodeit != nodeId.end(); ++nodeit) {
      nodeit->sigma =  max(0.9965*nodeit->sigma,2.0);
    }
  } else {
    if (Temperature > 10e-8) {
      Temperature = (0.885) * Temperature;
    } else {
      Temperature = 0.0000000001;
    }
    for (nodeit = nodeId.begin(); nodeit != nodeId.end(); ++nodeit) {
      nodeit->sigma =  max(0.855*nodeit->sigma,1.0);
    }
  }
}

/*
modified_lam_update
Update the SA parameters according to modified lam schedule
*/
void GridBasedPlacer::modified_lam_update(int i) {
  if(do_hplace){
    int outer_loop_iter = (H.levels.size()-1) * outer_loop_iter;
  }

  if ((double)i/(double)outer_loop_iter < 0.15){
    LamRate = 0.44 + 0.56 * pow(560.0, -i/outer_loop_iter/0.15);
  } else if (0.15 <= (double)i/(double)outer_loop_iter && (double)i/(double)outer_loop_iter <= 0.65) {
    LamRate = 0.44;
    //wl_normalization.second = wl_hist.back();
    //area_normalization.second = oa_hist.back();
  } else if (0.65 <= (double)i/(double)outer_loop_iter) {
    LamRate = 0.44 * pow(440.0, -((double)i/(double)outer_loop_iter - 0.65)/0.35);
    //wl_normalization.second = wl_hist.back();
    //area_normalization.second = oa_hist.back();
  }

  if (AcceptRate > LamRate) {
    Temperature = Temperature * lamtemp_update;
    l1 = 0.95*l1;
  } else {
    Temperature = min(Temperature / lamtemp_update, 1.2);
    l1 = 0.96*l1;
  }

  if (!do_hplace) {
    vector < Node > ::iterator nodeit = nodeId.begin();
    for (nodeit = nodeId.begin(); nodeit != nodeId.end(); ++nodeit) {
      nodeit->sigma =  max(0.98*nodeit->sigma,0.5);
    }
  } else {
    vector < Module* > ::iterator nodeit = moduleId.begin();
    for (nodeit = moduleId.begin(); nodeit != moduleId.end(); ++nodeit) {
      (*nodeit) -> sigma =  0.2;
    }
  }
}

/*
check_move
either accept or reject the move based on current & previous temperature & cost
*/
bool GridBasedPlacer::check_move(double prevCost, double newCost) {
  double delCost = 0;
  boost::uniform_real<> uni_dist(0,1);
  boost::variate_generator<boost::mt19937&, boost::uniform_real<> > uni(rng, uni_dist);
  double prob = uni();
  if(debug > 1) {
    cout << "new cost: " << newCost << endl;
  }
  delCost = newCost - prevCost;
  //cout << delCost << " " << Temperature << " " << prob << " " << exp(-delCost/Temperature) << " " << -delCost << " " << Temperature << " " << -delCost/Temperature << " " << prob << " " <<(prob <= (exp(-delCost/Temperature))) << endl;
  if (delCost <= 0 || prob <= (exp(-delCost/Temperature))) {
    prevCost = newCost;
    return true;
  } else {
    return false;
  }
}

/*
initialize_temperature
initialize temperature according to varanelli cahoon
*/
double GridBasedPlacer::initialize_temperature(map<int, vector<pPin> > &netToCell) {
  double t = 0.0;
  double emax = 0.0;
  double emin = 0.0;
  double xt = 1.0;
  double x0 = 0.84;
  double p = 2.0;
  random_initial_placement();
  for(int i=1; i<=10; i++){
    for(int j=1; j<=10; j++){
      random_initial_placement();
      emax += exp(cost(netToCell)/t);
      initiate_move(0.0, netToCell);
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
                                 map<int, vector<pPin> > &netToCell) {
    vector < double > cost_hist = report["cost_hist"];
    vector < double > wl_hist   = report["wl_hist"];
    vector < double > oa_hist   = report["oa_hist"];

    time_t t = time(0);
    struct tm * now = localtime( & t );

    char buf [80];
    strftime (buf,80,"%Y-%m-%d-%H-%M-%S",now);
    string buffer = string(buf);

    double wl = wirelength(netToCell);
    double oa = cell_overlap();
    double routability = rudy(netToCell);
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
}

/*
annealer
main loop for sa algorithm
*/
float GridBasedPlacer::annealer(map<int, vector<pPin> > &netToCell, string initial_pl) {
  Temperature = t_0;
  int num_components = 0;

  vector < Node > ::iterator itNode;
  map < string, vector < double > > report;
  report["cost_hist"] = cost_hist;
  report["wl_hist"] = wl_hist;
  report["oa_hist"] = oa_hist;

  cout << "calculating initial params..." << endl;
  initialize_params(netToCell);
  if (initial_pl != "") {
    readPlFile(initial_pl);
  } else {
    random_initial_placement();
  }

  cout << "calculating initial cost estimate..." << endl;  
  double cst = cost(netToCell,-1);
  if(var) {
    Temperature = initialize_temperature(netToCell);
  }
  int idx = 0;

  if(rt) {
    rtree.clear();
    for (itNode = nodeId.begin(); itNode != nodeId.end(); ++itNode) {
      rtree.insert(std::make_pair(itNode -> envelope, idx));
      num_components += 1;
      idx+=1;
    }
  } else {
    for (itNode = nodeId.begin(); itNode != nodeId.end(); ++itNode) { 
      num_components += 1;
      idx+=1;
    }
  }

  high_resolution_clock::time_point t1 = high_resolution_clock::now();
  long long int ii = 0; // outer loop iterator
  int i = 0; // inner loop iterator
  cout << "beginning optimization..." << endl;
  while (ii < outer_loop_iter) {
    i = inner_loop_iter*num_components; 

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast< duration<double> >(t2 - t1);

    cout << "=====" << ii << "=====" << endl;
    cout << "iteration: " << ii << " time: " <<  time_span.count() << " (s)" << " updates/time: " <<  ii/time_span.count() << 
    " time remaining: " <<  time_span.count()/ii * (outer_loop_iter-ii) << " (s)" << " temperature: " << Temperature << " wl weight: " << l1 << " s samp: " << ssamp <<
    " acceptance rate: " << AcceptRate << " lam rate: " << LamRate << endl;

    //gen_report(report,
    //           accept_ratio_history,
    //           netToCell);

      if (ii % 10 == 0) {
          cst = cost(netToCell);
      }

    while (i > 0) {
      cst = initiate_move(cst, netToCell);
      cost_hist.push_back(cst);
      i -= 1;
    }

    // convergence criterion
    if(eps > 0 && abs(cost_hist.end()[-1] - cost_hist.end()[-2]) < eps) {
      break;
    }
    if (lam) {
      modified_lam_update(ii);
    } else {
      update_temperature();
    }
    ii += 1;
  }
  vector < Node > ::iterator nodeit;
  for (nodeit = nodeId.begin(); nodeit != nodeId.end(); ++nodeit) {
    nodeit->sigma =  0.2;
  }

  i = 5 * inner_loop_iter * num_components;
  Temperature = 0.0;
  cst = cost(netToCell);
  while (i > 0) {
    cst = initiate_move(cst, netToCell);
    if (i%10==0) {
        cst = cost(netToCell);
    }
    i -= 1;
  }

  return cost(netToCell);
}

/*
hannealer
main loop for hierarchical sa algorithm
*//*
float GridBasedPlacer::hannealer(map<int, vector<Module *> > &netToCell, string initial_pl, int level) {
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

  int idx = 0;

  for (itNode = moduleId.begin(); itNode != moduleId.end(); ++itNode) {
    idx+=1;
    if(!(*itNode) -> terminal && !(*itNode) -> fixed) {
      num_components += 1;
    }
  }

  high_resolution_clock::time_point t1 = high_resolution_clock::now();
  long long int ii = 0; // outer loop iterator
  int i = 0; // inner loop iterator
  cout << "beginning optimization..." << endl;
  while (ii < outer_loop_iter) {
    i = inner_loop_iter*num_components; 

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast< duration<double> >(t2 - t1);

    cout << "=====" << ii << "=====" << endl;
    cout << "iteration: " << ii << " time: " <<  time_span.count() << " (s)" << " updates/time: " <<  ii/time_span.count() << 
    " time remaining: " <<  time_span.count()/ii * (outer_loop_iter-ii) << " (s)" << " temperature: " << Temperature << " wl weight: " << l1 << " s samp: " << ssamp <<
    " acceptance rate: " << AcceptRate << " lam rate: " << LamRate << endl;

    if (ii % 10 == 0) {
        cst = h_cost(netToCell);
        nodeId = H.update_cell_positions_at_level(nodeId, level);
        writePlFile("./cache/"+std::to_string( level )+"_"+std::to_string( ii )+".pl");
    }
    while (i > 0) {
      cst = h_initiate_move(cst, netToCell);
      cost_hist.push_back(cst);
      i -= 1;
    }
    if (lam) {
      modified_lam_update(ii);
    } else {
      update_temperature();
    }
    ii += 1;
  }

  return h_cost(netToCell);
}
*/


//////////////////////////////////////

/*
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

/*
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

/*
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

/*
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

/*
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

/*
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

/*
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

/*
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

  double wl = h_wirelength(netToCell);
  double oa = 0.0;
  double rd = 0.0;//rudy(netToCell);

  if (densityFlag == 0) {
    oa = h_cell_overlap();
  } else {
    oa = h_cellDensity();
  }

  double normalized_wl = (wl - wl_normalization.first)/(wl_normalization.second - wl_normalization.first);
  double normalized_oa = (oa - area_normalization.first)/(area_normalization.second - area_normalization.first);
  double normalized_rd = 0.0;//(rd - routability_normalization.first)/(routability_normalization.second - routability_normalization.first);
  wl_hist.push_back(normalized_wl);
  oa_hist.push_back(normalized_oa);

  double wirelength_cost = l1 * normalized_wl;
  double overlap_cost = l2 * 0.85 * normalized_oa;
  double routability_cost = l2 * 0.15 * normalized_rd;

  //double routability_cost = l1 * 0.1 * (this->rudy(netToCell) - routability_normalization.first)/(routability_normalization.second - routability_normalization.first);
  double total_cost = wirelength_cost + overlap_cost;// + routability_cost;
  cout << "cost: " << total_cost << " wirelength: " << wl << " " << wirelength_cost << " overlap: " << oa << " " << overlap_cost << endl;
  
  if (oa < 1 && wl < best_wl) {
      bestSol = nodeId;
      best_wl = wl;
      best_overlap = oa;
  } else if (oa < best_overlap) {
     bestSol = nodeId;
     best_wl = wl;
     best_overlap = oa;
  }

  return total_cost;
}

double GridBasedPlacer::h_cost_partial(vector < Module *> &nodes, map<int, vector<Module *> > &netToCell) {
  double l2 = 1-l1;
  if(densityFlag == 0) {
    return l1 * (h_wirelength_partial(nodes, netToCell) - wl_normalization.first)/(wl_normalization.second - wl_normalization.first) +
           l2 * 0.85 * (h_cell_overlap_partial(nodes) - area_normalization.first)/(area_normalization.second - area_normalization.first) +
           l2 * 0.15 * 0.0;//(rudy(netToCell) - routability_normalization.first)/(routability_normalization.second - routability_normalization.first);

  } else {
    return l1 * (h_wirelength_partial(nodes, netToCell) - wl_normalization.first)/(wl_normalization.second - wl_normalization.first) +
           l2 * 0.85 * (h_cellDensity() - area_normalization.first)/(area_normalization.second - area_normalization.first) +
           l2 * 0.15 * 0.0;//(rudy(netToCell) - routability_normalization.first)/(routability_normalization.second - routability_normalization.first);
  }
}

/*
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

/*
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

//double c = 0.0;
double GridBasedPlacer::h_initiate_move(double current_cost, map<int, vector<Module *> > &netToCell) {
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

  if (i< swap_proba) { // swap
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
  } else if (swap_proba <= i && i < shift_proba + swap_proba) { // shift
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
  } else if (shift_proba + swap_proba <= i && i < shift_proba + swap_proba + rotate_proba) { // rotate
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

  bool accept = this->h_check_move(current_cost, updated_cost, Temperature);
  if (!accept) {
    AcceptRate = 1.0/500.0 *(499.0*AcceptRate);
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
    AcceptRate = 1.0/500.0 *(499.0*AcceptRate + 1.0);
    if(debug > 1) {
      cout << "accept" << endl;
    }
    return updated_cost;
  }
}

/*
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
      this->h_initiate_move(0.0, netToCell);
      emin += exp(h_cost(netToCell)/t);
    }
    xt = emax/emin;
    t = t * pow(log(xt),1/p)/log(x0);
  }
  return t;
}

/*
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
    cout << (*itNode) -> width << " " << (*itNode) -> height << endl;
  }

  high_resolution_clock::time_point t1 = high_resolution_clock::now();
  long long int ii = 0; // outer loop iterator
  int i = 0; // inner loop iterator
  cout << "beginning iterations..." << endl;
  while (ii < outer_loop_iter) {
    i = inner_loop_iter*num_components; 
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast< duration<double> >(t2 - t1);

    wl_hist.push_back((h_wirelength(netToCell) - wl_normalization.first)/(wl_normalization.second - wl_normalization.first));
    if (densityFlag == 0) {
      oa_hist.push_back((h_cell_overlap() - area_normalization.first)/(area_normalization.second - area_normalization.first));
    } else {
      oa_hist.push_back(h_cellDensity());
    }
    l_hist.push_back(l1);
    temp_hist.push_back(Temperature);

    cout << "=====" << ii << "=====" << endl;
    cout << "iteration: " << ii << " time: " <<  time_span.count() << " (s)" << " updates/time: " <<  ii/time_span.count() << 
    " time remaining: " <<  time_span.count()/ii * (outer_loop_iter-ii) << " (s)" << " temperature: " << Temperature << " wl weight: " << l1 << " s samp: " << ssamp <<
    " acceptance rate: " << AcceptRate << " lam rate: " << LamRate << endl;

    nodeId = H.update_cell_positions_at_level(nodeId, level);
    writePlFile("./cache/"+std::to_string( level )+"_"+std::to_string( ii )+".pl");

    if (ii % 10 == 0) {
        cst = h_cost(netToCell,-1);
    }
    while (i > 0) {
      cst = this->h_initiate_move(cst, netToCell);
      report["cost_hist"].push_back(cst);
      i -= 1;
    }

    // convergence criterion
    if(eps > 0 && abs(cost_hist.end()[-1] - cost_hist.end()[-2]) < eps) {
      break;
    }
    if (lam) {
      modified_lam_update(level*outer_loop_iter + ii);
    } else {
      update_temperature();
    }
    ii += 1;
  }
//////////
  for (itNode = moduleId.begin(); itNode != moduleId.end(); ++itNode) {
    (*itNode) -> sigma =  0.2;
  }

  i = 5 * inner_loop_iter * num_components;
  Temperature = 0.0;
  cst = h_cost(netToCell);
  while (i > 0) {
    cst = h_initiate_move(cst, netToCell);
    if (i%10==0) {
        cst = h_cost(netToCell);
    }
    i -= 1;
  }
  return this->h_cost(netToCell);
}