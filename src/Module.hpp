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

#include <string>

#include <boost/regex.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim_all.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/tuple/tuple.hpp>

using namespace std;
namespace trans = boost::geometry::strategy::transform;
namespace bg = boost::geometry;
typedef boost::geometry::model::d2::point_xy<double> Point;

class Module {
  public:
    string name;
    int level;
    int idx;
    model::polygon<model::d2::point_xy<double> > poly;
    boost::geometry::model::box< model::d2::point_xy<double> > envelope;
    double width;
    double height;
    double sigma;
    double xCoordinate;
    double yCoordinate;
    double xBy2;
    double yBy2;
    double initialX;
    double initialY;
    int orientation = 0;

    double x_offset = 0.0;
    double y_offset = 0.0; // virtual pin offsets


    vector < int > Netlist; // vector of netlist ids this module is associated with [for online wl comp.]
    vector < int > cells;

    vector < Module *> children;
    Module *parent = nullptr;
    bool leaf = false;
    bool root = false;
    bool macroModule = false;
    int fixed = 0;


  void setNetList(int NetId) {
    // Sets parameters given an entry in Nets file
    Netlist.push_back(NetId);
  }

  void setParameterNodes(double width, double height) {
    // Sets parameters given an entry in Nodes file
    this -> xCoordinate = 0.0;
    this -> yCoordinate = 0.0;
    this -> width = width;
    this -> height = height;
    double points[][2] = {{0.0, 0.0}, {width, 0.0}, {width, height}, {0.0, height}};
    model::polygon< model::d2::point_xy<double> > poly;
    append(poly, points);
    boost::geometry::correct(poly);
    this -> poly = poly;
    boost::geometry::envelope(poly, envelope);
  }

  void setParameterPl(double xCoordinate, double yCoordinate) {
    // Sets parameters given an entry in Pl file
    this -> setPos(xCoordinate, yCoordinate);
    this -> sigma = 50.0;
  }

  void setPos(double x, double y) {
    // Sets position of a node object (lower-left corner)
    trans::translate_transformer<double, 2, 2> translate(x - this->xCoordinate, y - this->yCoordinate);
    model::polygon<model::d2::point_xy<double> > tmp;
    boost::geometry::transform(this->poly, tmp, translate);
    this->poly = tmp;
    updateCoordinates();
  }

  /*
  upateCoordinates
  Updates parameters of Node class from a geometry object
  */
  void updateCoordinates() {
    boost::geometry::model::d2::point_xy<double> centroid;
    boost::geometry::centroid(this->poly, centroid);
    this -> xBy2 = centroid.get<0>();
    this -> yBy2 = centroid.get<1>();

    boost::geometry::model::box< model::d2::point_xy<double> > envtmp;
    boost::geometry::envelope(this->poly, envtmp);
    this->envelope = envtmp;
    Point minCorner = this->envelope.min_corner();
    Point maxCorner = this->envelope.max_corner();
    this->width  = maxCorner.get<0>() - minCorner.get<0>();
    this->height = maxCorner.get<1>() - minCorner.get<1>();
    this->xCoordinate = bg::get<bg::min_corner, 0>(this->envelope);
    this->yCoordinate = bg::get<bg::min_corner, 1>(this->envelope);
  }

  /*
  printExterior
  print polygon vertices
  */
  void printExterior() const{
    for(auto it = boost::begin(boost::geometry::exterior_ring(this->poly)); it != boost::end(boost::geometry::exterior_ring(this->poly)); ++it) {
        double x = bg::get<0>(*it);
        double y = bg::get<1>(*it);
        cout << x << " " << y << endl;
    }
  }

  /*
  printParameter
  print node params
  */
  void printParameter() {
    cout << "name      " << this -> name << endl;
    cout << "Width         " << this -> width << endl;
    cout << "Height        " << this -> height << endl;
    cout << "X_Co-ordinate " << this -> xCoordinate << endl;
    cout << "Y_Co-ordinate " << this -> yCoordinate << endl;
    cout << "X/2           " << xBy2 << endl;
    cout << "Y/2           " << yBy2 << endl;
    cout << "NetList       ";
    vector < int > ::iterator it2;
    for (it2 = Netlist.begin(); it2 != Netlist.end(); ++it2) {
      cout << * it2 << " ";
    }
    this -> printExterior();
    cout << "\n" << endl;
  }

  void init_module(int id, int lev, bool r) {
    root = r;
    idx = id;
    level = lev;
  }
  void add_child(Module* m) {
    children.push_back(m);
  }
  void insert_cell(int i) {
    cells.push_back(i);
  }
  bool isleaf() {
    return leaf;
  }
};

class Level {
public:
  vector <Module *> modules; // mapping from module id to module pointer
  map<int, vector<Module *> > netToModule; // mapping from netid to vpin/module
  int level;
};

class Hierarchy {
public:
  int num_levels;
  vector < int > num_modules_per_layer;
  //vector < vector < Module * > > levels;
  vector < Level > levels;
  Module *root;
  vector < Node > id2cell;

  void init_hierarchy(int nl, vector <int> nm) {
    num_levels = nl;
    num_modules_per_layer = nm;

    for (int l=0; l<=num_levels; ++l) {
      Level tmp;
      tmp.level = l;
      levels.push_back(tmp);
    }
    root = create_hierarchy(0,0,true);
    root->root = true;
  }

  void insert_cell(int cell_id, vector < int> cluster_id_vec, Module *m, int level) {
    if(num_levels <= level) {
      return;
    }
    if(cluster_id_vec.size() < 1) {
      // macro cell - just insert into every level
      vector < Module * > macroModulePath;
      for(int l=0; l<=num_levels; ++l) {
        Module *tmp = new Module;
        tmp->init_module(-cell_id, l, false);
        tmp->macroModule = true;
        tmp->fixed = 0;
        tmp->insert_cell(cell_id);
        levels[l].modules.push_back(tmp);
        macroModulePath.push_back(tmp);
      }
      for(int l=0; l<num_levels; ++l) {
        Module *tmp = macroModulePath[l];
        if (l == 0) {
          root->children.push_back(tmp);
          tmp->parent = root;
        } else {
          tmp->parent = macroModulePath[l-1];
        }
        tmp->children.push_back(macroModulePath[l+1]);
      }
      macroModulePath[num_levels]->parent = macroModulePath[num_levels-1];
      macroModulePath[num_levels]->leaf = true;
    } else {
      int cluster_id = cluster_id_vec[level];
      m->children[cluster_id-1]->insert_cell(cell_id);
      insert_cell(cell_id, cluster_id_vec, m->children[cluster_id-1], level+1);
    }
  }

  Module *get_leaf_module_from_id(int idx, Module *m) {
    if(m->leaf) {
      return m;
    } else {
      for (auto &child : m->children) {
        if(std::find(child->cells.begin(), child->cells.end(), idx) != child->cells.end()) {
          return get_leaf_module_from_id(idx, child);
        }
      }
    }
  }

  /*
  for each level in the hierarchy, instantiate a netlist for that level:
  associated nets for each module (list of net ids)
  netlist-module mapping id -> list of modules
  */
  void set_netlist_hierarchy(map<int, vector<Pin> > netToCell) {
    int netidx = 1;
    for (auto &net : netToCell) { 
      int netid = net.first;
      vector<Pin> pvec = net.second; //pvec_i.idx -> cell id
      vector<Module *> ms; // convert vector of pins to vector of clusters for that net add to netToModule in leaf-level

      for (auto &pin : pvec) {
        int cell_id = pin.idx;
        Module *m = get_leaf_module_from_id(cell_id, root);

        if(!m->leaf || !(std::find(m->cells.begin(), m->cells.end(), cell_id) != m->cells.end())) {
          cout << "cell not found in leaf nodes: " << cell_id << " " << m->leaf << endl; // COULD BE MACRO
          return;
        }

        if(!(std::find(ms.begin(), ms.end(), m) != ms.end())) { // uniqueify module-level netlist
          m->setNetList(netidx);
          ms.push_back(m);
        }

      }
      levels.back().netToModule.insert(pair < int, vector < Module * > > (netidx, ms));
      netidx ++;
    }

    // NOW need to propagate up the hierarchy
    cout << "propagating up.." << endl;
    propagate_netlist(levels.back().level);
  }

  void propagate_netlist(int lev) {
    if (lev <= 0) {
      return;
    }

    map<int, vector< Module * > >netToModules = levels[lev].netToModule;
    int netidx = 1;
    for (auto &net : netToModules) { 
      int netid = net.first;
      vector<Module *> pvec = net.second; //pvec_i.idx -> cell id
      vector<Module *> ms; // convert vector of pins to vector of clusters for that net add to netToModule in leaf-level

      for (auto &module : pvec) {
        Module *m = module->parent;

        if(!(std::find(ms.begin(), ms.end(), m) != ms.end())) { // uniqueify module-level netlist
          m->setNetList(netidx); // crash here, module is root for some reason so parent m is garbage
          ms.push_back(m);
        }
      }
      if (ms.size() > 0) {
        levels[lev-1].netToModule.insert(pair < int, vector < Module* > > (netidx, ms));
      }
      netidx ++;
    }
    propagate_netlist(lev - 1);
  }

  void set_module_geometries(vector < Node > nodeId) { // should fix to be bottom up
    // first get leaf params
    for (auto &m : levels.back().modules) {
      double xmin = 99999999999999.0;
      double ymin = 99999999999999.0;
      double xmax = 0.0;
      double ymax = 0.0;
      for (auto &cellid : m->cells) {
        double cellx = nodeId[cellid].xCoordinate;
        double celly = nodeId[cellid].yCoordinate;

        xmin = min(xmin, cellx);
        ymin = min(ymin, celly);
        xmax = max(xmax, cellx);
        ymax = max(ymax, celly);
      }
      m->setParameterNodes(xmax - xmin, ymax - ymin);
      m->setParameterPl(xmin, ymin);
    }
    propagate_geometries(num_levels-1);
  }

  void propagate_geometries(int lev) {
    if (lev <= 0) {
      return;
    }

    vector <Module *> modules = levels[lev].modules;
    for (auto &module : modules) {
      double xmin = 99999999999999.0;
      double ymin = 99999999999999.0;
      double xmax = 0.0;
      double ymax = 0.0;
      for (auto &m : module->children) {
        double cellx = m->xCoordinate;
        double celly = m->yCoordinate;

        xmin = min(xmin, cellx);
        ymin = min(ymin, celly);
        xmax = max(xmax, cellx);
        ymax = max(ymax, celly);
      }
      module->setParameterNodes(xmax - xmin, ymax - ymin);
      module->setParameterPl(xmin, ymin);      
    }
    propagate_geometries(lev - 1);
  }

  Module* create_hierarchy(int id, int lev, bool r) {
    Module *m = new Module;
    m->init_module(id, lev, r);
    levels[lev].modules.push_back(m);

    if (num_levels <= lev) {
      m->leaf = true;
      return m;
    }

    int num_children = num_modules_per_layer[lev];
    for (int i =1; i<=num_children; ++i) {
      Module* mc;
      mc = create_hierarchy(i, lev+1, false);
      mc->parent = m;
      if(mc) {
        m->add_child(mc);
      }
    }
    
    return m;
  }
};
