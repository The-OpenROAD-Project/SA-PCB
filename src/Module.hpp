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
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctype.h>
#include <map>
#include <vector>


#include <boost/regex.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim_all.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/tuple/tuple.hpp>

#include "Node.hpp"
#include "Pin.hpp"

using namespace std;
using namespace boost::geometry;
namespace trans = boost::geometry::strategy::transform;
namespace bg = boost::geometry;
typedef boost::geometry::model::d2::point_xy<double> Point;

/**
Module
Base module class for hierarchy
*/
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

    int terminal = 0;

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

    void setNetList(int NetId);
    void setParameterNodes(double width, double height);
    void setParameterPl(double xCoordinate, double yCoordinate);
    void setPos(double x, double y);
    void updateCoordinates();
    void printExterior() const;
    void printParameter();
    void init_module(int id, int lev, bool r);
    void add_child(Module* m);
    void insert_cell(int i);
    bool isleaf();
};

/**
Level
contains an instance/level of the hierarchal board representation to place
*/
/*class Level {
public:
vector <Module *> modules; // mapping from module id to module pointer
map<int, vector<Module *> > netToModule; // mapping from netid to vpin/module
int level;
};
*/
/**
Hierarchy
Hierarchical tree datastructure
*/
/*class Hierarchy {
public:
int num_levels;
vector < int > num_modules_per_layer;
//vector < vector < Module * > > levels;
vector < Level > levels;
Module *root;
//vector < Node > id2cell;

void init_hierarchy(int nl, vector <int> nm);
Module* create_hierarchy(int id, int lev, bool r);
void print_param(int i);
void insert_cell(int cell_id, vector < int> cluster_id_vec, Module *m, int level);
Module* get_leaf_module_from_id(int idx, Module *m);
Module* get_level_module_from_id(int idx, int level, Module *m);
void set_netlist_hierarchy(map<int, vector<pPin> > netToCell);
void propagate_netlist(int lev);
void set_module_geometries(vector < Node > nodeId);
void propagate_geometries(int lev);
vector <Node> update_cell_positions_at_level(vector < Node > nodeId, int level);
vector <Node> update_cell_positions(vector < Node > nodeId);
};*/
