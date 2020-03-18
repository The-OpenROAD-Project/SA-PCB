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

//#include "Node.hpp"
//#include "Pin.hpp"
#include "Module.hpp"

using namespace std;
using namespace boost::geometry;
namespace trans = boost::geometry::strategy::transform;
namespace bg = boost::geometry;
typedef boost::geometry::model::d2::point_xy<double> Point;

/**
Level
contains an instance/level of the hierarchal board representation to place
*/
class Level {
public:
vector <Module *> modules; // mapping from module id to module pointer
map<int, vector<Module *> > netToModule; // mapping from netid to vpin/module
int level;
};

/**
Hierarchy
Hierarchical tree datastructure
*/
class Hierarchy {
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
};
