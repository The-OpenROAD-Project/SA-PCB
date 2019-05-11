#include <iostream>
#include <vector>
#include <fstream>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/adapted/boost_tuple.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/adapted/c_array.hpp>
#include <string>
#include <iostream>
#include <regex>
#include <map>

using namespace std;
using namespace boost::geometry;
namespace trans = boost::geometry::strategy::transform;
namespace bg = boost::geometry;
typedef boost::geometry::model::d2::point_xy<double> point_type;

void readNodesFile(string fname);
void readWtsFile(string fname);
void readPlFile(string fname);
void readNetsFile(string fname);
void writePlFile(string fname);
void printMap();

class Node;
class Pin;
extern map < string, Node > nodeId;
//extern map<int, vector<string> > netToCell; // TODO netToNode consider going net -> node reference
extern map<int, vector<Pin> > netToCell;

class Pin {
  public:
  string name;
  double x_offset;
  double y_offset;

  void set_params(string nn, double xoffset, double yoffset) {
    this -> name = nn;
    this -> x_offset = xoffset;
    this -> y_offset = yoffset;
  }
};

class Node {
  public:
  string name;
  model::polygon<model::d2::point_xy<double> > poly;
  double width;
  double height;
  int weight;
  bool terminal;
  bool fixed;
  bool overlap;
  double sigma;
  double xCoordinate;
  double yCoordinate;
  int xBy2;
  int yBy2; // change to floats
  string orientation_str;
  int init_orientation;
  int orientation;
  vector < int > Netlist;

  void setParameterNodes(string name, double width, double height, int terminal) {
    // Sets parameters given an entry in Nodes file
    this -> name = name;
    this -> width = width;
    this -> height = height;
    this -> terminal = terminal;
    this->orientation_str = "N";
    if (terminal == 1) {
      double points[][2] = {{0.0, 0.0}, {width, 0.0}, {width, height}, {0.0, height}};
      model::polygon< model::d2::point_xy<double> > poly;
      append(poly, points);
      boost::geometry::correct(poly);
      this -> poly = poly;
    }
  }

  void setParameterWts(int weight) {
    // Sets parameters given an entry in weights file.
    // weights can correspond to size, importance, etc.
    this -> weight = weight;
  }

  void setParameterPl(int xCoordinate, int yCoordinate, string orientation_str, int fixed) {
    // Sets parameters given an entry in Pl file
    this -> xCoordinate = xCoordinate;
    this -> yCoordinate = yCoordinate;
    this -> orientation_str = orientation_str;
    this -> init_orientation = str2orient(orientation_str);
    this -> orientation = this -> init_orientation;
    this -> fixed = fixed;
    this -> sigma = 10.0;
  }

  void setNetList(int NetId) {
    // Sets parameters given an entry in Nets file
    Netlist.push_back(NetId);
  }

  void setPos(int x, int y) {
    // Sets position of a node object (lower-left corner)
    trans::translate_transformer<double, 2, 2> translate(x - this->xCoordinate, y - this->yCoordinate);
    model::polygon<model::d2::point_xy<double> > tmp;
    boost::geometry::transform(this->poly, tmp, translate);
    this->poly = tmp;
    updateCoordinates();
  }

  int str2orient(string o) {
    // convert a string to an orientation int
    if(o == "N") {
      return 0;
    } else if(o == "E") {
      return 1;
    } else if(o == "S") {
      return 2;
    } else if(o == "W") {
      return 3;
    }
    return -1;
  }

  string orient2str(int o) {
    // converts an orientation int to a string
    if(o == 0) {
      return "N";
    } else if(o == 1) {
      return "E";
    } else if(o == 2) {
      return "S";
    } else if(o == 3) {
      return "W";
    }
    return "";
  }
/*
  int orient2degree(int o) {
    // converts an orientation int to degrees
    if(o == 0) {
      return 0;
    } else if(o == 1) {
      return 90;
    } else if(o == 2) {
      return 180;
    } else if(o == 3) {
      return 270;
    }
    return -1;
  }
*/
  int wrap_orientation(int kX, int const kLowerBound, int const kUpperBound) {
    // wraps an orientation int post rotation
    int range_size = kUpperBound - kLowerBound + 1;
    if (kX < kLowerBound)
        kX += range_size * ((kLowerBound - kX) / range_size + 1);

    return kLowerBound + (kX - kLowerBound) % range_size;
  }

  // TODO fix rotation origin?
  void setRotation(int r) {
    // rotates Node
    int o = this->orientation;
    int rot_deg = 90*r;//orient2degree(r);

    // rotate polygon about its origin
    model::polygon<model::d2::point_xy<double> > tmp;
    trans::rotate_transformer<boost::geometry::degree, double, 2, 2> rotate(rot_deg);
    boost::geometry::transform(this->poly, tmp, rotate);
    updateCoordinates();
    this->orientation = wrap_orientation(o + r, 0, 4);
  }

  void updateCoordinates() {
    // Updates parameters of Node class from a geometry object
    auto it = boost::begin(boost::geometry::exterior_ring(this->poly));
    this -> xCoordinate = bg::get<0>(*it);
    this -> yCoordinate = bg::get<1>(*it);

    // centroid
    boost::geometry::model::d2::point_xy<double> centroid;
    boost::geometry::centroid(this->poly, centroid);
    this -> xBy2 = centroid.get<0>();
    this -> yBy2 = centroid.get<1>();
  }

  void printExterior() {
    // print polygon vertices
    for(auto it = boost::begin(boost::geometry::exterior_ring(this->poly)); it != boost::end(boost::geometry::exterior_ring(this->poly)); ++it) {
        double x = bg::get<0>(*it);
        double y = bg::get<1>(*it);
        cout << x << " " << y << endl;
    }
  }

  void printParameter() {
    // print parameters of Node class
    cout << "name      " << this -> name << endl;
    cout << "Width         " << this -> width << endl;
    cout << "Height        " << this -> height << endl;
    cout << "Weight        " << this -> weight << endl;
    cout << "X_Co-ordinate " << this -> xCoordinate << endl;
    cout << "Y_Co-ordinate " << this -> yCoordinate << endl;
    cout << "X/2           " << xBy2 << endl;
    cout << "Y/2           " << yBy2 << endl;
    cout << "Orientation   " << this -> orientation << endl;
    cout << "terminal      " << this -> terminal << endl;
    cout << "NetList       ";
    vector < int > ::iterator it2;
    for (it2 = Netlist.begin(); it2 != Netlist.end(); ++it2) {
      cout << * it2 << " ";
    }
    cout << "\n" << endl;
  }
};
