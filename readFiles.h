#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <regex>
#include <map>

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/adapted/boost_tuple.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/adapted/c_array.hpp>

using namespace std;
using namespace boost::geometry;
namespace trans = boost::geometry::strategy::transform;
namespace bg = boost::geometry;
typedef boost::geometry::model::d2::point_xy<double> point_type;


class Node;
class Pin;
extern map < string, int > name2id;
extern vector < Node > nodeId;

int readNodesFile(string fname);
int readShapesFile(string fname);
int readWtsFile(string fname);
int readPlFile(string fname);
map<int, vector<Pin> > readNetsFile(string fname);
int writePlFile(string fname);

class Pin {
  public:
  string name;
  double x_offset;
  double y_offset;
  int idx;
  double depth;

  void set_params(string nn, double xoffset, double yoffset, int idx) {
    this -> name = nn;
    this -> x_offset = xoffset;
    this -> y_offset = yoffset;
    this -> idx = idx;
  }
};

class Node {
  public:
  string name;
  model::polygon<model::d2::point_xy<double> > poly;
  boost::geometry::model::box< model::d2::point_xy<double> > envelope;
  int idx;
  double width;
  double height;
  int weight;
  bool terminal;
  bool fixed;
  bool overlap;
  double sigma;
  double xCoordinate;
  double yCoordinate;
  double xBy2;
  double yBy2; // change to floats
  string orientation_str;
  int init_orientation;
  int orientation;
  vector < int > Netlist;

  void setParameterNodes(string name, double width, double height, int terminal, int idx) {
    // Sets parameters given an entry in Nodes file
    this -> name = name;
    this -> width = width;
    this -> height = height;
    this -> terminal = terminal;
    this -> orientation_str = "N";
    this -> idx = idx;
    if (!terminal) {
      double points[][2] = {{0.0, 0.0}, {width, 0.0}, {width, height}, {0.0, height}};
      model::polygon< model::d2::point_xy<double> > poly;
      append(poly, points);
      boost::geometry::correct(poly);
      this -> poly = poly;
      boost::geometry::envelope(poly, envelope);
    }
  }

  void setParameterShapes(string wkt) {
    // Sets parameters given an entry in Shapes file
    model::polygon< model::d2::point_xy<double> > poly;
    boost::geometry::read_wkt(wkt, poly);
    boost::geometry::correct(poly);
    this -> poly = poly;
  }

  void setParameterWts(int weight) {
    // Sets parameters given an entry in weights file.
    // weights can correspond to size, importance, etc.
    this -> weight = weight;
  }

  void setParameterPl(double xCoordinate, double yCoordinate, string orientation_str, int fixed) {
    // Sets parameters given an entry in Pl file
    this -> xCoordinate = 0.0;
    this -> yCoordinate = 0.0;
    this -> setPos(xCoordinate, yCoordinate);
    this -> orientation_str = orientation_str;
    this -> init_orientation = str2orient(orientation_str);
    this -> orientation = 0;
    this -> setRotation(this->init_orientation);
    this -> fixed = fixed;
    this -> sigma = 10.0;
  }

  void setNetList(int NetId) {
    // Sets parameters given an entry in Nets file
    Netlist.push_back(NetId);
  }

  void setPos(double x, double y) {
    // Sets position of a node object (lower-left corner)
    if(!this->terminal) {
        trans::translate_transformer<double, 2, 2> translate(x - this->xCoordinate, y - this->yCoordinate);
        model::polygon<model::d2::point_xy<double> > tmp;
        boost::geometry::transform(this->poly, tmp, translate);
        this->poly = tmp;
        updateCoordinates();
    } else {
        this->xCoordinate = x;
        this->yCoordinate = y;
      }
  }

  int str2orient(string o) const{
    if(o == "N") {
      return 0;
    } else if(o == "NE") {
      return 1;
    } else if(o == "E") {
      return 2;
    } else if(o == "SE") {
      return 3;
    } else if (o == "S") {
      return 4;
    } else if (o == "SW") {
      return 5;
    } else if (o == "W") {
      return 6;
    } else if (o == "NW") {
      return 7;
    }
    return -1;
  }

  string orient2str(int o) const{
    if(o == 0) {
      return "N";
    } else if(o == 1) {
      return "NE";
    } else if(o == 2) {
      return "E";
    } else if(o == 3) {
      return "SE";
    } else if(o == 4) {
      return "S";
    } else if(o == 5) {
      return "SW";
    } else if(o == 6) {
      return "W";
    } else if(o == 7) {
      return "NW";
    }
    return "";
  }

  int wrap_orientation(int kX) {
    return kX % 8;
  }

  void setRotation(int r) {
    int rot_deg = 45*r;
    double tmpx = this->xCoordinate;
    double tmpy = this->yCoordinate;

    // rotate polygon about its origin
    model::polygon<model::d2::point_xy<double> > tmp;
    trans::rotate_transformer<boost::geometry::degree, double, 2, 2> rotate(rot_deg);
    boost::geometry::transform(this->poly, tmp, rotate);
    this->poly = tmp;
    updateCoordinates();
    this->orientation = wrap_orientation(this->orientation + r);
    this->setPos(tmpx, tmpy);
  }

  /*
  upateCoordinates
  Updates parameters of Node class from a geometry object
  */
  void updateCoordinates() {
    if(this->terminal) {
        this -> xBy2 = this->xCoordinate;
        this -> yBy2 = this->yCoordinate;
    } else {
        boost::geometry::model::d2::point_xy<double> centroid;
        boost::geometry::centroid(this->poly, centroid);
        this -> xBy2 = centroid.get<0>();
        this -> yBy2 = centroid.get<1>();
        boost::geometry::envelope(poly, envelope);
        this->xCoordinate = bg::get<bg::min_corner, 0>(this->envelope);
        this->yCoordinate = bg::get<bg::min_corner, 1>(this->envelope);
    }
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
