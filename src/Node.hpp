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
using namespace boost::geometry;
namespace trans = boost::geometry::strategy::transform;
namespace bg = boost::geometry;
typedef boost::geometry::model::d2::point_xy<double> Point;

/**
Node
Base Node class for board components
*/
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
    double initialX;
    double initialY;
    double xBy2;
    double yBy2;
    string orientation_str;
    int init_orientation;
    int orientation;
    int layer=1; // 1 or -1 for 2-sided placement
    int mirror;
    int flipped=-1;
    vector < int > Netlist;
/*
    void setParameterNodes(string _name, double _width, double _height, bool _terminal, int _idx, int _mirror=0);
    void setParameterShapes(string wkt);
    void setParameterWts(int _weight);
    void setParameterPl(double xCoordinate, double yCoordinate, string _orientation_str, bool _fixed);
    void setNetList(int NetId);
    void setPos(double x, double y);
    void layerChange();
    int wrap_orientation(int kX);
    void setRotation(int r);
    void updateCoordinates();
    int str2orient(string o) const;
    string orient2str(int o) const;
    void printExterior() const;
    void printParameter();
*/
/**
setParameterNodes
Sets parameters given an entry in Nodes file
*/
void Node::setParameterNodes(string _name, double _width, double _height, bool _terminal, int _idx, int _mirror=0) {
  name = _name;
  xCoordinate = 0.0;
  yCoordinate = 0.0;
  orientation = 0;
  width = _width;
  height = _height;
  terminal = _terminal;
  mirror = _mirror;
  if (_mirror == -1) {
    layer = -1;
  }
  orientation_str = "N";
  idx = _idx;
  if (!terminal) {
    double points[][2] = {{0.0, 0.0}, {_width, 0.0}, {_width, _height}, {0.0, _height}};
    model::polygon< model::d2::point_xy<double> > _poly;
    append(_poly, points);
    boost::geometry::correct(_poly);
    poly = _poly;
    boost::geometry::envelope(_poly, envelope);
  }
}

/**
setParameterShapes
Sets parameters given an entry in Shapes file
*/
void Node::setParameterShapes(string wkt) {
  model::polygon< model::d2::point_xy<double> > _poly;
  boost::geometry::read_wkt(wkt, _poly);
  boost::geometry::correct(_poly);
  poly = _poly;
}

/**
setParameterWts
Sets parameters given an entry in weights file. 
*/
void Node::setParameterWts(int _weight) {
  weight = _weight;
}

/**
setParameterPl
Sets parameters given an entry in Pl file  
*/
void Node::setParameterPl(double xCoordinate, double yCoordinate, string _orientation_str, bool _fixed) {
  setPos(xCoordinate, yCoordinate);
  this -> initialX = xCoordinate;
  this -> initialY = yCoordinate;
  orientation_str = _orientation_str;
  init_orientation = str2orient(_orientation_str);
  setRotation(str2orient(_orientation_str));
  fixed = _fixed;
  sigma = 20.0;
}

/**
setNetList
Sets parameters given an entry in Nets file  
*/
void Node::setNetList(int NetId) {
  Netlist.push_back(NetId);
}

/**
setPos
Sets position of a node object (lower-left corner)
*/
void Node::setPos(double x, double y) {
  if(!terminal) {
      trans::translate_transformer<double, 2, 2> translate(x - xCoordinate, y - yCoordinate);
      model::polygon<model::d2::point_xy<double> > tmp;
      boost::geometry::transform(poly, tmp, translate);
      poly = tmp;
      updateCoordinates();
  } else {
      xCoordinate = x;
      yCoordinate = y;
    }
}

void Node::layerChange() {
  if(!terminal) {
    //model::polygon<model::d2::point_xy<double> > tmp;
    //boost::geometry::transform(poly, tmp, translate);
    //poly = tmp;
    layer=-1*layer;
    //updateCoordinates();  
  }  else {
    layer=-1*layer;
  }
  flipped = -1 * flipped;
}

int Node::wrap_orientation(int kX) {
  return kX % 8;
}

/**
setRotation
Rotate polygon about global origin and transform back to local origin
*/
void Node::setRotation(int r) {
  int rot_deg = 45*r;
  double tmpx = xCoordinate;
  double tmpy = yCoordinate;

  setPos(-width/2,-height/2);
  model::polygon<model::d2::point_xy<double> > tmp;
  trans::rotate_transformer<boost::geometry::degree, double, 2, 2> rotate(rot_deg);
  boost::geometry::transform(poly, tmp, rotate);
  poly = tmp;
  double otmp = wrap_orientation(orientation + r);
  orientation = otmp;
  setPos(tmpx, tmpy);
}

/**
upateCoordinates
Updates parameters of Node class from a geometry object
*/
void Node::updateCoordinates() {
  if(terminal) {
      xBy2 = xCoordinate;
      yBy2 = yCoordinate;
  } else {
      boost::geometry::model::d2::point_xy<double> centroid;
      boost::geometry::centroid(poly, centroid);
      xBy2 = centroid.get<0>();
      yBy2 = centroid.get<1>();

      boost::geometry::model::box< model::d2::point_xy<double> > envtmp;
      boost::geometry::envelope(poly, envtmp);
      envelope = envtmp;
      Point minCorner = envelope.min_corner();
      Point maxCorner = envelope.max_corner();
      width  = maxCorner.get<0>() - minCorner.get<0>();
      height = maxCorner.get<1>() - minCorner.get<1>();
      xCoordinate = bg::get<bg::min_corner, 0>(envelope);
      yCoordinate = bg::get<bg::min_corner, 1>(envelope);
  }
}

int Node::str2orient(string o) const{
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

string Node::orient2str(int o) const{
  if(o == 0 || o == 8) {
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

/**
printExterior
print polygon vertices
*/
void Node::printExterior() const{
  for(auto it = boost::begin(boost::geometry::exterior_ring(poly)); it != boost::end(boost::geometry::exterior_ring(poly)); ++it) {
      double x = bg::get<0>(*it);
      double y = bg::get<1>(*it);
      cout << x << " " << y << endl;
  }
}

/**
printParameter
print node params
*/
void Node::printParameter() {
  cout << "name      " << name << endl;
  cout << "Width         " << width << endl;
  cout << "Height        " << height << endl;
  cout << "Weight        " << weight << endl;
  cout << "X_Co-ordinate " << xCoordinate << endl;
  cout << "Y_Co-ordinate " << yCoordinate << endl;
  cout << "X/2           " << xBy2 << endl;
  cout << "Y/2           " << yBy2 << endl;
  cout << "Orientation   " << orientation << endl;
  cout << "terminal      " << terminal << endl;
  cout << "fixed         " << fixed << endl;
  cout << "layer         " << layer << endl;
  cout << "NetList       ";
  vector < int > ::iterator it2;
  for (it2 = Netlist.begin(); it2 != Netlist.end(); ++it2) {
    cout << *it2 << " ";
  }
  cout << "\n" << endl;
}
};
