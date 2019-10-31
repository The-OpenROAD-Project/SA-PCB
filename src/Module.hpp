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
#include "node.hpp"

using namespace std;
namespace trans = boost::geometry::strategy::transform;
namespace bg = boost::geometry;
typedef boost::geometry::model::d2::point_xy<double> Point;

class Module {
  public:
    string name;
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
    bool root_module;
    string orientation_str;
    int orientation;
    vector < int > Netlist;


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

  int wrapOrientation(int kX) {
    return kX % 8;
  }

  /*
  setRotation
  Rotate polygon about global origin and transform back to local origin
  */
  void setRotation(int r) {
    int rot_deg = 45*r;
    double tmpx = this->xCoordinate;
    double tmpy = this->yCoordinate;

    this->setPos(-this->width/2,-this->height/2);
    model::polygon<model::d2::point_xy<double> > tmp;
    trans::rotate_transformer<boost::geometry::degree, double, 2, 2> rotate(rot_deg);
    boost::geometry::transform(this->poly, tmp, rotate);
    this->poly = tmp;
    double otmp = wrapOrientation(this->orientation + r);
    this->orientation = otmp;
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
    this -> printExterior();
    cout << "\n" << endl;
  }
};
