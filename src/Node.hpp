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
    vector < int > Netlist;

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

};
