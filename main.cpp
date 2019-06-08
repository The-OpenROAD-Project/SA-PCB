

#include "main.h"

using namespace std::chrono;

BOOST_GEOMETRY_REGISTER_BOOST_TUPLE_CS(cs::cartesian)

namespace bg = boost::geometry;
namespace bnu = boost::numeric::ublas;
typedef boost::geometry::model::polygon<boost::geometry::model::d2::point_xy<double> > polygon;

double Temperature;
double sigma = 25.0;
int xLimit;
map < int, row > rowId;
int RowWidth;
map<int, vector<Pin> > netToCell;
map < string, Node > nodeId;
boundaries b;

float l1 = 0.8;
float l2 = 1-l1;
std::pair <double,double> wl_normalization;
std::pair <double,double> area_normalization;
std::pair <double,double> routability_normalization;

long long int iii = 0;

int debug = 1;
long long int idx = -1;

vector< int > accept_history;
float accept_ratio = 0;
vector< double > accept_ratio_history;

boost::mt19937 rng;

int main(int argc, char *argv[]) {
  srand(time(NULL));
  int opt;
  int outer_loop_iter = 2000;
  int inner_loop_iter = 20;
  double eps = -1.0;
  double t_0 = 1.0;

  string parg = "";
  int bound_inc = 0;
  while ((opt = getopt(argc,argv,"i:x:j:t:f:p:b:")) != EOF) {
      switch(opt) {
          case 'x': idx = atoi(optarg); break;
          case 'p': parg=string(optarg); break;
          case 'd': debug=atoi(optarg); break;
          case 'i': outer_loop_iter = atoi(optarg); break;
          case 'j': inner_loop_iter = atoi(optarg); break;
          case 't': t_0 = atoi(optarg); break;
          case 'e': eps = atof(optarg); break;
          case 'b':
            bound_inc += 1;
            switch(bound_inc) {
              case 1: b.minX = atof(optarg);
              case 2: b.minY = atof(optarg);
              case 3: b.maxX = atof(optarg);
              case 4: b.maxY = atof(optarg);
            }
          case 'h': fprintf(stderr, "./sa usuage is \n \
                                    -i <value> : for denoting # outer iterations PER SA INSTANCE \n \
                                    -j <value> : for denoting 'j'*#nodes inner iterations \n \
                                    -t <value> : for denoting initial temperature \n \
                                    -f <str>   : for denoting number output pl \n \
                                    -e <float> : convergence epsilon \n \
                                    EXAMPLE: ./sa -i 20000 -j 20 -t 1 -f output.pl for the standard single instance timberwolf algorithm");
          default: cout<<endl; abort();
      }
  }
  if(string(parg) == "")   {
    cout << "./main: option requires an argument -- p\n";
    exit(1);
  }
  string p = string(parg);
  cout << "circuit: " << p << endl;

  string nodesfname = p + ".nodes";
  string netsfname  = p + ".nets";
  string plfname    = p + ".pl";
  string wtsfname   = p + ".wts";

  cout << nodesfname << endl;
  cout << netsfname << endl;
  cout << plfname << endl;

  if(debug) { cout << "reading nodes..." << endl; }
  readNodesFile(nodesfname);
  //readWtsFile(wtsfname);
  if(debug) { cout << "reading pl..." << endl; }
  readPlFile(plfname);
  if(debug) { cout << "reading nets..." << endl; }
  readNetsFile(netsfname);
  //readSclFile();
  if(debug) { cout << "calculating boundaries..." << endl; }
  CalcBoundaries();
  if(debug) { cout << "annealing" << endl; }
  float cost = timberWolfAlgorithm(outer_loop_iter, inner_loop_iter, eps, t_0);
  writePlFile("./"+std::to_string( idx )+".pl");
  return cost;
}

/*
SetInitParameters
Empirically finds normalization parameters for scaling cost terms.
Currently, we scale by 1/(f) where f is the cost of the inital placement.
*/
void SetInitParameters() {
  vector < std::pair <double,double> > normalization_terms;

  int num_components = 0;
  map < string, Node > ::iterator itNode;
  for (itNode = nodeId.begin(); itNode != nodeId.end(); ++itNode) {
    if(!itNode -> second.terminal) {
      num_components += 1;
    }
  }

  // min/max wl
  std::pair <double,double> wl;
  double min_wl = std::numeric_limits<double>::max();
  double max_wl = 0.0;

  // min/max overlap area
  std::pair <double,double> area;
  double min_area = 0.0;
  double max_area = 0.0;

  map < string, Node > ::iterator nodeit1 = nodeId.begin();
  map < string, Node > ::iterator nodeit2 = nodeId.begin();
  vector < string > hist;

  for (nodeit1 = nodeId.begin(); nodeit1 != nodeId.end(); ++nodeit1) {
    if(nodeit1->second.terminal) {
      continue;
    } else {
        hist.push_back(nodeit1->second.name);
        double cell_wl = nodeit1->second.width + nodeit1->second.height;
        if (cell_wl < min_wl) {
          min_wl = cell_wl;
        }

        for (nodeit2 = nodeId.begin(); nodeit2 != nodeId.end(); ++nodeit2) {
          if (!nodeit2->second.terminal && !(find(hist.begin(), hist.end(), nodeit2->second.name) != hist.end()) ) {
            if(nodeit1->second.fixed && nodeit2->second.fixed) {
                if(intersects(nodeit1->second.poly, nodeit2->second.poly)) {
                  double oa = 0.0;
                  std::deque<polygon> intersect_poly;
                  boost::geometry::intersection(nodeit1->second.poly, nodeit2->second.poly, intersect_poly);

                  BOOST_FOREACH(polygon const& p, intersect_poly) {
                      oa += bg::area(p);
                  }
                  min_area += pow(oa,2);
                  max_area += pow(oa,2);
                }
              } else {
                max_area += pow(min(nodeit1->second.width, nodeit2->second.width) * min(nodeit1->second.height, (nodeit2->second.width)), 2);
              }
          }
        }
      }
  }

  min_wl = min_wl*num_components;
  max_wl = ((b.maxX - b.minX) + (b.maxY - b.minY))*netToCell.size();

  wl.first = 0.0;
  wl.second = wireLength();

  area.first = 0.0;
  area.second = cellOverlap();

/*
  wl.first = min_wl;
  wl.second = max_wl;

  area.first = min_area;
  area.second = max_area;
*/
  normalization_terms.push_back(wl);
  normalization_terms.push_back(area);

  wl_normalization = normalization_terms[0];
  area_normalization = normalization_terms[1];

  routability_normalization.first = 0.0;
  routability_normalization.second = rudy();

}

/*
CalcBoundaries
Calculate board boundaries by finding the maximum-area
rectangular envelop encompasing terminals & fixed modules.
*/
void CalcBoundaries() {
  double xval, yval, width, height;
  map < string, Node > ::iterator itNode;
  for (itNode = nodeId.begin(); itNode != nodeId.end(); ++itNode) {
    if (itNode -> second.terminal || itNode->second.fixed) {
      xval = itNode -> second.xCoordinate;
      yval = itNode -> second.yCoordinate;
      width = itNode -> second.width;
      height =itNode -> second.height;
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

  if(debug) {
    cout << "min: " << b.minX << " " << b.minY << endl;
    cout << "max: " << b.maxX << " " << b.maxY << endl;
  }
}

/*
macroPlacement [DEPRECIATED]
Preplace macros/large components
*/
int macroPlacement() {
  // macro placement [depreceiated]
  int xValue = b.minX, yValue = 0;
  xLimit = b.minX;
  map < int, row > ::iterator itRow;
  itRow = rowId.begin();
  int rowHeight = itRow -> second.height;
  for (map < string, Node > ::iterator itNode = nodeId.begin(); itNode != nodeId.end(); ++itNode) {
    if (!itNode -> second.terminal && itNode -> second.height > rowHeight) {
      if (xValue + itNode -> second.width > xLimit) {
        xLimit = xValue + itNode -> second.width + 1;
      }
      if ((yValue + itNode -> second.height) < b.maxY) {
        itNode -> second.yCoordinate = yValue;
        itNode -> second.xCoordinate = xValue;
      } else {

        yValue = 0;
        xValue = xLimit;
        itNode -> second.yCoordinate = yValue;
        itNode -> second.xCoordinate = xValue;
      }
      yValue = yValue + itNode -> second.height;
    }
  }
  return xLimit;
}

/*
randomPlacement
Randomly place and orient a single component within a bounded region
*/
void randomPlacement(int xmin, int xmax, int ymin, int ymax, Node n) {
  boost::uniform_int<> uni_distx(xmin,xmax);
  boost::variate_generator<boost::mt19937&, boost::uniform_int<> > unix(rng, uni_distx);
  int rx = unix();

  boost::uniform_int<> uni_disty(ymin,ymax);
  boost::variate_generator<boost::mt19937&, boost::uniform_int<> > uniy(rng, uni_disty);
  int ry = uniy();

  boost::uniform_int<> uni_disto(0,4);
  boost::variate_generator<boost::mt19937&, boost::uniform_int<> > unio(rng, uni_disto);
  int ro = unio();

  string ostr = n.orient2str(ro);
  n.setParameterPl(rx,ry,ostr, n.fixed);
}

/*
initialPlacement
Randomly place and orient all movable components in the board area
*/
void initialPlacement() {
  map < string, Node > ::iterator itNode;
  for (itNode = nodeId.begin(); itNode != nodeId.end(); ++itNode) {
    Node n = itNode -> second;
    if (!n.fixed) {
      randomPlacement(b.minX, b.maxX - n.width, b.minY, b.maxY - n.height, n);
    }
  }
}

/*
project_soln
projects the final solution to enforce 0 overlap by iteratively
shifting overlapped components until cumulative overlap reduces below an eps.
*/
void project_soln() {
  // ensure zero overlap
  double oa = cellOverlap();
  double eps = 10.0e-5;
  map < string, Node > ::iterator nodeit = nodeId.begin();
  map < string, Node > ::iterator nodeit2 = nodeId.begin();

  boost::geometry::model::box< model::d2::point_xy<double> > env1;
  boost::geometry::model::box< model::d2::point_xy<double> > env2;

  while (oa > eps) {
    for (nodeit = nodeId.begin(); nodeit != nodeId.end(); ++nodeit) {
      if (!nodeit->second.terminal) {
        env1 = nodeit->second.envelope;
        for (nodeit2 = nodeId.begin(); nodeit2 != nodeId.end(); ++nodeit2) {
          if(!intersects(nodeit->second.poly, nodeit2->second.poly) || (nodeit->second.fixed && nodeit2->second.fixed)) {
            continue;
          } else {
            // Project out of envelope
            env2 = nodeit2->second.envelope;

            double dx = 0.0;
            double dy = 0.0;
            double rx = nodeit->second.xCoordinate;
            double ry = nodeit->second.yCoordinate;
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
              validateMove(&nodeit->second,rx + dx,ry);
            } else {
              // project in y direction
              validateMove(&nodeit->second,rx,ry + dy);
            }
          }
        }
      }
    }
  }
}

/*
wireLength
Computes HPWL for all nets
*/
double wireLength() {
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
      int orient = nodeId[itCellList->name].orientation;
      xVal = nodeId[itCellList->name].xBy2;
      yVal = nodeId[itCellList->name].yBy2;
      if(orient == 0) {
        xVal = xVal + itCellList->x_offset;
        yVal = yVal + itCellList->y_offset;
      } else if(orient == 1) {
        xVal = xVal + itCellList->y_offset;
        yVal = yVal - itCellList->x_offset;
      } else if(orient == 2) {
        xVal = xVal - itCellList->x_offset;
        yVal = yVal - itCellList->y_offset;
      } else if(orient == 3) {
        xVal = xVal - itCellList->y_offset;
        yVal = yVal + itCellList->x_offset;
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
cellOverlap
Compute sum squared overlap for all components
*/
double cellOverlap() {
  double overlap = 0.0;
  map < string, Node > ::iterator nodeit = nodeId.begin();
  map < string, Node > ::iterator nodeit2 = nodeId.begin();

  vector < string > hist;
  for (nodeit = nodeId.begin(); nodeit != nodeId.end(); ++nodeit) {
    hist.push_back(nodeit->second.name);
    if (!nodeit->second.terminal) {
      for (nodeit2 = nodeId.begin(); nodeit2 != nodeId.end(); ++nodeit2) {
        if (!nodeit2->second.terminal && !(find(hist.begin(), hist.end(), nodeit2->second.name) != hist.end()) ) {
          if(!intersects(nodeit->second.poly, nodeit2->second.poly) || (nodeit->second.fixed && nodeit2->second.fixed)) {
            continue;
          } else {
            double oa = 0.0;
            std::deque<polygon> intersect_poly;
            boost::geometry::intersection(nodeit->second.poly, nodeit2->second.poly, intersect_poly);

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
  }
  return overlap;
}

/*
wireLength_partial
Compute HPWL for select nets
*/
double wireLength_partial(map < string, Node > nodes) {
  map<int, vector < Pin > > ::iterator itNet;
  vector < Pin > ::iterator itCellList;
  double xVal, yVal, wireLength = 0;
  for (itNet = netToCell.begin(); itNet != netToCell.end(); ++itNet) {
    if (nodes.find("f") == nodes.end()) {
      continue;
    }
    double minXW = b.maxX, minYW = b.maxY, maxXW = b.minX, maxYW = b.minY;
    for (itCellList = itNet -> second.begin(); itCellList != itNet -> second.end(); ++itCellList) {
      if(itCellList->name == "") {
        continue;
      }
      int orient = nodeId[itCellList->name].orientation;
      xVal = nodeId[itCellList->name].xBy2;
      yVal = nodeId[itCellList->name].yBy2;
      if(orient == 0) {
        xVal = xVal + itCellList->x_offset;
        yVal = yVal + itCellList->y_offset;
      } else if(orient == 1) {
        xVal = xVal + itCellList->y_offset;
        yVal = yVal - itCellList->x_offset;
      } else if(orient == 2) {
        xVal = xVal - itCellList->x_offset;
        yVal = yVal - itCellList->y_offset;
      } else if(orient == 3) {
        xVal = xVal - itCellList->y_offset;
        yVal = yVal + itCellList->x_offset;
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
cellOverlap_partial
Compute sum squared overlap for select components
*/
double cellOverlap_partial(map < string, Node > nodes) {
  double overlap = 0.0;
  map < string, Node > ::iterator nodeit = nodes.begin();
  map < string, Node > ::iterator nodeit2 = nodeId.begin();

  vector < string > hist;
  for (nodeit = nodeId.begin(); nodeit != nodeId.end(); ++nodeit) {
    hist.push_back(nodeit->second.name);
    if (!nodeit->second.terminal) {
      for (nodeit2 = nodeId.begin(); nodeit2 != nodeId.end(); ++nodeit2) {
        if (!nodeit2->second.terminal && !(find(hist.begin(), hist.end(), nodeit2->second.name) != hist.end()) ) {
          if(!intersects(nodeit->second.poly, nodeit2->second.poly) || (nodeit->second.fixed && nodeit2->second.fixed)) {
            continue;
          } else {
            double oa = 0.0;
            std::deque<polygon> intersect_poly;
            boost::geometry::intersection(nodeit->second.poly, nodeit2->second.poly, intersect_poly);

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
  }
  return overlap;
}

/*
rudy
Computes a routability score
*/
double rudy() {
  bnu::matrix<double> D (static_cast<int>(abs(b.maxY)+abs(b.minY)+1), static_cast<int>(abs(b.maxX)+abs(b.minX)+1), 0);
  bnu::matrix<double> D_route_sup (static_cast<int>(abs(b.maxY)+abs(b.minY)+1), static_cast<int>(abs(b.maxX)+abs(b.minX)+1), 1);

  map<int, vector < Pin > > ::iterator itNet;
  vector < Pin > ::iterator itCellList;
  for (itNet = netToCell.begin(); itNet != netToCell.end(); ++itNet) {
    double xVal, yVal, hpwl = 0.0;
    double minXW = b.maxX, minYW = b.maxY, maxXW = b.minX, maxYW = b.minY;
    double rudy = 0.0;
    bnu::matrix<double> D_net (static_cast<int>(abs(b.maxY)+abs(b.minY)+1), static_cast<int>(abs(b.maxX)+abs(b.minX)+1), 0);
    for (itCellList = itNet -> second.begin(); itCellList != itNet -> second.end(); ++itCellList) {
      if(itCellList->name == "") {
        continue;
      }
      int orient = nodeId[itCellList->name].orientation;
      xVal = nodeId[itCellList->name].xBy2;
      yVal = nodeId[itCellList->name].yBy2;
      if(orient == 0) {
        xVal = xVal + itCellList->x_offset;
        yVal = yVal + itCellList->y_offset;
      } else if(orient == 1) {
        xVal = xVal + itCellList->y_offset;
        yVal = yVal - itCellList->x_offset;
      } else if(orient == 2) {
        xVal = xVal - itCellList->x_offset;
        yVal = yVal - itCellList->y_offset;
      } else if(orient == 3) {
        xVal = xVal - itCellList->y_offset;
        yVal = yVal + itCellList->x_offset;
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
    // now have boundary of net
    hpwl = (abs((maxXW - minXW)) + abs((maxYW - minYW)));
    rudy = hpwl / (max((maxXW - minXW)*(maxYW - minYW), 1.0));

    // set read_net
    for (unsigned i = minYW; i < maxYW; ++ i) {
        for (unsigned j = minXW; j < maxXW; ++ j) {
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
  if(iii % 10 == 0) {
      ofstream dat("cache_rudy/"+std::to_string(iii) + ".txt");
      for (unsigned i = 0; i < D.size1() ; i++) {
          for (unsigned j = 0; j < D.size2(); j++) {
              dat << D(i, j) << "\t";
          }
          dat << endl;
      }
  }

  r = r/(wl_normalization.second * D.size1() * D.size2()); // normalize r
  return r;
}

double cost(int temp_debug) {
  if(debug > 1 || temp_debug == -1) {
    cout << "wirelength: " << wireLength() << endl;
    cout << "overlap: " << cellOverlap() << endl;
    cout << "l1: " << l1 << " l2: " << l2 << endl;
    cout << "wirelength_cost: " << l1*(wireLength() - wl_normalization.first)/(wl_normalization.second - wl_normalization.first) << endl;
    cout << "overlap_cost: " << l2  * 0.99 * (cellOverlap() - area_normalization.first)/(area_normalization.second - area_normalization.first) << endl;
    cout << "routability_cost: " << l2  * 0.01  * rudy() << endl;
    cout << "cost: " << l1*(wireLength() - wl_normalization.first)/(wl_normalization.second - wl_normalization.first) +
                        l2  * 0.99 * (cellOverlap() - area_normalization.first)/(area_normalization.second - area_normalization.first)  +
                        l2 * 0.01*rudy() << endl;
  }
  return l1 * (wireLength() - wl_normalization.first)/(wl_normalization.second - wl_normalization.first) +
         l2 * 0.9 * (cellOverlap() - area_normalization.first)/(area_normalization.second - area_normalization.first) +
         l2 * 0.1 * rudy();
         //l2 * 0.1 * (rudy() - routability_normalization.first)/(routability_normalization.second - routability_normalization.first);
}

double cost_partial(int temp_debug, map < string, Node > nodes) {
  return l1 * (wireLength_partial(nodes) - wl_normalization.first)/(wl_normalization.second - wl_normalization.first) +
         l2 * 0.9 * (cellOverlap_partial(nodes) - area_normalization.first)/(area_normalization.second - area_normalization.first) +
		     l2 * 0.1 * rudy();
}

/*
random_node
Select random node from set of nodes
*/
map < string, Node > ::iterator random_node() {
  map < string, Node > ::iterator itNode = nodeId.begin();
  int size = nodeId.size();
  boost::uniform_int<> uni_dist(0,size-1);
  boost::variate_generator<boost::mt19937&, boost::uniform_int<> > uni(rng, uni_dist);
  int randint = uni();
  std::advance(itNode, randint);
  return itNode;
}

/*
validateMove
validates a shift, project within board boundary
*/
void validateMove(Node* node, double rx, double ry) {
  int orient = node->orientation;
  double width = 0.0;
  double height = 0.0;

  if (orient == 0 || orient == 2) {
    width = node->width;
    height = node->height;
  } else {
    width = node->height;
    height = node->width;
  }
  rx = max(rx, b.minX);
  ry = max(ry, b.minY);

  rx = min(rx, b.maxX - width);
  ry = min(ry, b.maxY - height);

  node->setPos(rx,ry);
}

double c = 0.0;
double initiateMove() {
  // Initate a transition
  int state = -1;
  double prevCost = cost();
  if(debug > 1) {
    cout << "prevCost: " << prevCost << endl;
  }
  map < string, Node > ::iterator rand_node1;
  map < string, Node > ::iterator rand_node2;
  int r = 0;
  rand_node1 = random_node();
  while(rand_node1->second.terminal || rand_node1->second.name == "" || rand_node1->second.fixed) {
    rand_node1 = random_node();
  }
  double rand_node1_orig_x = rand_node1->second.xCoordinate;
  double rand_node1_orig_y = rand_node1->second.yCoordinate;
  double rand_node2_orig_x = 0.0;
  double rand_node2_orig_y = 0.0;

  if(debug > 1) {
    cout << "=======" << endl;
    cout <<  "name: " << rand_node1->second.name << endl;
  }

  boost::uniform_real<> uni_dist(0,1);
  boost::variate_generator<boost::mt19937&, boost::uniform_real<> > uni(rng, uni_dist);
  double i = uni();

  if (i<0.15) { // swap
    state = 0;
    if(debug > 1) {
      cout << "swap" << endl;
    }
    rand_node2 = random_node();
    while(rand_node2->second.terminal || rand_node2->first == rand_node1->first || rand_node2->second.name == "" || rand_node2->second.fixed) {
      rand_node2 = random_node();
    }
    rand_node2_orig_x = rand_node2->second.xCoordinate;
    rand_node2_orig_y = rand_node2->second.yCoordinate;

    validateMove(&rand_node1->second, rand_node2_orig_x, rand_node2_orig_y);
    validateMove(&rand_node2->second, rand_node1_orig_x, rand_node1_orig_y);
  } else if (i < 0.75) { // shift
    state = 1;
    if(debug > 1) {
      cout << "shift" << endl;
    }
    //double sigma = rand_node1->second.sigma;

    boost::normal_distribution<> nd(0.0, sigma);
    boost::variate_generator<boost::mt19937&,
                             boost::normal_distribution<> > var_nor(rng, nd);

    double dx = var_nor();
    double dy = var_nor();

    double rx = rand_node1_orig_x + dx;
    double ry = rand_node1_orig_y + dy;

    validateMove(&rand_node1->second, rx, ry);
  } else if (i < 0.2) { // rotate
    state = 2;
    if(debug > 1) {
      cout << "rotate" << endl;
    }
    boost::uniform_int<> uni_dist(0,3);
    boost::variate_generator<boost::mt19937&, boost::uniform_int<> > uni(rng, uni_dist);
    int r = uni();
    rand_node1->second.setRotation(r);
  }
  bool accept = checkMove(prevCost);
  accept_ratio_history.push_back(accept_ratio);
  if (!accept) {
    if(debug > 1) {
      cout << "reject" << endl;
    }
    // revert state
    if (state == 0) {
      rand_node1->second.setPos(rand_node1_orig_x,rand_node1_orig_y);
      rand_node2->second.setPos(rand_node2_orig_x,rand_node2_orig_y);
    } else if (state == 1) {
      rand_node1->second.setPos(rand_node1_orig_x,rand_node1_orig_y);
    } else if (state == 2) {
      rand_node1->second.setRotation(-r);
    }
    return prevCost;
  }
  else {
    if(debug > 1) {
      cout << "accept" << endl;
    }
    return c;
  }
}

/*
update_Temperature
Update the SA parameters according to annealing schedule
*/
void update_temperature() {
  map < string, Node > ::iterator nodeit = nodeId.begin();
  if (Temperature > 0.05) {
    Temperature = (0.99986) * Temperature;
    if (l1 > 10e-2) {
      l1 -= 10e-5;
      l2 += 10e-5;
    }
    sigma = max(0.9999*sigma,0.5);
    for (nodeit = nodeId.begin(); nodeit != nodeId.end(); ++nodeit) {
      nodeit->second.sigma =  max(0.9999*nodeit->second.sigma,0.5);
    }
  } else if (Temperature <= 0.05 && Temperature > 0.005) {
    Temperature = (0.999) * Temperature;
    if (l1 > 25e-4) {
      l1 -= 25e-4;
      l2 += 25e-4;
    }
    sigma = max(0.999*sigma,0.25);
    for (nodeit = nodeId.begin(); nodeit != nodeId.end(); ++nodeit) {
      nodeit->second.sigma =  max(0.999*nodeit->second.sigma,0.25);
    }
  } else if (Temperature < 0.005) {
    Temperature = (0.9985) * Temperature;
    if (l1 > 20e-4) {
      l1 -= 20e-4;
      l2 += 20e-4;
    }
    sigma = max(0.999*sigma,0.18);
    for (nodeit = nodeId.begin(); nodeit != nodeId.end(); ++nodeit) {
      nodeit->second.sigma =  max(0.999*nodeit->second.sigma,0.18);
    }
  } else if (Temperature < 0.001) {
    Temperature = (0.998) * Temperature;
    if (l1 > 10e-5) {
      l1 -= 10e-5;
      l2 += 10e-5;
    }
    sigma = max(0.999*sigma,0.08);
    for (nodeit = nodeId.begin(); nodeit != nodeId.end(); ++nodeit) {
      nodeit->second.sigma =  max(0.999*nodeit->second.sigma,0.08);
    }
  }
}

/*
checkMove
either accept or reject the move based on current & previous temperature & cost
*/
bool checkMove(double prevCost) {
  double newCost = cost();
  c = newCost;
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
    if (accept_history.size() > 100) {
      accept_ratio = ((accept_ratio*100) - (accept_history[accept_history.size()-100]) + 1)/100;
    } else {
      accept_ratio = accumulate(accept_history.begin(), accept_history.end(),0.0) / accept_history.size();
    }
    return true;
  } else {
    accept_history.push_back(0);
    if (accept_history.size() > 100) {
      accept_ratio = ((accept_ratio*100) - (accept_history[accept_history.size()-100]) + 0)/100;
    } else {
      accept_ratio = accumulate(accept_history.begin(), accept_history.end(),0.0) / accept_history.size();
    }
    return false;
  }
}

//void varanelli_cohoon() {
  //return
//}

/*
gen_report
generates a report and outputs files to ./reports/ directory
*/
void gen_report(map<string, vector<double>* > report) {
    vector < double > cost_hist = *report["cost_hist"];
    vector < double > wl_hist   = *report["wl_hist"];
    vector < double > oa_hist   = *report["oa_hist"];

    time_t t = time(0);
    struct tm * now = localtime( & t );

    char buf [80];
    strftime (buf,80,"%Y-%m-%d-%H-%M-%S",now);
    string buffer = string(buf);

    double wl = wireLength();
    double oa = cellOverlap();
    double routability = rudy();
    double normalized_wl = (wl - wl_normalization.first)/(wl_normalization.second - wl_normalization.first);
    double normalized_oa = (oa - area_normalization.first)/(area_normalization.second - area_normalization.first);
    double normalized_routability = routability;
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

    writePlFile("./reports/"+buffer+"_pl.pl");
}

/*
timberWolfAlgorithm
main loop for sa algorithm
*/
float timberWolfAlgorithm(int outer_loop_iter, int inner_loop_iter, double eps, double t_0) {
  //limit = macroPlacement();

  //varanelli_cohoon();

  initialPlacement();
  SetInitParameters();
  Temperature = t_0;

  int num_components = 0;
  map < string, Node > ::iterator itNode;
  map < string, vector < double >* > report;
  vector < double > cost_hist;
  vector < double > wl_hist;
  vector < double > oa_hist;
  report["cost_hist"] = &cost_hist;
  report["wl_hist"] = &wl_hist;
  report["oa_hist"] = &oa_hist;
  double cst = cost();

  for (itNode = nodeId.begin(); itNode != nodeId.end(); ++itNode) {
    if(!itNode -> second.terminal && !itNode -> second.fixed) {
      num_components += 1;
    }
  }

  high_resolution_clock::time_point t1 = high_resolution_clock::now();

  int i; // inner loop iterator
  long long int ii = 0; // outer loop iterator
  while (ii < outer_loop_iter) {
    i = 2*num_components;
    if(ii % 100 == 0 && debug) {
      high_resolution_clock::time_point t2 = high_resolution_clock::now();
      duration<double> time_span = duration_cast< duration<double> >(t2 - t1);

      cout << "******" << ii << "******" << endl;
      cout << "iteration: " << ii << endl;
      cout << "time: " <<  time_span.count() << " (s)" << endl;
      cout << "move/time: " <<  i*ii/time_span.count() << endl;
      cout << "time remaining: " <<  time_span.count()/ii * (1000-ii) << " (s)" << endl;
      cout << "temperature: " << Temperature << endl;
      cout << "acceptance ratio: " << accept_ratio << endl;
      cost(-1);
    }

    while (inner_loop_iter > 0) {
      cst = initiateMove();
      cost_hist.push_back(cst);

      wl_hist.push_back((wireLength() - wl_normalization.first)/(wl_normalization.second - wl_normalization.first));
      oa_hist.push_back((cellOverlap() - area_normalization.first)/(area_normalization.second - area_normalization.first));

      i -= 1;
    }

    // convergence criterion
    if(eps > 0 && abs(cost_hist.end()[-1] - cost_hist.end()[-2]) < eps) {
      break;
    }

    update_temperature();
    iii += 1;
    ii += 1;
    if (ii % 10 == 0) {
      writePlFile("./cache/"+std::to_string( ii )+".pl");
    }
  }
  gen_report(report);
  return cost();
}
