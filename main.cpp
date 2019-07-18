#include "main.h"
#define PI 3.14159265

using namespace std::chrono;

BOOST_GEOMETRY_REGISTER_BOOST_TUPLE_CS(cs::cartesian)

namespace bg = boost::geometry;
namespace bnu = boost::numeric::ublas;
namespace bgi = boost::geometry::index;
//typedef bg::model::point<float, 2, bg::cs::cartesian> point;
typedef bg::model::box< bg::model::d2::point_xy<double> > box;
typedef boost::geometry::model::polygon<boost::geometry::model::d2::point_xy<double> > polygon;

int debug = 1;
boundaries b;
boost::mt19937 rng;

//map < string, Node > nodeId;
vector < Node > nodeId;
map < string, int > name2id;
//bgi::rtree<std::pair<box, unsigned>, bgi::quadratic<16>> rtree;

float l1 = 1.0;
long long int iii = 0;

int main(int argc, char *argv[]) {
  srand(time(NULL));
  int opt;
  int outer_loop_iter = 105;
  int inner_loop_iter = 2;
  double eps = -1.0;
  double t_0 = 1.0;
  string out_file;

  long long int idx = -1;

  string parg = "";
  int bound_inc = 0;
  bool var = false;
  // TODO: replace with config
  while ((opt = getopt(argc,argv,"i:x:j:t:f:p:b::h")) != EOF) {
      switch(opt) {
          case 'x': idx = atoi(optarg); break;
          case 'p': parg=string(optarg); break;
          case 'd': debug=atoi(optarg); break;
          case 'i': outer_loop_iter = atoi(optarg); break;
          case 'j': inner_loop_iter = atoi(optarg); break;
          case 't': t_0 = atof(optarg); break;
          case 'e': eps = atof(optarg); break;
          case 'v': var = true; break;
          case 'f': out_file=string(optarg); break;
          case 'b':
            bound_inc += 1;
            switch(bound_inc) { // TODO: have 4 different argument flags for each point
              case 1: b.minX = atof(optarg);
              case 2: b.minY = atof(optarg);
              case 3: b.maxX = atof(optarg);
              case 4: b.maxY = atof(optarg);
            }
          case 'h': fprintf(stderr, "./sa usuage is \n \
                                    -i <optional, value> : for denoting # outer iterations PER SA INSTANCE \n \
                                    -j <optional, value> : for denoting 'j'*#nodes inner iterations \n \
                                    -t <optional, value> : for denoting initial temperature \n \
                                    -f <optional, str>   : for output filename \n \
                                    -e <optional, float> : convergence epsilon \n \
                                    -v <optional>        : ben-amur flag \n \
                                    -x <optional, int>   : simulated annealing instance index \n \
                                    -p <required, string>: input placement board \n \
                                    -d <optional, {0-3}> : debug verbosity \n \
                                    EXAMPLE: ./sa -i 20000 -j 20 -t 1 -p input -f output");
          default: cout<<endl; abort();
      }
  }
  if(parg == "")   {
    cout << "./main: option requires an argument -- p\n";
    exit(1);
  }
  if(parg == "")   {
    out_file = std::to_string( idx );
  }
  string p = string(parg);
  cout << "circuit: " << p << endl;

  string nodesfname  = p + ".nodes";
  string shapesfname = p + ".shapes";
  string netsfname   = p + ".nets";
  string plfname     = p + ".pl";
  string wtsfname    = p + ".wts";

  cout << nodesfname << endl;
  cout << shapesfname << endl;
  cout << netsfname << endl;
  cout << plfname << endl;

  map<int, vector<Pin> > netToCell;

  if(debug) { cout << "reading nodes..." << endl; }
  readNodesFile(nodesfname);
  //readShapesFile(shapesfname);
  //readWtsFile(wtsfname);
  if(debug) { cout << "reading pl..." << endl; }
  readPlFile(plfname);
  if(debug) { cout << "reading nets..." << endl; }
  netToCell = readNetsFile(netsfname);
  //readSclFile();
  if(debug) { cout << "calculating boundaries..." << endl; }
  set_boundaries();
  if(debug) { cout << "annealing" << endl; }
  float cost = timberWolfAlgorithm(outer_loop_iter, inner_loop_iter, eps, t_0, var, netToCell);
  writePlFile("./"+out_file+".pl");
  cout << " " << idx << " " << cost << endl;
  return cost;
}

/*
initialize_params
Empirically finds normalization parameters for scaling cost terms.
Currently, we scale by 1/(f) where f is the cost of the inital placement.
*/
void initialize_params(std::pair <double,double> *wl_normalization,
                       std::pair <double,double> *area_normalization,
                       std::pair <double,double> *routability_normalization,
                       map<int, vector<Pin> > &netToCell) {

  vector < std::pair <double,double> > normalization_terms;

  int num_components = 0;
  vector < Node > ::iterator itNode;
  for (itNode = nodeId.begin(); itNode != nodeId.end(); ++itNode) {
    if(!itNode -> terminal) {
      num_components += 1;
    }
  }

  // min/max wl
  std::pair < double, double > wl;
  double min_wl = std::numeric_limits<double>::max();
  double max_wl = 0.0;

  // min/max overlap area
  std::pair < double, double > area;
  double min_area = 0.0;
  double max_area = 0.0;

  vector < Node > ::iterator nodeit1;
  vector < Node > ::iterator nodeit2;

  for (nodeit1 = nodeId.begin(); nodeit1 != nodeId.end(); ++nodeit1) {
    if(nodeit1->terminal) {
      continue;
    } else {
        double cell_wl = nodeit1->width + nodeit1->height;
        if (cell_wl < min_wl) {
          min_wl = cell_wl;
        }
        for (nodeit2 = nodeit1++; nodeit2 != nodeId.end(); ++nodeit2) {
          if(nodeit2->idx == nodeit1->idx){continue;}
          if(nodeit1->fixed && nodeit2->fixed) {
              if(intersects(nodeit1->poly, nodeit2->poly)) {
                double oa = 0.0;
                std::deque<polygon> intersect_poly;
                boost::geometry::intersection(nodeit1->poly, nodeit2->poly, intersect_poly);

                BOOST_FOREACH(polygon const& p, intersect_poly) {
                    oa += bg::area(p);
                }
                min_area += pow(oa,2);
                max_area += pow(oa,2);
              }
            } else {
              max_area += pow(min(nodeit1->width, nodeit2->width) * min(nodeit1->height, (nodeit2->width)), 2);
            }
        }
      }
  }

  min_wl = min_wl*num_components;
  max_wl = ((b.maxX - b.minX) + (b.maxY - b.minY))*netToCell.size();

  wl.first = 0.0; // min_area
  wl.second = max(wirelength(netToCell),1.0); // max_area

  area.first = 0.0; // min_area;
  area.second = max(cell_overlap(),1.0); // max_area;

  normalization_terms.push_back(wl);
  normalization_terms.push_back(area);

  *wl_normalization = normalization_terms[0];
  *area_normalization = normalization_terms[1];

  routability_normalization->first = 0.0;
  routability_normalization->second = max(rudy(netToCell),1.0);
}

/*
set_boundaries
Calculate board boundaries by finding the maximum-area
rectangular envelop encompasing terminals & fixed modules.
*/
void set_boundaries() {
  double xval, yval, width, height;
  vector < Node > ::iterator itNode;
  for (itNode = nodeId.begin(); itNode != nodeId.end(); ++itNode) {
    if (itNode -> terminal || itNode->fixed) {
      xval = itNode -> xCoordinate;
      yval = itNode -> yCoordinate;
      width = itNode -> width;
      height =itNode -> height;
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

  if (b.minX == 0 && b.maxX == 0 && b.minY == 0 && b.maxY == 0) {
    cout << "no boundary found" << endl;
    b.maxX = 50;
    b.maxY = 50;
  }

  if(debug) {
    cout << "min: " << b.minX << " " << b.minY << endl;
    cout << "max: " << b.maxX << " " << b.maxY << endl;
  }
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

  boost::uniform_int<> uni_disto(0,7);
  boost::variate_generator<boost::mt19937&, boost::uniform_int<> > unio(rng, uni_disto);
  int ro = unio();

  string ostr = n.orient2str(ro);
  n.setParameterPl(rx,ry,ostr, n.fixed);
}

/*
initial_placement
Randomly place and orient all movable components in the board area
*/
void initial_placement() {
  vector < Node > ::iterator itNode;
  for (itNode = nodeId.begin(); itNode != nodeId.end(); ++itNode) {
    if (!itNode -> fixed) {
      randomPlacement(b.minX, b.maxX - itNode -> width, b.minY, b.maxY - itNode -> height, *itNode);
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
  double oa = cell_overlap();
  double eps = 10.0e-5;

  vector < Node > ::iterator nodeit1;
  vector < Node > ::iterator nodeit2;

  box env1;
  box env2;

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
              validate_move(&(*nodeit1),rx + dx,ry);
            } else {
              // project in y direction
              validate_move(&(*nodeit1),rx,ry + dy);
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
double wirelength(map<int, vector<Pin> > &netToCell) {
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
double cell_overlap() {
  double overlap = 0.0;

  for(int i = 0; i < nodeId.size(); i++) {
    if(nodeId[i].terminal) { continue; }
      for(int j = i++; j < nodeId.size(); j++) {
        if (i == j) {continue;}
        if(nodeId[j].terminal) { continue; }

        if(!intersects(nodeId[i].poly, nodeId[j].poly) || (nodeId[i].fixed && nodeId[j].fixed)) {
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
double wirelength_partial(vector < Node > &nodes, map<int, vector<Pin> > &netToCell) {
  map<int, vector < Pin > > ::iterator itNet;
  vector < Pin > ::iterator itCellList;
  double xVal, yVal, wireLength = 0;
  for (itNet = netToCell.begin(); itNet != netToCell.end(); ++itNet) {
    //if (nodes.find("f") == nodes.end()) {
    //  continue;
    //}
    double minXW = b.maxX, minYW = b.maxY, maxXW = b.minX, maxYW = b.minY;
    for (itCellList = itNet -> second.begin(); itCellList != itNet -> second.end(); ++itCellList) {
      if(itCellList->name == "") {
        continue;
      }
      int orient = nodeId[itCellList->idx].orientation;
      xVal = nodeId[itCellList->idx].xBy2;
      yVal = nodeId[itCellList->idx].yBy2;
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
cell_overlap_partial
Compute sum squared overlap for select components
*/
double cell_overlap_partial(vector < Node > &nodes) {
  double overlap = 0.0;
  vector < Node > ::iterator nodeit1;
  vector < Node > ::iterator nodeit2;

  for (nodeit1 = nodeId.begin(); nodeit1 != nodeId.end(); ++nodeit1) {
    if (!nodeit1->terminal) {
      for (nodeit2 = nodeit1++; nodeit2 != nodeId.end(); ++nodeit2) {
        if(nodeit2->idx == nodeit1->idx){continue;}
        if(!intersects(nodeit1->poly, nodeit2->poly) || (nodeit1->fixed && nodeit2->fixed)) {
          continue;
        } else {
          double oa = 0.0;
          std::deque<polygon> intersect_poly;
          boost::geometry::intersection(nodeit1->poly, nodeit2->poly, intersect_poly);

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
double rudy(map<int, vector<Pin> > &netToCell) {
  return 0;
  static bnu::matrix<double> D (static_cast<int>(abs(b.maxY)+abs(b.minY)+1), static_cast<int>(abs(b.maxX)+abs(b.minX)+1), 0);
  static bnu::matrix<double> D_route_sup (static_cast<int>(abs(b.maxY)+abs(b.minY)+1), static_cast<int>(abs(b.maxX)+abs(b.minX)+1), 1);

  D.clear();

  map<int, vector < Pin > > ::iterator itNet;
  vector < Pin > ::iterator itCellList;

  // For each net
  for (itNet = netToCell.begin(); itNet != netToCell.end(); ++itNet) {
    double xVal, yVal, hpwl = 0.0;
    double minXW = b.maxX, minYW = b.maxY, maxXW = b.minX, maxYW = b.minY;
    double rudy = 0.0;
    for (itCellList = itNet -> second.begin(); itCellList != itNet -> second.end(); ++itCellList) {
      if(itCellList->name == "") {
        continue;
      }
      int orient = nodeId[itCellList->idx].orientation;
      xVal = nodeId[itCellList->idx].xBy2;
      yVal = nodeId[itCellList->idx].yBy2;

      // compute pin position from orientation & offsets
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

      // set read_net
      for (unsigned i = xVal-5; i < xVal+5; ++ i) {
          for (unsigned j = yVal-5; j < yVal+5; ++ j) {
            D (i,j) += 5;
          }
      }
    }
    // now have boundary of net
    hpwl = (abs((maxXW - minXW)) + abs((maxYW - minYW)));
    rudy = hpwl / (max((maxXW - minXW)*(maxYW - minYW), 1.0)); // rudy density

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

  return r;
}

double cost(
            std::pair <double,double> &wl_normalization,
            std::pair <double,double> &area_normalization,
            std::pair <double,double> &routability_normalization,
            map<int, vector<Pin> > &netToCell,
            int temp_debug) {
  double l2 = 1 - l1;
  if(debug > 1 || temp_debug == -1) {
    cout << "wirelength: " << wirelength(netToCell) << endl;
    cout << "overlap: " << cell_overlap() << endl;
    cout << "l1: " << l1 << " l2: " << l2 << endl;
    cout << "wirelength_cost: " << l1*(wirelength(netToCell) - wl_normalization.first)/(wl_normalization.second - wl_normalization.first) << endl;
    cout << "overlap_cost: " << l2  * 0.9 * (cell_overlap() - area_normalization.first)/(area_normalization.second - area_normalization.first) << endl;
    cout << "routability_cost: " << l2  * 0.1  * rudy(netToCell) << endl;
    cout << "cost: " << l1*(wirelength(netToCell) - wl_normalization.first)/(wl_normalization.second - wl_normalization.first) +
                        l2  * 0.9 * (cell_overlap() - area_normalization.first)/(area_normalization.second - area_normalization.first)  +
                        l2 * 0.1*rudy(netToCell) << endl;
  }
  return l1 * (wirelength(netToCell) - wl_normalization.first)/(wl_normalization.second - wl_normalization.first) +
         l2 * 0.9 * (cell_overlap() - area_normalization.first)/(area_normalization.second - area_normalization.first) +
         l2 * 0.1 * rudy(netToCell);
         //l2 * 0.1 * (rudy(netToCell) - routability_normalization.first)/(routability_normalization.second - routability_normalization.first);
}

double cost_partial(int temp_debug,
                    vector < Node > nodes,
                    std::pair <double,double> &wl_normalization,
                    std::pair <double,double> &area_normalization,
                    std::pair <double,double> &routability_normalization,
                    map<int, vector<Pin> > &netToCell) {
  double l2 = 1-l1;
  return l1 * (wirelength_partial(nodes, netToCell) - wl_normalization.first)/(wl_normalization.second - wl_normalization.first) +
         l2 * 0.9 * (cell_overlap_partial(nodes) - area_normalization.first)/(area_normalization.second - area_normalization.first) +
		     l2 * 0.1 * rudy(netToCell);
}

/*
random_node
Select random node from set of nodes
*/
vector < Node >::iterator random_node() {
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
void validate_move(Node* node, double rx, double ry) {
  int orient = node->orientation;
  double width = 0.0;
  double height = 0.0;
  if (orient == 0 || orient == 4) {
    width = node->width;
    height = node->height;
  } else if (orient == 2 || orient == 6) {
    width = node->height;
    height = node->width;
  }else {
    width = max(node->width,node->height);
    height = max(node->width,node->height);
  }
  rx = max(rx, b.minX);
  ry = max(ry, b.minY);

  rx = min(rx, b.maxX - width);
  ry = min(ry, b.maxY - height);

  node->setPos(rx,ry);
}

double c = 0.0;
double initiate_move(vector< int > *accept_history,
                     double & Temperature,
                     std::pair <double,double> &wl_normalization,
                     std::pair <double,double> &area_normalization,
                     std::pair <double,double> &routability_normalization,
                     map<int, vector<Pin> > &netToCell) {
  // Initate a transition
  int state = -1;
  double prevCost = cost(wl_normalization, area_normalization, routability_normalization, netToCell);
  if(debug > 1) {
    cout << "prevCost: " << prevCost << endl;
  }
  vector < Node > ::iterator rand_node1;
  vector < Node > ::iterator rand_node2;
  int r = 0;
  rand_node1 = random_node();
  while(rand_node1->terminal || rand_node1->name == "" || rand_node1->fixed) {
    rand_node1 = random_node();
  }
  //rtree.remove(std::make_pair(rand_node1->envelope, rand_node1->idx));
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

  if (i<0.20) { // swap
    state = 0;
    if(debug > 1) {
      cout << "swap" << endl;
    }
    rand_node2 = random_node();
    while(rand_node2->terminal || rand_node2->idx == rand_node1->idx || rand_node2->name == "" || rand_node2->fixed) {
      rand_node2 = random_node();
    }
    //rtree.remove(std::make_pair(rand_node2->envelope, rand_node2->idx));
    rand_node2_orig_x = rand_node2->xCoordinate;
    rand_node2_orig_y = rand_node2->yCoordinate;

    validate_move(&(*rand_node1), rand_node2_orig_x, rand_node2_orig_y);
    validate_move(&(*rand_node2), rand_node1_orig_x, rand_node1_orig_y);
    //rtree.insert(std::make_pair(rand_node1->envelope, rand_node1->idx));
    //rtree.insert(std::make_pair(rand_node2->envelope, rand_node2->idx));
  } else if (i < 0.85) { // shift
    state = 1;
    if(debug > 1) {
      cout << "shift" << endl;
    }
    double sigma = rand_node1->sigma;

    boost::normal_distribution<> nd(0.0, sigma);
    boost::variate_generator<boost::mt19937&,
                             boost::normal_distribution<> > var_nor(rng, nd);

    double dx = var_nor();
    double dy = var_nor();

    double rx = rand_node1_orig_x + dx;
    double ry = rand_node1_orig_y + dy;

    validate_move(&(*rand_node1), rx, ry);
    //rtree.insert(std::make_pair(rand_node1->envelope, rand_node1->idx));
  } else { // rotate
    state = 2;
    if(debug > 1) {
      cout << "rotate" << endl;
    }
    boost::uniform_int<> uni_dist(0,7);
    boost::variate_generator<boost::mt19937&, boost::uniform_int<> > uni(rng, uni_dist);
    r = uni();
    rand_node1->setRotation(r);
    validate_move(&(*rand_node1), rand_node1_orig_x, rand_node1_orig_y);
    //rtree.insert(std::make_pair(rand_node1->envelope, rand_node1->idx));
  }
  bool accept = check_move(prevCost,
                           accept_history,
                           Temperature,
                           wl_normalization,
                           area_normalization,
                           routability_normalization,
                           netToCell);
  if (!accept) {
    if(debug > 1) {
      cout << "reject" << endl;
    }
    // revert state
    //rtree.remove(std::make_pair(rand_node1->envelope, rand_node1->idx));
    if (state == 0) {
      //rtree.remove(std::make_pair(rand_node2->envelope, rand_node2->idx));
      rand_node1->setPos(rand_node1_orig_x,rand_node1_orig_y);
      rand_node2->setPos(rand_node2_orig_x,rand_node2_orig_y);
      //rtree.insert(std::make_pair(rand_node2->envelope, rand_node2->idx));
    } else if (state == 1) {
      rand_node1->setPos(rand_node1_orig_x,rand_node1_orig_y);
    } else if (state == 2) {
      rand_node1->setRotation(8 - r);
      rand_node1->setPos(rand_node1_orig_x,rand_node1_orig_y);
    }
    //rtree.insert(std::make_pair(rand_node1->envelope, rand_node1->idx));
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
void update_temperature(double* Temperature) {
  vector < Node > ::iterator nodeit = nodeId.begin();
  if (*Temperature > 0.1) {
    *Temperature = (0.985) * *Temperature;
    if (l1 > 92e-2) {
      l1 -= 2e-4;
    }
    for (nodeit = nodeId.begin(); nodeit != nodeId.end(); ++nodeit) {
      nodeit->sigma =  max(0.9999*nodeit->sigma,0.5);
    }
  } else if (*Temperature > 0.01) {
    *Temperature = (0.9995) * *Temperature;
    if (l1 > 88.5e-2) {
      l1 -= 4e-4;
    }
    for (nodeit = nodeId.begin(); nodeit != nodeId.end(); ++nodeit) {
      nodeit->sigma =  max(0.999*nodeit->sigma,0.25);
    }
  } else if (*Temperature > 0.005) {
    *Temperature = (0.9955) * *Temperature;
    if (l1 > 85.5e-2) {
      l1 -= 8e-4;
    }
    for (nodeit = nodeId.begin(); nodeit != nodeId.end(); ++nodeit) {
      nodeit->sigma =  max(0.999*nodeit->sigma,0.18);
    }
  } else if (*Temperature > 0.001) {
    *Temperature = (0.9965) * *Temperature;
    if (l1 > 82e-2) {
      l1 -= 16e-4;
    }
    for (nodeit = nodeId.begin(); nodeit != nodeId.end(); ++nodeit) {
      nodeit->sigma =  max(0.999*nodeit->sigma,0.1);
    }
  } else {
    if (*Temperature > 0.000001) {
      *Temperature = (0.855) * *Temperature;
    } else {
      *Temperature = 0.0000000001;
    }
    if (l1 > 1e-4) {
      l1 -= 10e-4;
    }
    for (nodeit = nodeId.begin(); nodeit != nodeId.end(); ++nodeit) {
      nodeit->sigma =  max(0.999*nodeit->sigma,0.1);
    }
  }
  l1 = 0.0;
}

void update_accept_history(vector< int > &accept_history, vector< double > *accept_ratio_history, float *accept_ratio) {
  if (accept_history.size() > 100) {
    *accept_ratio = ((*accept_ratio*100) - (accept_history[accept_history.size()-101]) + accept_history[accept_history.size()-1])/100;
  } else {
    *accept_ratio = accumulate(accept_history.begin(), accept_history.end(),0.0) / accept_history.size();
  }
  accept_ratio_history->push_back(*accept_ratio);
}

/*
check_move
either accept or reject the move based on current & previous temperature & cost
*/
bool check_move(double prevCost,
                vector< int > *accept_history,
                double & Temperature,
                std::pair <double,double> &wl_normalization,
                std::pair <double,double> &area_normalization,
                std::pair <double,double> &routability_normalization,
                map<int, vector<Pin> > &netToCell) {
  double newCost = cost(wl_normalization, area_normalization, routability_normalization, netToCell);
  c = newCost;
  double delCost = 0;
  boost::uniform_real<> uni_dist(0,1);
  boost::variate_generator<boost::mt19937&, boost::uniform_real<> > uni(rng, uni_dist);
  double prob = uni();
  if(debug > 1) {
    cout << "new cost: " << newCost << endl;
  }
  delCost = newCost - prevCost;
  if (delCost <= 0  || prob <= (exp(-delCost/Temperature))) {
    prevCost = newCost;
    accept_history->push_back(1);
    return true;
  } else {
    accept_history->push_back(0);
    return false;
  }
}

double initialize_temperature(vector< int > &accept_history,
                              double & Temperature,
                              std::pair <double,double> &wl_normalization,
                              std::pair <double,double> &area_normalization,
                              std::pair <double,double> &routability_normalization,
                              map<int, vector<Pin> > &netToCell) {
  double t = 0.0;
  double emax = 0.0;
  double emin = 0.0;
  double xt = 1.0;
  double x0 = 0.84;
  double p = 2.0;
  initial_placement();
  for(int i=1; i<=10; i++){
    for(int j=1; j<=10; j++){
      initial_placement();
      emax += exp(cost(wl_normalization, area_normalization, routability_normalization, netToCell)/t);
      initiate_move(&accept_history, Temperature, wl_normalization, area_normalization, routability_normalization, netToCell);
      emin += exp(cost(wl_normalization, area_normalization, routability_normalization, netToCell)/t);
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
void gen_report(map<string, vector<double> > &report,
                vector< double > &accept_ratio_history,
                std::pair <double,double> &wl_normalization,
                std::pair <double,double> &area_normalization,
                std::pair <double,double> &routability_normalization,
                map<int, vector<Pin> > &netToCell) {
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
float timberWolfAlgorithm(int outer_loop_iter,
                          int inner_loop_iter,
                          double eps,
                          double t_0,
                          bool var,
                          map<int, vector<Pin> > &netToCell) {
  double Temperature = t_0;
  int num_components = 0;

  vector < Node > ::iterator itNode;
  map < string, vector < double > > report;
  vector < double > cost_hist;
  vector < double > wl_hist;
  vector < double > oa_hist;
  report["cost_hist"] = cost_hist;
  report["wl_hist"] = wl_hist;
  report["oa_hist"] = oa_hist;

  std::pair <double,double> wl_normalization;
  std::pair <double,double> area_normalization;
  std::pair <double,double> routability_normalization;

  vector< int > accept_history;
  float accept_ratio = 0;
  vector< double > accept_ratio_history;
  initial_placement();
  initialize_params(&wl_normalization, &area_normalization, &routability_normalization, netToCell);
  double cst = cost(wl_normalization, area_normalization, routability_normalization, netToCell);
  if(var) {
    Temperature = initialize_temperature(accept_history,
                                         Temperature,
                                         wl_normalization,
                                         area_normalization,
                                         routability_normalization,
                                         netToCell);
  }
  int idx = 0;
  for (itNode = nodeId.begin(); itNode != nodeId.end(); ++itNode) {
    //rtree.insert(std::make_pair(itNode -> envelope, idx));
    idx+=1;
    if(!itNode -> terminal && !itNode -> fixed) {
      num_components += 1;
    }
  }

  high_resolution_clock::time_point t1 = high_resolution_clock::now();
  long long int ii = 0; // outer loop iterator
  int i = 0; // inner loop iterator
  while (ii < outer_loop_iter) {
    i = inner_loop_iter*num_components;
    if(ii > 1 && ii % 10 == 0 && debug) {
      high_resolution_clock::time_point t2 = high_resolution_clock::now();
      duration<double> time_span = duration_cast< duration<double> >(t2 - t1);

      cout << "******" << ii << "******" << endl;
      cout << "iteration: " << ii << endl;
      cout << "time: " <<  time_span.count() << " (s)" << endl;
      cout << "move/time: " <<  i*ii/time_span.count() << endl;
      cout << "time remaining: " <<  time_span.count()/ii * (1000-ii) << " (s)" << endl;
      cout << "temperature: " << Temperature << endl;
      cout << "acceptance ratio: " << accept_ratio << endl;
      cost(wl_normalization, area_normalization, routability_normalization,netToCell,-1);
    }
    while (i > 0) {
      cst = initiate_move(&accept_history, Temperature, wl_normalization, area_normalization, routability_normalization, netToCell);
      update_accept_history(accept_history, &accept_ratio_history, &accept_ratio);
      cost_hist.push_back(cst);

      wl_hist.push_back((wirelength(netToCell) - wl_normalization.first)/(wl_normalization.second - wl_normalization.first));
      oa_hist.push_back((cell_overlap() - area_normalization.first)/(area_normalization.second - area_normalization.first));

      i -= 1;
    }

    // convergence criterion
    if(eps > 0 && abs(cost_hist.end()[-1] - cost_hist.end()[-2]) < eps) {
      break;
    }
    update_temperature(&Temperature);
    iii += 1;
    ii += 1;
    if (ii % 10 == 0) {
      writePlFile("./cache/"+std::to_string( ii )+".pl");
    }
  }
  gen_report(report,
             accept_ratio_history,
             wl_normalization,
             area_normalization,
             routability_normalization,
             netToCell);

  return cost(wl_normalization, area_normalization, routability_normalization, netToCell);
}
