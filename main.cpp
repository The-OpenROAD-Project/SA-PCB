#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctype.h>
#include <map>
#include <vector>
#include <string>
#include <numeric>
#include <cmath>
#include <unistd.h>
#include <ctime>
#include <ratio>
#include <chrono>

// boost version 1.69.0
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/geometry.hpp>
#include <boost/assign.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/io/io.hpp>
#include <boost/geometry/algorithms/area.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/foreach.hpp>

//#include "taskflow/taskflow.hpp"

#include "main.h"
#include "time.h"
#include "readScl.h"

using namespace std::chrono;

BOOST_GEOMETRY_REGISTER_BOOST_TUPLE_CS(cs::cartesian)

namespace bg = boost::geometry;
typedef boost::geometry::model::polygon<boost::geometry::model::d2::point_xy<double> > polygon;

double Temperature;
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

int debug = 0;
long long int idx = -1;

vector< int > accept_history;
float accept_ratio = 0;

boost::mt19937 rng;

int main(int argc, char *argv[]) {
  // parser and entry point into the program
  // TODO add wirelength flag
  srand(time(NULL));
  int opt;
  int i, j, t, k, s, w;
  i=j=t=k=1;
  string f = "";
  char *parg;
  /*while ((opt = getopt(argc,argv,"i:x:j:t:k:s:w:f:h:p::d")) != EOF) {
      switch(opt) {
          case 'p': parg=strdup(optarg); break;
          case 'd': debug=1; break;
          case 'i': i = 1; break;
          case 'j': j = 1; break;
          case 't': t = 1; break;
          case 'k': k = 1; break;
          case 's': s = 1; break;
          case 'w': w = 1; break;
          case 'f': f = ""; break;
          case 'x': idx = atoi(optarg); break;
          case 'h': fprintf(stderr, "./sa usuage is \n \
                                    -i <value> : for denoting # outer iterations PER SA INSTANCE \n \
                                    -j <value> : for denoting 'j'*#nodes inner iterations \n \
                                    -t <value> : for denoting initial temperature \n \
                                    -k <value> : for denoting parallel instances \n \
                                    -s <value> : for denoting splits every floor(i/s) iterations. Use s=0 for 'k' independent instances of SA \n \
                                    -w <value> : for denoting number of winners per splitting \n \
                                    -f <str> : for denoting number output pl \n \
                                    EXAMPLE: ./sa -i 20000 -j 20 -t 40000 -k 0 -s -w 0 -f output.pl for the standard single instance timberwolf algorithm");
          default: cout<<endl; abort();
      }*/
      /*if (idx == -1) {
        cout << "./main: option requires an argument -- x\n";
        exit(1);
      }
      if (i == -1) {
         cout << "./main: option requires an argument -- i\n";
         exit(1);
      }
      if (j == -1) {
         cout << "./main: option requires an argument -- j\n";
         exit(1);
      }
      if (t == -1) {
         cout << "./main: option requires an argument -- t\n";
         exit(1);
      }
      if (k == -1) {
         cout << "./main: option requires an argument -- k\n";
         exit(1);
      }
      if (f == "") {
        cout << "./main: option requires an argument -- f\n";
        exit(1);
      }*/
      //if (parg == "") {
      /*if(strncmp(parg,"",1) == 0) {
        cout << "./main: option requires an argument -- p\n";
        exit(1);
      }
    }

  // optind is for the extra arguments which are not parsed
  for(; optind < argc; optind++){
    printf("extra arguments: %s\n", argv[optind]);
  }*/
  cout << "idx: " << idx << endl;
  //std::string p = std::string(parg);
  string p = "bm3";
  cout << "circuit: " << p << endl;

  debug = 1;

  string nodesfname = p +".nodes";
  string netsfname = p +".nets";
  string plfname = p +".pl";
  string wtsfname = p +".wts";

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
  SetInitParameters();
  if(debug) { cout << "annealing" << endl; }
  float cost = timberWolfAlgorithm();
  writePlFile("./"+std::to_string( idx )+".pl");
  return cost;
}

//vector < std::pair <double,double> > SetInitParameters() { // (x - min) / (max - min)
void SetInitParameters() { // (x - min) / (max - min)
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
            max_area += pow(min(nodeit1->second.width, nodeit2->second.width) * min(nodeit1->second.height, (nodeit2->second.width)), 2);
          }
        }
      }
  }
  min_wl = min_wl*num_components;
  max_wl = ((b.maxX - b.minX) + (b.maxY - b.minY))*netToCell.size();
  wl.first = min_wl;
  wl.second = max_wl;

  area.first = min_area;
  area.second = max_area;

  normalization_terms.push_back(wl);
  normalization_terms.push_back(area);

  wl_normalization = normalization_terms[0];
  area_normalization = normalization_terms[1];

  //return normalization_terms;
}

void CalcBoundaries() {
  // calculate boundaries from terminal nodes
  double xval, yval;
  map < string, Node > ::iterator itNode;
  for (itNode = nodeId.begin(); itNode != nodeId.end(); ++itNode) {
    if (itNode -> second.terminal) {
      xval = itNode -> second.xCoordinate;
      yval = itNode -> second.yCoordinate;
      if (xval < b.minX) {
        b.minX = xval;
      }
      if (xval > b.maxX) {
        b.maxX = xval;
      }
      if (xval < b.minY) {
        b.minY = yval;
      }
      if (xval > b.maxY) {
        b.maxY = yval;
      }
    }
  }
  b.minX = -10.0;
  b.maxX = 100.0;
  b.minY = -15.0;
  b.maxY = 60.0;
  if(debug) {
    cout << "min: " << b.minX << " " << b.minY << endl;
    cout << "max: " << b.maxX << " " << b.maxY << endl;
  }
}

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

void randomPlacement(int xmin, int xmax, int ymin, int ymax, Node n) {
  // randomly place node
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
  n.setParameterPl(rx,ry,ostr , n.fixed);
}

void initialPlacement() {
  // initially randomly place components on the board
  map < string, Node > ::iterator itNode;
  for (itNode = nodeId.begin(); itNode != nodeId.end(); ++itNode) {
    Node n = itNode -> second;
    if (!n.fixed) {
      randomPlacement(b.minX, b.maxX - n.width, b.minY, b.maxY - n.height, n);
    }
  }
}

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
    wireLength += (abs((maxXW - minXW)) + abs((maxYW - minYW)))/max(static_cast<int>(itNet -> second.size() - 1) ,1);
  }
  return wireLength;
}

double cellOverlap() {
  // compute sum squared overlap
  double overlap = 0.0;
  map < string, Node > ::iterator nodeit = nodeId.begin();
  map < string, Node > ::iterator nodeit2 = nodeId.begin();

  vector < string > hist;
  for (nodeit = nodeId.begin(); nodeit != nodeId.end(); ++nodeit) {
    hist.push_back(nodeit->second.name);
    if (!nodeit->second.terminal) {
      for (nodeit2 = nodeId.begin(); nodeit2 != nodeId.end(); ++nodeit2) {
        if (!nodeit2->second.terminal && !(find(hist.begin(), hist.end(), nodeit2->second.name) != hist.end()) ) {
          if(!intersects(nodeit->second.poly, nodeit2->second.poly)) {
            continue;
          } else {
            std::deque<polygon> intersect_poly;
            boost::geometry::intersection(nodeit->second.poly, nodeit2->second.poly, intersect_poly);

            BOOST_FOREACH(polygon const& p, intersect_poly) {
                overlap +=  pow(bg::area(p),2);
            }
          }
        }
      }
    }
  }
  return overlap;
}

double rudy() {
  // compute magnitude of routability matrix
  map<int, vector < Pin > > ::iterator itNet;
  vector < Pin > ::iterator itCellList;
  double rudy = 0;
  for (itNet = netToCell.begin(); itNet != netToCell.end(); ++itNet) {
    double xVal, yVal, wireLength = 0;
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
      if (xVal > maxXW) {
        if (xVal > 10000) {
        cout << maxXW << " " << xVal << endl;
        nodeId[itCellList->name].printParameter();}
        maxXW = xVal;}
      if (yVal < minYW)
        minYW = yVal;
      if (yVal > maxYW)
        maxYW = yVal;
    }
    wireLength = (abs((maxXW - minXW)) + abs((maxYW - minYW)))/(itNet -> second.size() - 1);
    rudy += wireLength / ((maxXW - minXW)*(maxYW - minYW));
  }
  return rudy;
}

double cost(int temp_debug) {
  if(debug > 1 || temp_debug == -1) {
    cout << "wirelength: " << wireLength() << endl;
    cout << "overlap: " << cellOverlap() << endl;
    cout << "l1: " << l1 << " l2: " << l2 << endl;
    cout << "wirelength_cost: " << l1*(wireLength() - wl_normalization.first)/(wl_normalization.second - wl_normalization.first) << endl;
    cout << "overlap_cost: " << l2* (cellOverlap() - area_normalization.first)/(area_normalization.second - area_normalization.first) << endl;
    cout << "cost: " << l1*(wireLength() - wl_normalization.first)/(wl_normalization.second - wl_normalization.first) +
                        l2* (cellOverlap() - area_normalization.first)/(area_normalization.second - area_normalization.first) << endl;
  }
  return l1*(wireLength() - wl_normalization.first)/(wl_normalization.second - wl_normalization.first) +
         l2* (cellOverlap() - area_normalization.first)/(area_normalization.second - area_normalization.first);
}

double sigma = 40.0;
map < string, Node > ::iterator random_node() {
  map < string, Node > ::iterator itNode = nodeId.begin();
  int size = nodeId.size();
  boost::uniform_int<> uni_dist(0,size-1);
  boost::variate_generator<boost::mt19937&, boost::uniform_int<> > uni(rng, uni_dist);
  int randint = uni();
  std::advance(itNode, randint);
  return itNode;
}

void validateMove(Node* node, double rx, double ry) {
  double width = node->width;
  double height = node->height;

  rx = max(rx, b.minX);
  ry = max(ry, b.minX);

  rx = min(rx, b.maxX - width);
  ry = min(ry, b.maxY - height);

  node->setPos(rx,ry);
}

void initiateMove() {
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
  }
  else {
    if(debug > 1) {
      cout << "accept" << endl;
    }
  }
}

void update_temperature() {
  // update the temperature according to annealing schedule
  if (Temperature > 50000) {
    Temperature = (0.9992) * Temperature;
  } else if (Temperature <= 50000 && Temperature > 5000) {
    Temperature = (0.9992) * Temperature;
  } else if (Temperature < 5000) {
    Temperature = (0.999) * Temperature;
  } else if (Temperature < 1) {
    Temperature = (0.999) * Temperature;
  }
  if (l1 > 0.01) {
    l1 -= 0.00005;
    l2 += 0.00005;
  }
  sigma = max(sigma  - 0.002,0.01);
}

bool checkMove(double prevCost) {
  // either accept or reject the move based on temperature & cost
  double newCost = cost();
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

float timberWolfAlgorithm() {
  // main driver of the anealer
  //xLimit = macroPlacement();

  //varanelli_cohoon();

  initialPlacement();

  //Temperature = pow(8,5);
  Temperature = 1.0;
  int i; // n * 20, n is number of nodes
  int cnt = 0;
  long long int ii = 0;

  int num_components = 0;
  map < string, Node > ::iterator itNode;
  for (itNode = nodeId.begin(); itNode != nodeId.end(); ++itNode) {
    if(!itNode -> second.terminal && !itNode -> second.fixed) {
      num_components += 1;
    }
  }

  high_resolution_clock::time_point t1 = high_resolution_clock::now();
  //while (Temperature > 0.1) {
  while (ii < 25000) {
    //i = num_components;
    i = 20;
    if(ii % 100 == 0 && debug) {
      high_resolution_clock::time_point t2 = high_resolution_clock::now();
      duration<double> time_span = duration_cast< duration<double> >(t2 - t1);

      cout << "******" << ii << "******" << endl;
      cout << "iteration: " << ii << endl;
      cout << "time: " <<  time_span.count() << " (s)" << endl;
      cout << "move/time: " <<  i*ii/time_span.count() << endl;
      cout << "time remaining: " <<  time_span.count()/ii * (20000-ii) << " (s)" << endl;
      cout << "temperature: " << Temperature << endl;
      cout << "acceptance ratio: " << accept_ratio << endl;
      cost(-1);
    }

    while (i > 0) {
      initiateMove();
      i -= 1;
      cnt += 1;
    }
    update_temperature();
    ii += 1;
    if (ii % 10 == 0) {
      writePlFile("./cache/"+std::to_string( ii )+".pl");
    }
  }
  writePlFile("./cache/"+std::to_string( idx )+".pl");
  return cost();
}
/*
float multistart() {
  auto [A, B, C, D] = taskflow.emplace(
    [] () { timberWolfAlgorithm(); },
    [] () { timberWolfAlgorithm(); },
    [] () { timberWolfAlgorithm(); },
    [] () { timberWolfAlgorithm(); }
  );

  taskflow.wait_for_all();  // block until finish
  return 0;
}*/

void printCellMap() {
  map < string, Node > ::iterator itNode;
  for (itNode = nodeId.begin(); itNode != nodeId.end(); ++itNode) {
    itNode -> second.printParameter();
  }
}

void printRowMap() {
  map < int, row > ::iterator itRow;
  for (itRow = rowId.begin(); itRow != rowId.end(); ++itRow) {
    itRow -> second.printParameter();
  }
}
