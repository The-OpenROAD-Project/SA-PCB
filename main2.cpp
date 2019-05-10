#include <iostream>
#include <vector>
#include <fstream>
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
#include <boost/foreach.hpp>
#include <string>
#include <iostream>
#include <numeric>
#include <sstream>
#include <regex>
#include <map>
#include <cmath>
#include "main.h"
#include "time.h"
//#include "readNets.h"
#include "readScl.h"

#include <ctime>
#include <ratio>
#include <chrono>

using namespace std::chrono;

/*
TODO: accept/reject ratio - differen parameter
*/

BOOST_GEOMETRY_REGISTER_BOOST_TUPLE_CS(cs::cartesian)

namespace bg = boost::geometry;
typedef boost::geometry::model::polygon<boost::geometry::model::d2::point_xy<double> > polygon;

double Temperature;
int xLimit;
map < int, row > rowId;
int RowWidth;
map<int, vector<string> > netToCell;
map < string, Node > nodeId;
boundaries b;

int debug = 0;

vector< int > accept_history;
float accept_ratio = 0;

int main(int argc, char *argv[]) {
  // parser and entry point into the program
  // TODO add wirelength flag
  int idx = -1;
  int opt;
  int i, j, t, k, s, w;
  i=j=t=k=-1;
  string f = "";
  while ((opt = getopt(argc,argv,"i:j:t:k:s:w:f:h")) != EOF) {
      switch(opt) {
          case 'i': i = 1; cout << optarg <<endl; break;
          case 'j': j = 1; cout << optarg <<endl; break;
          case 't': t = 1; cout << optarg <<endl; break;
          case 'k': k = 1; cout << optarg <<endl; break;
          case 's': s = 1; cout << optarg <<endl; break;
          case 'w': w = 1; cout << optarg <<endl; break;
          case 'f': f = ""; cout << optarg <<endl; break;
          case 'x': idx = 1; cout << optarg <<endl; break;
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
      }
      if (idx == -1) {
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
      }
    }

  // optind is for the extra arguments which are not parsed
  for(; optind < argc; optind++){
    printf("extra arguments: %s\n", argv[optind]);
  }

  readNodesFile();
  //readWtsFile();
  readPlFile();
  readNetsFile();
  //readSclFile();
  CalcBoundaries();
  float cost = timberWolfAlgorithm();
  cout << cost << endl;
  return cost;
}

void CalcBoundaries() {
  // calculate boundaries from terminal nodes
  int xval, yval;
  map < string, Node > ::iterator itNode;
  for (itNode = nodeId.begin(); itNode != nodeId.end(); ++itNode) {
    if (itNode -> second.terminal == 0) {
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
}

int macroPlacement() {
  // macro placement [depreceiated]
  int xValue = b.minX, yValue = 0;
  xLimit = b.minX;
  map < int, row > ::iterator itRow;
  itRow = rowId.begin();
  int rowHeight = itRow -> second.height;
  for (map < string, Node > ::iterator itNode = nodeId.begin(); itNode != nodeId.end(); ++itNode) {
    if (itNode -> second.terminal == 1 && itNode -> second.height > rowHeight) {
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
  int nx = n.xCoordinate;
  int ny = n.yCoordinate;
  int rx = rand()%(xmax-xmin + 1) + xmin;
  int ry = rand()%(ymax-ymin + 1) + ymin;

  int dx = rx - nx;
  int dy = ry - ny;
  int ro = rand()%4;
  string ostr = n.orient2str(ro);
  n.setParameterPl(nx + dx, ny +dy,ostr , n.fixed);
}

void initialPlacement() {
  // initially randomly place components on the board
  map < string, Node > ::iterator itNode;
  for (itNode = nodeId.begin(); itNode != nodeId.end(); ++itNode) {
    Node n = itNode -> second;
    if (n.fixed == 0) {
      randomPlacement(b.minX, b.maxX - n.width, b.minY, b.maxY - n.height, n);
    }
    //n.printParameter();
  }
}

double wireLength() {
  // compute HPWL
  map < int, vector < string > > ::iterator itNet;
  vector < string > ::iterator itCellList;
  double xVal, yVal, wireLength = 0;
  int minXW = b.minX, minYW = b.minY, maxXW = b.maxX, maxYW = b.maxY;
  for (itNet = netToCell.begin(); itNet != netToCell.end(); ++itNet) {
    for (itCellList = itNet -> second.begin(); itCellList != itNet -> second.end(); ++itCellList) {
      xVal = nodeId[ * itCellList].xBy2;
      yVal = nodeId[ * itCellList].yBy2;
      if (xVal < minXW)
        minXW = xVal;
      if (xVal > maxXW)
        maxXW = xVal;
      if (yVal < minYW)
        minYW = yVal;
      if (yVal > maxYW)
        maxYW = yVal;
    }
    wireLength += abs((maxXW - minXW)) + abs((maxYW - minYW));
  }
  return wireLength;
}

double cellOverlap() {
  // compute sum squared overlap
  double overlap = 0;
  map < string, Node > ::iterator nodeit = nodeId.begin();
  map < string, Node > ::iterator nodeit2 = nodeId.begin();
  for (nodeit = nodeId.begin(); nodeit != nodeId.end(); ++nodeit) {
    if (nodeit->second.terminal == 1) {
      for (nodeit2 = nodeId.begin(); nodeit2 != nodeId.end(); ++nodeit2) {
        if (nodeit2->second.terminal == 1) {
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

double cost(int temp_debug) {
  //double wl = wireLength();
  //double oa = cellOverlap();
  //float cost = wl + oa;
  double cost = wireLength() + cellOverlap();
  if(debug || temp_debug == -1) {
    cout << "wirelength: " << wireLength() << endl;
    cout << "overlap: " << cellOverlap() << endl;
    cout << "cost: " << wireLength() + cellOverlap() << endl;
  }

  return wireLength() + cellOverlap();
}

map < string, Node > ::iterator random_node() {
  map < string, Node > ::iterator itNode = nodeId.begin();
  int size = nodeId.size();
  int randint = rand() % static_cast<int>(size + 1);
  std::advance(itNode, randint);
  return itNode;
}

void initiateMove() {
  // Initate a transition
  int state = -1;
  double prevCost = cost();
  if(debug) {
    cout << "prevCost: " << prevCost << endl;
  }
  map < string, Node > ::iterator rand_node1;
  int orig_x = rand_node1->second.xCoordinate;
  int orig_y =rand_node1->second.yCoordinate;
  map < string, Node > ::iterator rand_node2;
  int r = 0;
  rand_node1 = random_node();
  while(rand_node1->second.terminal != 1 || rand_node1->second.name == "") {
    rand_node1 = random_node();
  }
  if(debug) {
    cout << "=======" << endl;
    cout <<  "name: " << rand_node1->second.name << endl;
  }

  float i = ((double) rand() / (RAND_MAX));
  if (i<0.4) { // swap
    state = 0;
    if(debug) {
      cout << "swap" << endl;
    }
    rand_node2 = random_node();
    while(rand_node2->second.terminal != 1 || rand_node2->first == rand_node1->first || rand_node1->second.name == "") {
      rand_node2 = random_node();
    }
    int rand_node1_x = rand_node1->second.xCoordinate;
    int rand_node1_y = rand_node1->second.yCoordinate;
    int rand_node2_x = rand_node2->second.xCoordinate;
    int rand_node2_y = rand_node2->second.yCoordinate;

    rand_node1->second.setPos(rand_node2_x,rand_node2_y);
    rand_node2->second.setPos(rand_node1_x,rand_node1_y);

  } else if (i < 0.8) { // shift
    state = 1;
    if(debug) {
      cout << "shift" << endl;
    }
    int rx = b.minX + (rand() % static_cast<int>(b.maxX - b.minX + 1));
    int ry = b.minY + (rand() % static_cast<int>(b.maxY - b.minY + 1));
    rand_node1->second.setPos(rx,ry);
  } else if (i < 0) { // rotate
    state = 2;
    if(debug) {
      cout << "rotate" << endl;
    }
    r = (rand() % static_cast<int>(3 + 1));
    rand_node1->second.setRotation(r);
  }

  bool accept = checkMove(prevCost);
  if (!accept) {
    if(debug) {
      cout << "reject" << endl;
    }
    // revert state
    if (state == 0) {
      int rand_node1_x = rand_node1->second.xCoordinate;
      int rand_node1_y = rand_node1->second.yCoordinate;
      int rand_node2_x = rand_node2->second.xCoordinate;
      int rand_node2_y = rand_node2->second.yCoordinate;

      rand_node1->second.setPos(rand_node2_x,rand_node2_y);
      rand_node2->second.setPos(rand_node1_x,rand_node1_y);
    } else if (state == 1) {
      rand_node1->second.setPos(orig_x,orig_y);
    } else if (state == 2) {
      rand_node1->second.setRotation(-r);
    }
  }
  else {
    if(debug) {
      cout << "accept" << endl;
    }
  }
}

void update_temperature() {
  // update the temperature according to annealing schedule
  if (Temperature > 50000) {
    Temperature = 0.8 * Temperature;
  } else if (Temperature <= 50000 && Temperature > 5000) {
    Temperature = 0.94 * Temperature;
  } else if (Temperature < 5000) {
    Temperature = 0.8 * Temperature;
  } else if (Temperature < 1) {
    Temperature = 0.5 * Temperature;
  }
}

bool checkMove(long int prevCost) {
  // either accept or reject the move based on temperature & cost
  srand(time(NULL));
  double newCost = cost();
  int delCost = 0;
  double prob = drand48();
  if(debug) {
    cout << "new cost: " << newCost << endl;
  }
  delCost = newCost - prevCost;

  if (delCost <= 0 || prob > (exp(-delCost/Temperature))) {
    prevCost = newCost;
    accept_history.push_back(1);
    if (accept_history.size() > 100) {
      accept_ratio = ((accept_ratio*100) - (accept_history[accept_history.size()-100]) + 1)/100;
    } else {
      accept_ratio = accumulate(accept_history.begin(), accept_history.end(),0) / accept_history.size();
    }
    return true;
  } else {
    accept_history.push_back(0);
    if (accept_history.size() > 100) {
      accept_ratio = ((accept_ratio*100) - (accept_history[accept_history.size()-100]) + 0)/100;
    } else {
      accept_ratio = accumulate(accept_history.begin(), accept_history.end(),0) / accept_history.size();
    }
    return false;
  }
}

float timberWolfAlgorithm() {
  // main driver of the anealer
  //xLimit = macroPlacement();
  initialPlacement();

  Temperature = 500000;
  int i; // n * 20, n is number of nodes
  int cnt = 0;
  int ii = 0;

  high_resolution_clock::time_point t1 = high_resolution_clock::now();

  //while (Temperature > 0.1) {
  while (ii < 20000) {
    if(ii % 100 == 0) {
      high_resolution_clock::time_point t2 = high_resolution_clock::now();
      duration<double> time_span = duration_cast< duration<double> >(t2 - t1);

      cout << "******" << ii << "******" << endl;
      cout << "time: " <<  time_span.count() << " (s)" << endl;
      cout << "temperature: " << Temperature << endl;
      cout << "acceptance ratio: " << accept_ratio << endl;
      cost(-1);
    }

    i = 5 * nodeId.size();
    while (i > 0) {
      initiateMove();
      i--;
      cnt += 1;
    }
    update_temperature();
    ii += 1;
  }
  return cost();
}

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
