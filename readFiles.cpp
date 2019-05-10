#include <iostream>
#include <vector>
#include <fstream>
#include <boost/regex.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <string>
#include <iostream>
#include <regex>
#include <map>
#include "readFiles.h"
#include "readScl.h"

using namespace std;

void readNodesFile() {
  fstream file;
  string buf;
  int i = 0;
  vector < string > strVec;
  using boost::is_any_of;
  int value = 2;

  file.open("apte.nodes", ios:: in );

  while (getline(file, buf)) {
    i++;
    if (i > 7) {
      boost::algorithm::split(strVec, buf, is_any_of("\t, "), boost::token_compress_on);
      Node n;
      if (strVec[3] == "terminal" || strVec[3] == "terminal_NI") {
        value = 0;
      } else {
        value = 1;
      }
      n.setParameterNodes(strVec[0],atof(strVec[1].c_str()), atof(strVec[2].c_str()), value);
      nodeId.insert(pair < string, Node > (strVec[0], n));
    }
  }
  file.close();
}

void readWtsFile() {
  fstream file;
  string buf;
  int i = 0;
  vector < string > strVec;
  using boost::is_any_of;
  map < string, Node > ::iterator itr;

  file.open("ibm01.wts", ios:: in );

  while (getline(file, buf)) {
    i++;
    if (i > 5) {
      boost::algorithm::split(strVec, buf, is_any_of("\t,  "), boost::token_compress_on);
      nodeId[strVec[1]].setParameterWts(atof(strVec[2].c_str()));
    }
  }
  file.close();
}

void readPlFile() {
  fstream file;
  string buf;
  int i = 0;
  int value = 0;
  vector < string > strVec;
  using boost::is_any_of;

  file.open("apte.pl", ios:: in );

  while (getline(file, buf)) {
    i++;
    if (i > 4) {
      boost::algorithm::split(strVec, buf, is_any_of("\t,  "), boost::token_compress_on);
      if (strVec[0] == ""){
        continue;
      }
      if (strVec[4] == "/FIXED" || strVec[4] == "/FIXED_NI") {
        value = 1;
      } else {
        value = 0;
      }
      nodeId[strVec[0]].setParameterPl(atof(strVec[1].c_str()), atof(strVec[2].c_str()), strVec[4], value);
    }
  }
  file.close();
}

void readNetsFile() {
  fstream file;
  string buf;
  int i = 0, a = 0, j = 0, NetId = 1;
  vector < string > strVec;

  regex find("\\b(NetDegree : )");
  smatch match;
  string Out;

  file.open("apte.nets", ios:: in );

  while (getline(file, buf)) {
    i++;
    if (i > 7) {
      using boost::is_any_of;
      regex_search(buf, match, find);
      Out = match.suffix();
      a = atof(Out.c_str());
      vector < string > strTemp;
      for (j = 0; j < a; j++) {

        getline(file, buf);
        boost::algorithm::split(strVec, buf, is_any_of("\t,  "), boost::token_compress_on); // a1213	 I : 0.5 0.5
        strTemp.push_back(strVec[0]);
        nodeId[strVec[0]].setNetList(NetId);
      }
      netToCell.insert(pair < int, vector < string > > (NetId, strTemp));
      NetId++;
    }
  }
}

void writePlFile(int idx) {
  fstream file;
  string buf;
  int i = 0;
  vector < string > strVec;
  using boost::is_any_of;

  file.open("out.pl", ios:: in );
  file.close();
}

void printMap() {
  map < int, vector < string > > ::iterator itr;
  vector < string > ::iterator itr1;

  cout << "netToCellMap" << endl;
  for (itr = netToCell.begin(); itr != netToCell.end(); ++itr) {
    cout << itr -> first << "	";

    for (itr1 = itr -> second.begin(); itr1 != itr -> second.end(); ++itr1) {
      cout << * itr1 << "	";
    }

    cout << endl;
  }
  cout << "\n" << endl;
}
