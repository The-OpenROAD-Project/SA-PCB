#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <map>

#include <boost/regex.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim_all.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/tuple/tuple.hpp>

#include "readFiles.h"
//#include "readScl.h"

using namespace std;
using boost::is_any_of;

void readNodesFile(string fname) {
  fstream file;
  string buf;
  int i = 0;
  vector < string > strVec;
  int value = 2;
  int idx = 0;

  file.open(fname, ios:: in );
  while (getline(file, buf)) {
    i++;
    if (i > 7) {
      boost::trim(buf);
      boost::algorithm::split(strVec, buf, is_any_of("\t, "), boost::token_compress_on);
      if (strVec[0] == "" || strVec[0] == " "){
        continue;
      }
      Node n;

      if (strVec.size() > 3 && (strVec[3] == "terminal" || strVec[3] == "terminal_NI")) {
        value = 1;
      } else {
        value = 0;
      }
      n.setParameterNodes(strVec[0], atof(strVec[1].c_str()), atof(strVec[2].c_str()), value, idx);
      //nodeId.insert(pair < string, Node > (strVec[0], n));
      nodeId.push_back(n);
      name2id.insert(pair < string, int > (strVec[0], idx));
      idx += 1;
    }
  }
  file.close();
}

void readWtsFile(string fname) {
  fstream file;
  string buf;
  int i = 0;
  vector < string > strVec;
  map < string, Node > ::iterator itr;

  file.open(fname, ios:: in );
  while (getline(file, buf)) {
    i++;
    boost::trim(buf);
    if (i > 5) {
      boost::algorithm::split(strVec, buf, is_any_of("\t,  "), boost::token_compress_on);
      //nodeId[strVec[1]].setParameterWts(atof(strVec[2].c_str()));
      nodeId[name2id[strVec[1]]].setParameterWts(atof(strVec[2].c_str()));
    }
  }
  file.close();
}

void readPlFile(string fname) {
  fstream file;
  string buf;
  int i = 0;
  int value = 0;
  vector < string > strVec;

  file.open(fname, ios:: in );
  while (getline(file, buf)) {
    i++;
    if (i > 4) {
      boost::trim(buf);
      boost::algorithm::split(strVec, buf, is_any_of("\t,  "), boost::token_compress_on);
      if (strVec[0] == "" || strVec[0] == " "){
        continue;
      }
      if (strVec.size() > 5 && (strVec[5] == "/FIXED" || strVec[5] == "/FIXED_NI")) {
        value = 1;
      } else {
        value = 0;
      }
      //nodeId[strVec[0]].setParameterPl(atof(strVec[1].c_str()), atof(strVec[2].c_str()), strVec[4], value);
      nodeId[name2id[strVec[0]]].setParameterPl(atof(strVec[1].c_str()), atof(strVec[2].c_str()), strVec[4], value);
    }
  }
  file.close();
}

map<int, vector<Pin> > readNetsFile(string fname) {
  fstream file;
  string buf;
  int i = 0, a = 0, j = 0, NetId = 1;
  vector < string > strVec;
  map<int, vector<Pin> > netToCell;

  //boost::regex pattern("\\b(NetDegree : )");
  //boost::smatch match;
  string pat = "NetDegree : ";
  string Out;

  file.open(fname, ios:: in );
  while (getline(file, buf)) {
    i++;
    if (i > 7) {
      boost::trim_all(buf);
      //if(boost::regex_search(buf, match, pattern)) {
      //  Out = match.suffix();
      std::size_t found = buf.find(pat);
      if(found!=std::string::npos) {
        Out = buf.substr(buf.rfind(" ") + 1);
      } else {
        continue;
      }
      //a = atof(Out.c_str());
      a = stoi(Out);
      vector < string > strTemp;
      vector < Pin > pinTemp;
      for (j = 0; j < a; j++) {
        getline(file, buf);
        boost::trim_all(buf);
        boost::algorithm::split(strVec, buf, is_any_of("\t,  "), boost::token_compress_on); // a1213	 I : 0.5 0.5
        strTemp.push_back(strVec[0]);
        nodeId[name2id[strVec[0]]].setNetList(NetId);
        Pin p;
        if(strVec.size() > 2) {
          p.set_params(strVec[0].c_str(), atof(strVec[3].c_str()), atof(strVec[4].c_str()), nodeId[name2id[strVec[0]]].idx);
        } else {
          p.set_params(strVec[0].c_str(), 0, 0, -1);
        }
        pinTemp.push_back(p);
      }
      netToCell.insert(pair < int, vector< Pin > > (NetId, pinTemp));
      NetId++;
    }
  }
  return netToCell;
}

void writePlFile(string fname) {
  vector < string > strVec;
  fstream file;
  string buf;
  ofstream myfile (fname);
  if (myfile.is_open()) {
    myfile << "\n\n\n\n";
    vector <Node> ::iterator itNode;

    //components
    for (itNode = nodeId.begin(); itNode != nodeId.end(); ++itNode) {
      if(!itNode->terminal) {
        myfile << itNode->name << " " << itNode->xCoordinate << " " << itNode->yCoordinate <<  " : " << itNode->orient2str(itNode->orientation);
        if (itNode->fixed) {
          myfile << " /FIXED_NI\n";
        } else {
          myfile << "\n";
        }
      }
    }

      myfile << "\n";
      //terminals
      for (itNode = nodeId.begin(); itNode != nodeId.end(); ++itNode) {
        if(itNode->terminal) {
          myfile << itNode->name << " " << itNode->xCoordinate << " " << itNode->yCoordinate << " : " << itNode->orientation_str;
          if (itNode->fixed) {
            myfile << " /FIXED_NI\n";
          } else {
            myfile << "\n";
          }
        }
      }
    myfile.close();
  } else{ cout << "Unable to open file"; }
}
/*
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
}*/
