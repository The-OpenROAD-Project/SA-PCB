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

#include "readFiles.hpp"

using namespace std;
using boost::is_any_of;

int readNodesFile(string fname) {
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
      n.setParameterNodes(strVec[0], atof(strVec[1].c_str()), atof(strVec[2].c_str()), value, idx); // /1000000
      nodeId.push_back(n);
      name2id.insert(pair < string, int > (strVec[0], idx));
      idx += 1;
    }
  }
  file.close();
  return 0;
}

int readShapesFile(string fname) {
  fstream file;
  string buf;
  int i = 0;
  vector < string > strVec;

  file.open(fname, ios:: in );
  while (getline(file, buf)) {
    i++;
    if (i > 7) {
      boost::trim(buf);
      boost::algorithm::split(strVec, buf, is_any_of("\t,  "), boost::token_compress_on);
      if (strVec[0] == "" || strVec[0] == " "){
        continue;
      }
      strVec.erase(strVec.begin());
      nodeId[name2id[strVec[0]]].setParameterShapes(boost::algorithm::join(strVec,""));
    }
  }
  file.close();
  return 0;
}

int readWtsFile(string fname) {
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
      nodeId[name2id[strVec[1]]].setParameterWts(atof(strVec[2].c_str()));
    }
  }
  file.close();
  return 0;
}

int readPlFile(string fname) {
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
      nodeId[name2id[strVec[0]]].setParameterPl(atof(strVec[1].c_str()), atof(strVec[2].c_str()), strVec[4], value);
    }
  }
  file.close();
  return 0;
}

map<int, vector<Pin> > readNetsFile(string fname) {
  fstream file;
  string buf;
  int i = 0, a = 0, j = 0, NetId = 1;
  vector < string > strVec;
  map<int, vector<Pin> > netToCell;

  string pat = "NetDegree : ";
  string Out;

  file.open(fname, ios:: in );
  while (getline(file, buf)) {
    i++;
    if (i > 7) {
      boost::trim_all(buf);

      std::size_t found = buf.find(pat);
      if(found!=std::string::npos) {
        std::vector<std::string> results;
        boost::split(results, buf, [](char c){return c == ' ';});
        //Out = buf.substr(buf.rfind(" ") + 1);
        Out = results.at(2);
      } else {
        continue;
      }
      a = stoi(Out);
      vector < string > strTemp;
      vector < Pin > pinTemp;
      for (j = 0; j < a; j++) {
        getline(file, buf);
        boost::trim_all(buf);
        boost::algorithm::split(strVec, buf, is_any_of("\t,  "), boost::token_compress_on);
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

int writePlFile(string fname) {
  vector < string > strVec;
  fstream file;
  string buf;
  ofstream myfile (fname);
  if (myfile.is_open()) {
    myfile << "\n\n\n\n";
    vector <Node> ::iterator itNode;

    // print components
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
    // print terminals
    for (itNode = nodeId.begin(); itNode != nodeId.end(); ++itNode) {
      if(itNode->terminal) {
        myfile << itNode->name << " " << itNode->xCoordinate << " " << itNode->yCoordinate << " : " << itNode->orient2str(itNode->orientation);
        myfile << " /FIXED_NI\n";
      }
    }
    myfile.close();
    return 0;
  } else{ cout << "[ERR] Unable to open cache dir" <<endl; return 1;}
}
