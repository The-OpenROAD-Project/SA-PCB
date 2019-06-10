#include <iostream>
#include <vector>
#include <fstream>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <string>
#include <iostream>
#include <map>
#include "readFiles.h"
//#include "readScl.h"

using namespace std;

void readSclFile() {
  fstream file;
  string buf;
  int i = 0, j = 0, Id = 1;
  vector < string > strVec;
  using boost::is_any_of;

  int coOrdinate;
  int height;
  int siteWidth;
  int siteSpacing;
  string siteOrient;
  string siteSymmetry;
  int siteRowOrigin;
  int numSites;

  file.open("ibm01.scl", ios:: in );

  while (getline(file, buf)) {
    i++;
    if (i > 8) {
      boost::algorithm::split(strVec, buf, is_any_of("\t,  "), boost::token_compress_on);
      j = i % 9;
      if (j == 1) {
        coOrdinate = atoi(strVec[3].c_str());
      } else if (j == 2) {
        height = atoi(strVec[3].c_str());
      } else if (j == 3) {
        siteWidth = atoi(strVec[3].c_str());
      } else if (j == 4) {
        siteSpacing = atoi(strVec[3].c_str());
      } else if (j == 5) {
        siteOrient = strVec[3];
      } else if (j == 6) {
        siteSymmetry = strVec[3];
      } else if (j == 7) {
        siteRowOrigin = atoi(strVec[3].c_str());
        numSites = atoi(strVec[6].c_str());
      } else if (j == 8) {
        row r;
        r.setId(Id);
        rowId.insert(pair < int, row > (Id, r));
        rowId[Id].setParameterRows(coOrdinate, height, siteWidth, siteSpacing, siteOrient, siteSymmetry, siteRowOrigin, numSites);
        Id++;
      }

    }
  }
  file.close();
}
