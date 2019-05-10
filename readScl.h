#include <iostream>
#include <vector>
#include <fstream>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <string>
#include <iostream>
#include <map>

using namespace std;

void readSclFile();

class row;
extern map < int, row > rowId;
extern int RowWidth;
extern int xLimit;
class row {
  public:
    int Id;
  int coOrdinate;
  int height;
  int siteWidth;
  int siteSpacing;
  string siteOrient;
  string siteSymmetry;
  int siteRowOrigin;
  int numSites;
  int overlap;
  vector < string > cellList;

  void setId(int Id) {
    this -> Id = Id;
  }

  void setParameterRows(int coOrdinate, int height, int siteWidth, int siteSpacing, string siteOrient, string siteSymmetry, int siteRowOrigin, int numSites) {
    this -> coOrdinate = coOrdinate;
    this -> height = height;
    this -> siteWidth = siteWidth;
    this -> siteSpacing = siteSpacing;
    this -> siteOrient = siteOrient;
    this -> siteSymmetry = siteSymmetry;
    this -> siteRowOrigin = siteRowOrigin;
    this -> numSites = numSites;
  }

  void setCellList(string cellId) {
    cellList.push_back(cellId);
  }

  vector < string > sortByX() {
    int x = 0;
    map < int, string > sortx;
    map < int, string > ::iterator it;
    vector < string > ::iterator itl;
    vector < string > list;
    for (unsigned long i = 0; i < this -> cellList.size(); i++) {
      x = nodeId.find(cellList[i]) -> second.xCoordinate;
      sortx.insert(pair < int, string > (x, cellList[i]));
    }
    for (it = sortx.begin(); it != sortx.end(); ++it) {
      list.push_back(it -> second);
    }
    this -> cellList = list;
    return this -> cellList;
  }

  void calcRowOverlap() {
    vector < string > ::iterator it1;
    int xLast = 0, widthLast = 0;

    xLast = nodeId[cellList[cellList.size() - 1]].xCoordinate;
    widthLast = nodeId[cellList[cellList.size() - 1]].width;
    overlap = xLast + widthLast - (RowWidth + xLimit);

    this -> overlap = overlap;
  }

  void printParameter() {
    cout << "rowId          " << this -> Id << " " << endl;
    cout << "Row-Co-ordinate " << this -> coOrdinate << endl;
    cout << "height          " << this -> height << endl;
    cout << "siteWidth       " << this -> siteWidth << endl;
    cout << "siteSpacing     " << this -> siteSpacing << endl;
    cout << "siteOrientation " << this -> siteOrient << endl;
    cout << "siteSymmetry    " << this -> siteSymmetry << endl;
    cout << "siteRowOrigin   " << this -> siteRowOrigin << endl;
    cout << "numSites        " << this -> numSites << endl;
    cout << "Overlap of this row   " << overlap << endl;
    cout << "cellsInRow      " << " ";
    vector < string > ::iterator it1;
    for (it1 = cellList.begin(); it1 != cellList.end(); ++it1) {
      cout << * it1 << " ";
    }
    cout << "\n" << endl;
  }
};
