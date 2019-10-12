#include "plotter.h"

Plotter::Plotter(char *outputFileName)
{
    //_objCounter = 1;
    _outFile.open(outputFileName, ofstream::out);
    _outFile << fixed;
    if (!_outFile.is_open())
    {
        cout << "**Error. Cannot open file (" << outputFileName << ") in Plotter::Plotter" << endl;
    }
}
Plotter::~Plotter()
{
    if (_outFile.is_open())
    {
        _outFile.close();
    }
}

bool Plotter::isInit()
{
    return _outFile.is_open();
}

void Plotter::outputString(const string &s)
{
    _outFile << s;
}

void Plotter::outputDouble(double d, size_t digit)
{
    _outFile << std::fixed;
    _outFile << setprecision(digit) << d;
}

//void Plotter::outputString(char* s){
//    outFile << s;
//}

void Plotter::plotLine(double x1, double y1, double x2, double y2)
{
    _outFile << x1 << ", " << y1 << "\n"
             << x2 << ", " << y2 << "\n\n";
}

void Plotter::plotRec(double &x1, double &y1, double &x2, double &y2)
{
    _outFile << x1 << ", " << y1 << "\n"
             << x2 << ", " << y1 << "\n"
             << x2 << ", " << y2 << "\n"
             << x1 << ", " << y2 << "\n"
             << x1 << ", " << y1 << "\n\n";
}

void Plotter::plotRec(double x1, double y1, double x2, double y2)
{
    _outFile << x1 << ", " << y1 << "\n"
             << x2 << ", " << y1 << "\n"
             << x2 << ", " << y2 << "\n"
             << x1 << ", " << y2 << "\n"
             << x1 << ", " << y1 << "\n\n";
}
