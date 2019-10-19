//util.h
#ifndef PCBROUTER_UTIL_H
#define PCBROUTER_UTIL_H

#include <cstdlib>  // atof
#include <cassert>  // assert
#include <cmath>    // fabs
#include <iostream> // cout, ostream
#include <fstream>  // ifstream
#include <vector>
#include <sstream> // stringstream
#include <string>
#include <sys/resource.h>
#include <sys/time.h>
#include <sys/stat.h>  // Unix-only. Create Folder, change to <filesystem> when C++17 is ready
#include <sys/types.h> // Unix-only. Create Folder, change to <filesystem> when C++17 is ready

using namespace std;

namespace util
{
// =====================================================
// show system info. -----------------------------------
// =====================================================
inline void showSysInfoComdLine(int argc, char *argv[])
{
    cout << "Command line: ";
    for (int i = 0; i < argc; ++i)
    {
        cout << argv[i] << " ";
    }
    cout << endl;
    int systemret = 0;
    cout << "=================== SYSTEM INFORMATION ==================" << endl;
    systemret *= system("echo 'User:       '`whoami`@`hostname`");
    systemret *= system("echo 'Date:       '`date`");
    systemret *= system("echo 'System:     '`uname -a`");
    cout << "=========================================================" << endl;
    assert(!systemret);
}

// =====================================================
// filename --------------------------------------------
// =====================================================
inline string getFileDirName(const string filePathName)
{
    string retStr = filePathName;
    string::size_type pos = retStr.find_last_of("/\\");
    if (pos != string::npos)
        retStr = retStr.substr(0, pos);
    return retStr;
}
inline string getFileName(const string filePathName)
{
    string retStr = filePathName;
    string::size_type pos = retStr.find_last_of("/\\");
    if (pos != string::npos)
        retStr = retStr.substr(pos + 1);
    return retStr;
}
inline string getFileNameWoExtension(const string filePathName)
{
    string retStr = filePathName;
    string::size_type pos = retStr.find_last_of("/\\");
    if (pos != string::npos)
        retStr = retStr.substr(pos + 1);
    pos = retStr.find_last_of(".");
    if (pos != string::npos)
        retStr = retStr.substr(0, pos);
    return retStr;
}
inline string getFileExtension(const string filePathName)
{
    string retStr = filePathName;
    string::size_type pos = retStr.rfind(".");
    if (pos != string::npos)
        retStr = retStr.substr(pos + 1);
    return retStr;
}

// =====================================================
// Directory Unix-only (TODO: std::filesystem in C++17)-
// =====================================================
inline bool createDirectory(const string dirName)
{
    // Creating a directory
    struct stat info;

    if (stat(dirName.c_str(), &info) != 0)
    {
        cerr << "Cannot access directory " << dirName << endl;
    }
    else if (info.st_mode & S_IFDIR)
    {
        cout << dirName << " dirctory already exists" << endl;
    }
    else
    {
        if (mkdir(dirName.c_str(), 0777) == -1)
        {
            //cerr << "Error creating directory " << dirName << " :  " << strerror(errno) << endl;
            return false;
        }
        cout << "Directory created: " << dirName << endl;
    }
    return true;
}
inline string appendDirectory(const string dirName, const string fileName)
{
    string retStr = dirName + "/" + fileName;
    return retStr;
}

// =====================================================
// CPU time & memory usage -----------------------------
// =====================================================
#define TIME_SCALE 1000000.0
#define MEMORY_SCALE 1024.0
class TimeUsage
{
public:
    TimeUsage()
    {
        start(FULL);
        start(PARTIAL);
    }
    struct TimeState
    {
        TimeState(long r = 0, long u = 0, long s = 0) : rTime(r), uTime(u), sTime(s) {}
        long rTime, uTime, sTime; //real, user, system
    };

    enum TimeType
    {
        FULL,
        PARTIAL
    };
    void start(TimeType type) { (type == FULL) ? checkUsage(tStart_) : checkUsage(pStart_); }

    void showUsage(const string comment, TimeType type)
    {
        TimeState curSt;
        checkUsage(curSt);
        TimeState dur = (type == FULL) ? diff(tStart_, curSt) : diff(pStart_, curSt);
        if (type == FULL)
        {
            cout << "---------- " << comment << " total time usage -----------" << endl;
            cout << " Final Real:" << dur.rTime << "s;";
        }
        else
        {
            cout << "---------- " << comment << " period time usage -----------" << endl;
            cout << " Real:" << dur.rTime << "s;";
        }
        cout << " User:" << dur.uTime << "s;";
        cout << " System:" << dur.sTime << "s." << endl
             << endl;
    }

private:
    TimeState diff(TimeState &start, TimeState &end)
    {
        return TimeState(end.rTime - start.rTime, end.uTime - start.uTime, end.sTime - start.sTime);
    }

    void checkUsage(TimeState &st) const
    {
        rusage tUsg;
        getrusage(RUSAGE_SELF, &tUsg);
        timeval tReal;
        gettimeofday(&tReal, NULL);
        st.uTime = tUsg.ru_utime.tv_sec + tUsg.ru_utime.tv_usec / TIME_SCALE;
        st.sTime = tUsg.ru_stime.tv_sec + tUsg.ru_stime.tv_usec / TIME_SCALE;
        st.rTime = tReal.tv_sec + tReal.tv_usec / TIME_SCALE;
    }
    TimeState tStart_, pStart_; //total, period
};
// memory
inline double getPeakMemoryUsage()
{
#if defined(linux)
    char buf[1000];
    ifstream ifs("/proc/self/stat");
    for (int i = 0; i != 23; ++i)
        ifs >> buf;
    return (1.0 / (MEMORY_SCALE * MEMORY_SCALE) * atof(buf)); // GB
#else
    return -1;
#endif
}

// =====================================================
// math-related auxiliary ------------------------------
// =====================================================
inline bool is_equal(double a, double b)
{
    return (a - b) < 0.0001;
}
inline double interpolate(const double &a, const double &b, const double &ratio)
{
    return a + ((b - a) * ratio);
}

// =====================================================
// string-related conversion----------------------------
// =====================================================
enum Orient
{
    OR_N,
    OR_W,
    OR_S,
    OR_E,
    OR_FN,
    OR_FW,
    OR_FS,
    OR_FE,
    OR_OTHER
};
inline string orientStr(Orient orient)
{
    const char *orientString[] = {"N", "W", "S", "E", "FN", "FW", "FS", "FE"};
    return orientString[orient];
}
enum PinDir
{
    UNKNOWN,
    INPUT,
    OUTPUT
};
inline string dir2Str(PinDir pinDir)
{
    const char *dirString[] = {"unknown", "IN", "OUT"};
    return dirString[pinDir];
}
template <class T>
inline string num2str(const T &p)
{
    string num;
    std::stringstream out;
    out << p;
    num = out.str();
    return num;
}
inline string bool2Str(bool isTrue)
{
    if (isTrue)
        return "TRUE";
    else
        return "FALSE";
}

// =====================================================
// rectangle -------------------------------------------
// =====================================================
class Rect
{
public:
    Rect(double left = 0, double bottom = 0, double right = 0, double top = 0)
        : _left(left), _bottom(bottom), _right(right), _top(top)
    {
    }

    /////////////////////////////////////////////
    // get
    /////////////////////////////////////////////
    double left() const { return _left; }
    double bottom() const { return _bottom; }
    double right() const { return _right; }
    double top() const { return _top; }
    double width() const { return _right - _left; }
    double height() const { return _top - _bottom; }
    double centerX() const { return (_left + _right) / 2; }
    double centerY() const { return (_bottom + _top) / 2; }
    Rect scale(double s) const { return Rect(_left * s, _bottom * s, _right * s, _top * s); }
    Rect shift(double x, double y) const { return Rect(_left + x, _bottom + y, _right + x, _top + y); }

    /////////////////////////////////////////////
    // set
    /////////////////////////////////////////////
    void setBounds(const double &left, const double &bottom, const double &right, const double &top)
    {
        _left = left;
        _bottom = bottom;
        _right = right;
        _top = top;
    }
    void setTop(const double &top)
    {
        _top = top;
    }
    void setBottom(const double &bottom)
    {
        _bottom = bottom;
    }
    void setRight(const double &right)
    {
        _right = right;
    }
    void setLeft(const double &left)
    {
        _left = left;
    }

    /////////////////////////////////////////////
    // overlap area of two rectangles
    /////////////////////////////////////////////
    static double overlapArea(const Rect &rect1, const Rect &rect2)
    {
        double overlapH = min(rect1.right(), rect2.right()) - max(rect1.left(), rect2.left());
        double overlapV = min(rect1.top(), rect2.top()) - max(rect1.bottom(), rect2.bottom());
        if (overlapH < 0)
            overlapH = 0;
        if (overlapV < 0)
            overlapV = 0;
        return overlapH * overlapV;
    }
    void boundPosition(double &x, double &y)
    {
        x = max(min(x, _right), _left);
        y = max(min(x, _top), _bottom);
    }

    friend ostream &operator<<(ostream &, const Rect &);
    void operator=(const Rect &rect)
    {
        setBounds(rect.left(), rect.bottom(), rect.right(), rect.top());
    }
    friend bool operator==(const Rect &, const Rect &);

private:
    double _left;
    double _bottom;
    double _right;
    double _top;
};
inline ostream &operator<<(ostream &out, const Rect &rect)
{
    out << "rect: (" << rect.left() << "," << rect.bottom() << ")-(" << rect.right() << "," << rect.top() << ")";
    return out;
}
inline bool operator==(const Rect &r1, const Rect &r2)
{
    return (r1.left() == r2.left()) && (r1.bottom() == r2.bottom()) &&
           (r1.right() == r2.right()) && (r1.top() == r2.top());
}

// =====================================================
// Point -----------------------------------------------
// =====================================================
class Point
{
public:
    Point(double x = 0, double y = 0) : _x(x), _y(y) {}
    double x() const { return _x; }
    double y() const { return _y; }

    void set_x_y(double x, double y)
    {
        _x = x;
        _y = y;
    }
    void shift(double x, double y)
    {
        _x += x;
        _y += y;
    }
    void scale(double sx, double sy)
    {
        _x *= sx;
        _y *= sy;
    }

    static double dist(const Point &p1, const Point &p2)
    {
        return fabs(p1.x() - p2.x()) + fabs(p1.y() - p2.y());
    }
    friend ostream &operator<<(ostream &, const Point &);
    void operator=(const Point &p)
    {
        set_x_y(p.x(), p.y());
    }
    Point operator+(const Point &p)
    {
        return Point(_x + p.x(), _y + p.y());
    }
    Point operator-(const Point &p)
    {
        return Point(_x - p.x(), _y - p.y());
    }

private:
    double _x, _y;
};
inline ostream &operator<<(ostream &out, const Point &p)
{
    out << "(" << p.x() << "," << p.y() << ")";
    return out;
}

template <class T>
class Array2D
{
public:
    Array2D(int xDim, int yDim, T defultValue) : _xDim(xDim), _yDim(yDim)
    {
        values.resize(_xDim * _yDim, defultValue);
    }
    T g(int xId, int yId)
    {
        return values[yId * _xDim + xId];
    }
    void set(int xId, int yId, T v)
    {
        values[yId * _xDim + xId] = v;
    }

private:
    vector<T> values;
    int _xDim, _yDim;
};
} // namespace util

#endif // UTIL_H
