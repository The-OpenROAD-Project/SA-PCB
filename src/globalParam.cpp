#include "globalParam.h"

int GlobalParam::gLayerNum = 3;
double GlobalParam::gEpsilon = 0.00000000000001;
bool GlobalParam::g90DegreeMode = true;

float GlobalParam::HPWLWeight = 0.5;
float GlobalParam::DensityWeight = 0.5;
float GlobalParam::RoutabilityWeight = 0.0;

// Outputfile
int GlobalParam::gOutputPrecision = 5;
string GlobalParam::gOutputFolder = "output";
// logfile
string GlobalParam::gLogFolder = "log";

int GlobalParam::gSeed = 1470295829; //time(NULL);
const double GlobalParam::gSqrt2 = sqrt(2);
const double GlobalParam::gTan22_5 = tan(22.5 * PI / 180.0);

util::TimeUsage GlobalParam::runTime = util::TimeUsage();

void GlobalParam::showParam()
{
    /*
    cout << "=============PARAM============"
         << "\nDesign Related:"
         << "\ngLayerNum: " << gLayerNum
         << "\ngWireWidth: " << gWireWidth
         << "\ngWireSpace: " << gWireSpace
         << "\ngWirePadSpace: " << gWirePadSpace
         << "\n\nDesign Related (derived by info. above):"
         << "\ngHalfWireWidth: " << gHalfWireWidth
         << "\ngHWidthAndSpace: " << gHWidthAndSpace
         << "\ngHWidthAndWPSpace: " << gHWidthAndWPSpace
         << "\ngWire2Wire: " << gWire2Wire
         << "\ngDisBetweenWire: " << gDisBetweenWire
         << "\n\nAlgorithm Related:"
         << "\ngEpsilon: " << gEpsilon
         << "\ngSeed: " << gSeed << endl;
    */
}

void GlobalParam::setFolders()
{
    if(!util::createDirectory(gOutputFolder)){
        gOutputFolder = "";
    }
    if(!util::createDirectory(gLogFolder)){
        gLogFolder = "";
    }
}

void GlobalParam::showCurrentUsage(const string comment)
{
    runTime.showUsage(comment, util::TimeUsage::PARTIAL);
}

void GlobalParam::showFinalUsage(const string comment)
{
    runTime.showUsage(comment, util::TimeUsage::FULL);
}

void GlobalParam::setUsageStart()
{
    runTime.start(util::TimeUsage::PARTIAL);
}
