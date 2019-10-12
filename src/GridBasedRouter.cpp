//GridBasedRouter.cpp
#include "GridBasedRouter.h"

bool GridBasedRouter::outputResults2KiCadFile(std::vector<MultipinRoute> &nets)
{
  std::ifstream ifs;
  ifs.open(mDb.getFileName(), std::ifstream::in);
  if (!ifs)
  {
    std::cerr << "Cannot open input file: " << mDb.getFileName() << std::endl;
    return false;
  }

  std::vector<std::string> fileLines;
  while (!ifs.eof())
  {
    std::string inputLine;
    std::getline(ifs, inputLine);
    fileLines.push_back(inputLine);
    //std::cout << inputLine << std::endl;
  }

  ifs.close();

  if (fileLines.empty())
  {
    std::cerr << "Input file has no content: " << mDb.getFileName() << std::endl;
    return false;
  }

  // Handle output filename
  std::string fileExtension = util::getFileExtension(mDb.getFileName());
  std::string fileNameWoExtension = util::getFileNameWoExtension(mDb.getFileName());
  std::string outputFileName = fileNameWoExtension + ".routed.ours." + fileExtension;
  outputFileName = util::appendDirectory(GlobalParam::gOutputFolder, outputFileName);
  std::cout << __FUNCTION__ << "() outputFileName: " << outputFileName << std::endl;

  std::ofstream ofs;
  ofs.open(outputFileName, std::ofstream::out);
  if (!ofs)
  {
    std::cerr << "Cannot open output file: " << outputFileName << std::endl;
    return false;
  }
  // Remove the latest ")"
  for (int i = (int)fileLines.size() - 1; i >= 0; --i)
  {
    size_t found = fileLines.at(i).find_last_of(")");
    if (found != string::npos)
    {
      std::string newLine = fileLines.at(i).substr(0, found);
      fileLines.at(i) = newLine;
      break;
    }
  }
  // Paste inputs to outputs
  for (auto str : fileLines)
  {
    ofs << str << std::endl;
  }

  if (!writeNets(nets, ofs))
  {
    std::cerr << "Failed to write nets to the output file: " << outputFileName << std::endl;
    ofs.close();
    std::remove(outputFileName.c_str());
    return false;
  }

  ofs << ")" << std::endl;
  ofs.close();
  return true;
}

bool GridBasedRouter::writeNets(std::vector<MultipinRoute> &multipinNets, std::ofstream &ofs)
{
  if (!ofs)
    return false;

  // Set output precision
  ofs << std::fixed << std::setprecision(GlobalParam::gOutputPrecision);
  // Estimated total routed wirelength
  double totalEstWL = 0.0;

  // Multipin net
  for (auto &mpNet : multipinNets)
  {
    if (!mDb.isNetId(mpNet.netId))
    {
      std::cerr << __FUNCTION__ << "() Invalid net id: " << mpNet.netId << std::endl;
      continue;
    }

    auto &net = mDb.getNet(mpNet.netId);
    if (!mDb.isNetclassId(net.getNetclassId()))
    {
      std::cerr << __FUNCTION__ << "() Invalid netclass id: " << net.getNetclassId() << std::endl;
      continue;
    }

    auto &netclass = mDb.getNetclass(net.getNetclassId());
    Location last_location = mpNet.features[0];

    for (int i = 1; i < mpNet.features.size(); ++i)
    {
      auto &feature = mpNet.features[i];
      // check if close ???????
      if (
          abs(feature.x - last_location.x) <= 1 &&
          abs(feature.y - last_location.y) <= 1 &&
          abs(feature.z - last_location.z) <= 1)
      {
        // Print Through Hole Via
        if (feature.z != last_location.z)
        {
          ofs << "(via";
          ofs << " (at " << grid_factor * (last_location.x + mMinX * inputScale - enlargeBoundary / 2) << " " << grid_factor * (last_location.y + mMinY * inputScale - enlargeBoundary / 2) << ")";
          ofs << " (size " << netclass.getViaDia() << ")";
          ofs << " (drill " << netclass.getViaDrill() << ")";
          ofs << " (layers Top Bottom)";
          ofs << " (net " << mpNet.netId << ")";
          ofs << ")" << std::endl;
        }

        // Print Segment/Track/Wire
        if (feature.x != last_location.x || feature.y != last_location.y)
        {
          point_2d start{grid_factor * (last_location.x + mMinX * inputScale - enlargeBoundary / 2), grid_factor * (last_location.y + mMinY * inputScale - enlargeBoundary / 2)};
          point_2d end{grid_factor * (feature.x + mMinX * inputScale - enlargeBoundary / 2), grid_factor * (feature.y + mMinY * inputScale - enlargeBoundary / 2)};
          totalEstWL += point_2d::getDistance(start, end);

          ofs << "(segment";
          ofs << " (start " << start.m_x << " " << start.m_y << ")";
          ofs << " (end " << end.m_x << " " << end.m_y << ")";
          ofs << " (width " << netclass.getTraceWidth() << ")";
          ofs << " (layer " << mGridLayerToName.at(feature.z) << ")";
          ofs << " (net " << mpNet.netId << ")";
          ofs << ")" << std::endl;
        }
      }
      last_location = feature;
    }
  }

  std::cout << "=================" << __FUNCTION__ << "=================" << std::endl;
  std::cout << "\tEstimated Total WL: " << totalEstWL << std::endl;
  return true;
}

void GridBasedRouter::addPinCost(const pin &p, const float cost)
{
  // TODO: Id Range Checking?
  auto &comp = mDb.getComponent(p.m_comp_id);
  auto &inst = mDb.getInstance(p.m_inst_id);
  auto &pad = comp.getPadstack(p.m_padstack_id);

  addPinCost(pad, inst, cost);
}

void GridBasedRouter::addPinCost(const padstack &pad, const instance &inst, const float cost)
{
  point_2d pinDbLocation;
  mDb.getPinPosition(pad, inst, &pinDbLocation);
  double width = 0, height = 0;
  mDb.getPadstackRotatedWidthAndHeight(inst, pad, width, height);
  point_2d pinDbUR{pinDbLocation.m_x + width / 2.0, pinDbLocation.m_y + height / 2.0};
  point_2d pinDbLL{pinDbLocation.m_x - width / 2.0, pinDbLocation.m_y - height / 2.0};
  point_2d pinGridLL, pinGridUR;
  dbPointToGridPoint(pinDbUR, pinGridUR);
  dbPointToGridPoint(pinDbLL, pinGridLL);

  // TODO: Unify Rectangle to set costs
  const auto &layers = pad.getLayers();
  for (auto &layer : layers)
  {
    const auto &layerIte = mLayerNameToGrid.find(layer);
    if (layerIte != mLayerNameToGrid.end())
    {
      for (int x = pinGridLL.m_x; x < (int)pinGridUR.m_x; ++x)
      {
        for (int y = pinGridLL.m_y; y < (int)pinGridUR.m_y; ++y)
        {
          Location gridPt{x, y, layerIte->second};
          if (!mBg.validate_location(gridPt))
          {
            continue;
          }
          mBg.base_cost_add(cost, gridPt);
        }
      }
    }
  }
}

void GridBasedRouter::add_pin_cost_to_via_cost(const pin &p, const float cost)
{
  // TODO: Id Range Checking?
  auto &comp = mDb.getComponent(p.m_comp_id);
  auto &inst = mDb.getInstance(p.m_inst_id);
  auto &pad = comp.getPadstack(p.m_padstack_id);

  add_pin_cost_to_via_cost(pad, inst, cost);
}

void GridBasedRouter::add_pin_cost_to_via_cost(const padstack &pad, const instance &inst, const float cost)
{
  point_2d pinDbLocation;
  mDb.getPinPosition(pad, inst, &pinDbLocation);
  double width = 0, height = 0;
  mDb.getPadstackRotatedWidthAndHeight(inst, pad, width, height);
  point_2d pinDbUR{pinDbLocation.m_x + width / 2.0, pinDbLocation.m_y + height / 2.0};
  point_2d pinDbLL{pinDbLocation.m_x - width / 2.0, pinDbLocation.m_y - height / 2.0};
  point_2d pinGridLL, pinGridUR;
  dbPointToGridPoint(pinDbUR, pinGridUR);
  dbPointToGridPoint(pinDbLL, pinGridLL);
  //std::cout << __FUNCTION__ << "() cost:" << cost << ", inst:" << inst.getName() << "(" << inst.getId() << "), pad:" << pad.getName() << ", at(" << pinDbLocation.m_x << ", " << pinDbLocation.m_y << "), w:" << width << ", h:" << height << std::endl;

  // TODO: Unify Rectangle to set costs
  const auto &layers = pad.getLayers();
  for (auto &layer : layers)
  {
    const auto &layerIte = mLayerNameToGrid.find(layer);
    if (layerIte != mLayerNameToGrid.end())
    {
      for (int x = pinGridLL.m_x; x < (int)pinGridUR.m_x; ++x)
      {
        for (int y = pinGridLL.m_y; y < (int)pinGridUR.m_y; ++y)
        {
          Location gridPt{x, y, layerIte->second};
          if (!mBg.validate_location(gridPt))
          {
            //std::cout << "\tWarning: Out of bound, pin cost at " << gridPt << std::endl;
            continue;
          }
          //std::cout << "\tAdd pin cost at " << gridPt << std::endl;
          mBg.via_cost_add(cost, gridPt);
        }
      }
    }
  }
} 

bool GridBasedRouter::dbPointToGridPoint(const point_2d &dbPt, point_2d &gridPt)
{
  //TODO: boundary checking
  //TODO: consider integer ceiling or flooring???
  gridPt.m_x = dbPt.m_x * inputScale - mMinX * inputScale + enlargeBoundary / 2;
  gridPt.m_y = dbPt.m_y * inputScale - mMinY * inputScale + enlargeBoundary / 2;
  return true;
}

bool GridBasedRouter::gridPointToDbPoint(const point_2d &gridPt, point_2d &dbPt)
{
  //TODO: boundary checking
  //TODO: consider integer ceiling or flooring???
  dbPt.m_x = grid_factor * (gridPt.m_x + mMinX * inputScale - enlargeBoundary / 2);
  dbPt.m_y = grid_factor * (gridPt.m_y + mMinY * inputScale - enlargeBoundary / 2);
  return true;
}

void GridBasedRouter::test_router()
{
  std::vector<std::set<std::pair<double, double>>> routerInfo;
  mDb.getPcbRouterInfo(&routerInfo);

  std::cout << std::fixed << std::setprecision(5);
  std::cout << "=================test_placer==================" << std::endl;

  for (int i = 0; i < routerInfo.size(); ++i)
  {
    std::cout << "netId: " << i << std::endl;

    for (auto ite : routerInfo.at(i))
    {
      std::cout << " (" << ite.first << ", " << ite.second << ")";
      mMinX = std::min(ite.first, mMinX);
      mMaxX = std::max(ite.first, mMaxX);
      mMinY = std::min(ite.second, mMinY);
      mMaxY = std::max(ite.second, mMaxY);
    }
    std::cout << std::endl;
  }

  std::cout << "Routing Outline: (" << mMinX << ", " << mMinY << "), (" << mMaxX << ", " << mMaxY << ")" << std::endl;

  mDb.getBoardBoundaryByPinLocation(mMinX, mMaxX, mMinY, mMaxY);

  std::cout << "Routing Outline From DB: (" << mMinX << ", " << mMinY << "), (" << mMaxX << ", " << mMaxY << ")" << std::endl;

  // Initialize board grid
  const unsigned int h = int(std::abs(mMaxY * inputScale - mMinY * inputScale)) + enlargeBoundary;
  const unsigned int w = int(std::abs(mMaxX * inputScale - mMinX * inputScale)) + enlargeBoundary;
  const unsigned int l = 2;

  std::vector<MultipinRoute> multipinNets;

  std::cout << "BoardGrid Size: w:" << w << ", h:" << h << ", l:" << l << std::endl;

  mBg.initilization(w, h, l);
  mBg.base_cost_fill(0.0);

  // Prepare all the nets to route
  for (int i = 0; i < routerInfo.size(); ++i)
  {
    std::cout << "Routing netId: " << i << "..." << std::endl;
    size_t netDegree = routerInfo.at(i).size();
    if (netDegree < 2)
      continue;

    std::vector<Location> pins;
    for (auto ite : routerInfo.at(i))
    {
      std::cout << "\t(" << ite.first << ", " << ite.second << ")";
      pins.push_back(Location(ite.first * inputScale - mMinX * inputScale + enlargeBoundary / 2, ite.second * inputScale - mMinY * inputScale + enlargeBoundary / 2, 0));
      std::cout << " location in grid: " << pins.back() << std::endl;
    }

    //continue;
    //MultipinRoute mpRoute(pins);
    multipinNets.push_back(MultipinRoute(pins, i));
    mBg.add_route(multipinNets.back());
    std::cout << "Routing netId: " << i << "...done!" << std::endl;
  }

  // Routing has done
  // Print the final base cost
  //mBg.printGnuPlot();
  //mBg.printMatPlot();

  //outputResults2KiCadFile(multipinNets);
}
