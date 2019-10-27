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


#include "kicadPcbDataBase.h"
#include "GridBasedPlacer.hpp"
#include "util.h"
int test_flow()
{
  GlobalParam::setFolders();
  GlobalParam::setUsageStart();

  std::string designName = "bm1";
  std::cout << "Parsing design: " << designName << std::endl;
  kicadPcbDataBase db(designName);

  db.printLayer();
  db.printComp();
  db.printInst();
  db.printNetclass();
  db.printNet();
  db.printFile();
  db.printPcbRouterInfo();

  GlobalParam::showCurrentUsage("Parser");
  GlobalParam::setUsageStart();

  std::cout << "Starting placer..." << std::endl;
  GridBasedPlacer placer(db);
  placer.test_placer_flow();

  GlobalParam::showCurrentUsage("GridBasedPlacer");
  GlobalParam::showFinalUsage("End of Program");

  return 0;
}
int main(int argc, char *argv[])
{
  if (argc < 2)
  {
    std::cout << "Please provide input testcase filename." << std::endl;
    return 0;
  }

  //util::showSysInfoComdLine(argc, argv);
  GlobalParam::setFolders();
  GlobalParam::setUsageStart();

  std::string designName = argv[1];
  std::cout << "Parsing design: " << designName << std::endl;
  kicadPcbDataBase db(designName);

  db.printLayer();
  db.printComp();
  db.printInst();
  db.printNetclass();
  db.printNet();
  db.printFile();
  db.printPcbRouterInfo();

  GlobalParam::showCurrentUsage("Parser");
  GlobalParam::setUsageStart();

  std::cout << "Starting placer..." << std::endl;
  GridBasedPlacer placer(db);
  placer.test_placer_flow();

  GlobalParam::showCurrentUsage("GridBasedPlacer");
  GlobalParam::showFinalUsage("End of Program");

  return 0;
}
