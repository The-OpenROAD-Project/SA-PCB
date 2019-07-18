# RePlAce
*SA-PCB: Simulated Annealing-based Placement For PCB Layout*

## Getting Started

## Install on a bare-metal machine

### Pre-requisite
* GCC compiler and libstdc++ static library >= 4.8.5
* boost library >= 1.62
* Python >= 3.6
* shapely
* matplotlib
* docopt
* numpy
* Recommended OS: Centos6, Centos7 or Ubuntu 16.04

### Clone repo and submodules
    $ git clone --recursive https://github.com/choltz95/c-pcb-annealer
    $ cd ~/c-pcbannealer
    $ make install
    $ make

### Fix eagle2bookshelf errors
    You may have to fix some eagle2bookshelf errors. Download DRU.py from https://github.com/NVSL/Swoop/blob/master/Swoop/DRU.py
    and place into the eagle2bookshelf directory. Test the installation by running their testcases.

### Check your installation [in development]
    To make sure your installation is correct and the current tool version is stable enough,
    run a Hello World application:

    $ make test

### How to execute
    ./sa parameters
       -i <optional, value> : for denoting # outer iterations PER SA INSTANCE
       -j <optional, value> : for denoting 'j'*#nodes inner iterations
       -t <optional, value> : for denoting initial temperature
       -f <optional, str>   : for output filename
       -e <optional, float> : convergence epsilon
       -v <optional>        : ben-amur flag
       -x <optional, int>   : simulated annealing instance index
       -p <required, string>: input placement board
       -d <optional, {0-3}> : debug verbosity
       -r <optional, {0-3}> : rotation
       EXAMPLE: ./sa -i 20000 -j 20 -t 1 -p input -f output

### License
* BSD-3-clause License [[Link]](LICENSE)

## Description
 C++ annealer for simple pcb placement of polygonal components.
 Supports the following:
  - Analytical overlap for arbitrary polygons using boost geometries
  - HPWL cost term for wirelength
  - BEN-AMEUR et al. cost-normalization & automatic initial temperature
  - Variable placement shift window, smaller displacement with temperature
  - Geometric temperature updates (cooling schedule) based on timberwolf schedule
  - Bookshelf parser
  - Bookshelf version translator
  - Plotting & animations

## Currently has issues / in development
  - Parallel multistart [New issues with cpp-taskflow. Still support for multi-start via python script.]
  - rtree spatial indexing for fast overlap computation

### Authors
- Chester Holtz, Devon Merrill, James (Ting-Chou) Lin, Connie (Yen-Yi) Wu (respective Ph.D. advisors: Chung-Kuan Cheng, Andrew B. Kahng).

### Limitations
 - Algorithm very sensitive to parameters
 - Broken support for weighted modules/nets
