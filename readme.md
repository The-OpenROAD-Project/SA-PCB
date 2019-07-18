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
       -v <optional>        : Ben-amur flag
       -x <optional, int>   : simulated annealing instance index
       -p <required, string>: input placement board
       -d <optional, {0-3}> : debug verbosity
       -r <optional, {0-3}> : rotation
       EXAMPLE: ./sa -i 20000 -j 20 -t 1 -p input -f output

### License
* BSD-3-clause License [[Link]](LICENSE)

## Description
 C++ annealer for simple PCB placement of polygonal components.
 Supports the following:
  - Analytical overlap for arbitrary polygons using boost geometries
  - Support 90 deg., 45 deg., free rotation [free rotation in progress]
  - HPWL cost term for wirelength
  - BEN-AMEUR et al. cost-normalization & automatic initial temperature
  - Variable placement shift window, smaller displacement with temperature
  - Geometric temperature updates (cooling schedule) based on Timberwolf schedule
  - Bookshelf parser
  - Bookshelf version translator
  - Plotting & animations

## Currently has issues / in development
  - Parallel multistart [New issues with cpp-taskflow. Still support for multi-start via python script.]
  - R-Tree spatial indexing for fast overlap computation

### Authors
  - Chester Holtz, Devon Merrill, James (Ting-Chou) Lin, Connie (Yen-Yi) Wu (Ph.D. advisor: Chung-Kuan Cheng, Steven Swanson).

### Limitations
  - Algorithm very sensitive to parameters
  - Broken support for weighted modules/nets

### Citations
  - Sechen, C. and Sangiovanni-Vincentelli, A. L., "The Timber-Wolf placement and routing package", IEEE J. Solid-State Circuits 30:510–522 1985.
  - C.-K. Cheng, A. B. Kahng, I. Kang and L. Wang, "RePlAce: Advancing Solution Quality and Routability Validation in Global Placement", to appear in IEEE Transactions on Computer-Aided Design of Integrated Circuits and Systems, 2018. (Digital Object Identifier: 10.1109/TCAD.2018.2859220)
  - Ben-Ameur, Walid. "Computing the initial temperature of simulated annealing." Computational Optimization and Applications 29, no. 3 (2004): 369-385.
  - James M. Varanelli and James P. Cohoon. Two-stage simulated annealing methodology. In Proceedings of the 5th Great Lakes Symposium on VLSI, volume ?, pages 50–53, Buffalo, NY, 16.-18. March 1995. IEEE, Los Alamitos, CA. †EI M153001/95
