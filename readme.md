# SA-PCB
  *SA-PCB: Simulated Annealing-based Placement For PCB Layout*

## Getting Started

## Installation

### Pre-requisites
  * GCC compiler >= 4.8.5
  * boost library >= 1.62
  * cpp-taskflow >= 2.0
  * Python >= 3.6
  * shapely >= 1.6.4
  * matplotlib >= 3.0.2
  * docopt >= 0.6.2
  * numpy >= 1.15.4
  * argparse >= 1.1
  * tqdm >= 4.30.0 [python-multistart]
  * numba >= 0.42.1 [python-multistart]
  * multiprocessing [python-multistart]
  * yaml
  * json
  * Recommended OS: Centos6, Centos7 or Ubuntu 16.04

### Clone repo and submodules
    $ git clone --recursive https://github.com/choltz95/c-pcb-annealer
    $ cd ~/c-pcbannealer
    $ pip install -r requirements.txt
    $ make install
    $ make

### Fix eagle2bookshelf errors
    You may have to fix some eagle2bookshelf errors. Download DRU.py and place into the eagle2bookshelf
    directory.

     $ wget https://raw.githubusercontent.com/NVSL/Swoop/master/Swoop/DRU.py
     $ mv DRU.py eagle2bookshelf


### Check your installation
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

## Limitations / Current issues / In development
  - Parallel multistart [Issues with cpp-taskflow. Still support for multi-start via python script.]
  - R-Tree spatial indexing for fast overlap computation
  - Free rotation
  - Algorithm very sensitive to parameters
  - Broken support for weighted modules/nets
  - Set up Dockerfile
  - Support json configuration files for algorithm parameters

### Authors
  - Chester Holtz, Devon Merrill, James (Ting-Chou) Lin, Connie (Yen-Yi) Wu (Ph.D. advisor: Chung-Kuan Cheng, Steven Swanson).
  - Pull requests to improve the tool are very appreciated. 

### Citations
  - Sechen, C. and Sangiovanni-Vincentelli, A. L., "The Timber-Wolf placement and routing package", IEEE J. Solid-State Circuits 30:510–522 1985.
  - C.-K. Cheng, A. B. Kahng, I. Kang and L. Wang, "RePlAce: Advancing Solution Quality and Routability Validation in Global Placement", to appear in IEEE Transactions on Computer-Aided Design of Integrated Circuits and Systems, 2018. (Digital Object Identifier: 10.1109/TCAD.2018.2859220)
  - Ben-Ameur, Walid. "Computing the initial temperature of simulated annealing." Computational Optimization and Applications 29, no. 3 (2004): 369-385.
  - James M. Varanelli and James P. Cohoon. Two-stage simulated annealing methodology. In Proceedings of the 5th Great Lakes Symposium on VLSI, pages 50–53, Buffalo, NY, 16.-18. March 1995. IEEE, Los Alamitos, CA. †EI M153001/95
  - K. D. Boese A. B. Kahng and S. Muddu A new adaptive multistart technique for combinatorial global optimizations, Operations Research Letters, 16(2): 101-113, September, 1993.
