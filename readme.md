## Description
C++ annealer for simple pcb placement of polygonal components.
Supports the following:
 - Analytical overlap for arbitrary polygons using boost geometries
 - HPWL cost term for wirelength
 - BEN-AMEUR et al. cost-normalization & automatic initial temperature
 - Variable placement shift window, smaller displacement with temperature
 - Geometric temperature updates (cooling schedule) based on timberwolf schedule
 - Parallel multistart
 - Bookshelf parser
 - Bookshelf version translator
 - Plotting & animations
 - rtree spatial indexing for fast overlap computation

 ## Setup and requirements
 Setup:
 - Clone repository from github
 - Run make install to create necessary directories
 
 Requirements:
 - Boost version >=1.62
 - python3
 - shapely
 - matplotlib
 - docopt
 - numpy

 ## Files & directories
 - benchmarks/*: benchmarks in new-bookshelf format
 - cpp-taskflow/ & taskflow/: headerfiles for cpp-taskflow concurrency library [currently bugged]
 - main.h/cpp: Contains annealer logic and parameter parsing
 - readFiles.h/cpp: Contains pcb component, pin, terminal structs & bookshelf parser
 - readScl.h/cpp: SCL file parser
 - multistart.py: Multistart/GWTW driver for CPP annealer. Temporary replacement for cpp-taskflow which is broken due to recent updates.
 - load_bookshelf.py/utils.py: Implements bookshelf parser & utilities
 - tranlsate.py: Implements translator from older version of bookshelf to newer.
 - make_plots.py: Implements batch animation for placements and routability, pl plotter, and cost plotter

## Example
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
    
    EXAMPLE: ./sa -i 20000 -j 20 -t 1 -p input -f output

## Usage Details
 - Boardfile is the only mandatory option.
 - Partial placements in cache dir.
 - Cost over time in reports dir.
 - Full placement in base dir, uses index option.

## TODO
 - Reorganize repository to facilitate better management
 - Restructure annealer using better code practice - i.e. reduce usage of global variables (ABKCommon, UMPack, UTD-BoxRouter, REplace)
 - Fix experimental framework - i.e. Make it easier to generate plots of experiment statistics
 - Fix cpp-taskflow based multistart & GWTW
 - Figure out a better method for tuning parameters
 - Fix support for weighted nodes
