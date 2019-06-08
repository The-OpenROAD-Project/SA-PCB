## Description
Timberwolf-based C++ annealer for simple pcb placement of polygonal components.
Supports the following:
 - Analytical overlap for arbitrary polygons using boost geometries
 - Clique-weighted wirelength & HPWL
 - Varanelli-cohoon cost-normalization & initial temperature
 - Variable shift window facilitates convergence of placement
 - Geometric temperature updates based on timberwolf schedule
 - Parallel multistart
 - Bookshelf parser
 - Bookshelf version translator
 - Plotting & animations

 ## Files & directories
 - benchmarks/*: benchmarks in new-bookshelf format
 - cpp-taskflow/ & taskflow/: headerfiles for cpp-taskflow concurrency library [currently bugged]
 - main.h/cpp: Contains annealer logic and parameter parsing
 - readFiles.h/cpp: Contains pcb component, pin, terminal structs & bookshelf parser
 - readScl.h/cpp: SCL file parser
 - multistart.py: Multistart/GWTW driver for CPP annealer. Temporary replacement for cpp-taskflow which is broken due to recent updates.
 - load_bookshelf.py/utils.py: Implements bookshelf parser & utilities
 - tranlsate.py: Implements translator from older version of bookshelf to newer.
 - plot.py: Implements batch animation & pl plotter

## Example
    "./sa usuage is \n \
                                    -i <value> : for denoting # outer iterations PER SA INSTANCE \n \
									-x <value>: denotes idx for parallel multistart
                                    -j <value> : for denoting 'j'*#nodes inner iterations \n \
                                    -t <value> : for denoting initial temperature \n \
                                    -k <value> : for denoting parallel instances \n \
                                    -s <value> : for denoting splits every floor(i/s) iterations. Use s=0 for 'k' independent instances of SA \n \
                                    -w <value> : for denoting number of winners per splitting \n \
                                    -f <str> : for denoting number output pl \n \
                                    EXAMPLE: ./sa -i 20000 -j 20 -t 40000 -k 0 -s -w 0 -f output.pl for the standard single instance timberwolf algorithm"

## TODO
 - Reorganize repository to facilitate better management
 - Restructure annealer using better code practice - i.e. reduce usage of global variables (ABKCommon, UMPack, UTD-BoxRouter, REplace)
 - Fix experimental framework - i.e. Make it easier to generate plots of experiment statistics
 - Fix cpp-taskflow based multistart & GWTW
 - Figure out a better method for tuning parameters
