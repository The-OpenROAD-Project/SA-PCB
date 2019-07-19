# eagle2bookshelf
Utilities for converting between EAGLE and bookshelf formats.

## eagle2bookshelf.py

This program converts EAGLE board files (.brd) to bookshelf format files.

Specifically this program outputs block component files (.blocks), netlist files (.nets), and net weight files (.wts).
These files are sutible for academic IC placement programs.

The current does not account for net weights. All weights are set to '1'.

Accounting for non-90 degree rotations is in progress.
