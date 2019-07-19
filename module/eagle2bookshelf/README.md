# eagle2bookshelf
This is a project that converts between EAGLE PCB files to bookshelf fileformat.

## eagle2bookshelf.py

This program converts EAGLE board files (.brd) to bookshelf format files.

Specifically this program outputs block component files (.blocks), netlist files (.nets), and net weight files (.wts).
These files are sutible for academic IC placement programs.

The current version does not account for pin placement. All nets are connected to the center of each block.

The current does not account for net weights. All weights are set to '1'.

```
Usage:
  eagle2bookshelf.py -h | --help
  eagle2bookshelf.py --brd <BRD> --output_prfx <STEM_NAME> --userid <USERID>
-h --help                      Show this message.
-i --brd BRD                   The EAGLE .brd file to convert.
-o --output_prfx STEM_NAME     The stem name for the new files. Includes directory.
--userid USERID                Your name and contact.
```
