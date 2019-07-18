import pcbnew  
from pcbnew import IO_MGR, LoadBoard, PCB_LAYER_ID_COUNT
board = pcbnew.LoadBoard(filename, pcbnew.IO_MGR.EAGLE)
board.Save(fname)


