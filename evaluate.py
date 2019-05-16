import sys, math
import numpy as np
import shapely
from shapely.geometry import Point, Polygon

import load_bookshelf
import utils

def euclidean(components, board_pins, nets):
	"""
	Compute sum euclidean distance over all nets
	"""
	euwl = 0
	for net in nets:
		for i,pin1 in enumerate(net):
			if pin1[0] not in board_pins or 'cc' in pin1[0] or 'clk' in pin1[0]:
				pin1x,pin1y = utils.pin_pos(pin1,components)
			else:
				pin1x = pin1[1].x
				pin1y = pin1[1].y
			for j,pin2 in enumerate(net,i+1): # for each pin in connection
				if pin2[0] not in board_pins or 'cc' in pin2[0] or 'clk' in pin2[0]:
					pin2x,pin2y = utils.pin_pos(pin2,components)
				else:
					pin2x = pin2[1].x
					pin2y = pin2[1].y
				euwl += np.sqrt(np.square(pin1x - pin2x) + np.square(pin1y - pin2y))/(len(net)-1)
	return int(euwl)

def hpwl(components, board_pins, nets):
	"""
	Compute sum half-perimeter wirelength over all nets
	"""
	hpwl = 0
	for net in nets:
		plxs = []
		plys = []
		for pin in net: # for each pin in connection
			#if pin[0] not in board_pins or 'cc' in pin[0] or 'clk' in pin[0]:
			if pin[0] in components:
				pinx,piny = utils.pin_pos2(pin,components,comp2rot)
			else:
				pinx = pin[1].x
				piny = pin[1].y
			plxs.append(math.floor(pinx))
			plys.append(math.floor(piny))

		yd = max(plys) - min(plys)
		xd = max(plxs) - min(plxs)
		hpwl += (yd + xd)
		#print(str([p[0] for p in net]),': ', yd + xd)
	return int(hpwl)

def manhattan(components, board_pins, nets):
	"""
	Compute sum manhattan distance over all nets
	"""
	mwl = 0
	for net in nets:
		hist = set()
		for i,pin1 in enumerate(net):
			if pin1[0] not in board_pins or 'cc' in pin1[0] or 'clk' in pin1[0]:
				pin1x,pin1y = utils.pin_pos(pin1,components)
			else:
				pin1x = pin1[1].x
				pin1y = pin1[1].y
			for j,pin2 in enumerate(net,i+1): # for each pin in connection
				if pin2[0] not in board_pins or 'cc' in pin2[0] or 'clk' in pin2[0]:
					pin2x,pin2y = utils.pin_pos(pin2,components)
				else:
					pin2x = pin2[1].x
					pin2y = pin2[1].y
				mwl += np.abs(pin1x - pin2x) + np.abs(pin1y - pin2y)/(len(net)-1)
	return int(mwl)

#nodesfile = sys.argv[1]
#plfile = sys.argv[2]
#netsfile = sys.argv[3]

circuitname = './benchmarks/pcb_benchmark_devel-master/bm1'
plfile = circuitname + '.pl'
nodesfile = circuitname + '.nodes'
netsfile = circuitname + '.nets'
board_pins = {}

components = load_bookshelf.read_nodes(nodesfile)
components,comp2rot,_,board_pins = load_bookshelf.read_pl2(plfile,components)
nets,mod2net = load_bookshelf.read_nets2(netsfile,components,board_pins)
print(nets)
#eu_wl = euclidean(components, board_pins, nets)
hpwl = hpwl(components, board_pins, nets)
#mh_wl = manhattan(components, board_pins, nets)

#print('euclidean: ' + str(eu_wl))
print('hpwl: ' + str(hpwl))
#print('manhattan: ' + str(mh_wl))
