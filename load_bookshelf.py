import sys
import yaml
import json
import numpy as np
import math
import ast
import utils
import shapely
from shapely.geometry import Point, Polygon
from shapely.affinity import translate, rotate

import argparse

def read_pl(fname):
	"""
	Read & parse .pl (placement) file
	:param fname: .pl filename
	"""
	with open(fname,'r') as f:
		lines = f.read().splitlines()
	#lines = lines[5:]
	lines = lines[4:]
	components = {}
	board_pins = {}
	bp = 0
	for line in lines:
		if line == '':
			bp = 1
			continue
		if bp == 0:
			l = line.split()
			cname = l[0]
			dims = l[5]+l[6]
			dims = [float(i.strip()) for i in dims[1:-1].split(",")]
			x = float(l[1].strip())
			y = float(l[2].strip())
			poly = Polygon([[x,y], \
						   [x + dims[0],y], \
						   [x + dims[0],y+dims[1]], \
						   [x,y+dims[1]]])
			components[cname] = poly
		else:
			l = line.split()
			pname = l[0]
			coords = Point([float(l[1].strip()),float(l[2].strip())])
			board_pins[pname] = coords

	return components, board_pins

def read_pl2(fname,components):
	"""
	Read & parse .pl (placement) file
	:param fname: .pl filename
	"""
	with open(fname,'r') as f:
		lines = f.read().splitlines()

	lines = lines[5:]
	bp = 0
	comp2rot = {}
	board_pins = {}
	static_components = []
	for line in lines:
		if line.strip() == '':
			bp = 1
			continue
		else:
			if bp == 0:
				l = line.split()
				pname = l[0].strip()
				minx, miny, maxx, maxy = components[pname].bounds
				newx,newy = (float(l[1]),float(l[2]))
				r = l[4]
				if len(l) > 5:
					if '/FIXED' in l[5]:
						static_components.append(pname)
				comp2rot[pname] = r
				rot2deg = {'N':0,'S':180,'E':90,'W':270}
				components[pname] = rotate(components[pname],rot2deg[r])
				minx, miny, maxx, maxy = components[pname].bounds
				components[pname] = translate(components[pname],newx-minx,newy-miny)
			else:
				l = line.split()
				pname = l[0]
				coords = Point([float(l[1].strip()),float(l[2].strip())])
				board_pins[pname] = coords
	return components, comp2rot, board_pins, static_components


def read_nets(fname,components,board_pins):
	"""
	Read & parse .nets (netlist) file
	:param fname: .nets filename
	"""
	nets = []
	mod2net = {} # {component: [nets]}
	i = -1
	with open(fname,'r') as f:
		lines = f.read().splitlines()
	target_ibdex = [i for i, s in enumerate(lines) if 'NetDegree' in s][0]
	lines = lines[target_ibdex:]
	t = -1
	for line in lines:
		if line[0] == '#':
			continue
		if 'NetDegree' in line:
			i += 1
			nets.append([])
			continue

		l = line.split()
		if len(l) < 3:
			pin_name = l[0]
			pin = board_pins[pin_name]
			nets[i].append([pin_name,pin])
		else:
			pin_name = l[0]
			#if t == 0:
			local_pin_loc = [float(l[3].replace('%',''))/100.0 + 0.5, float(l[4].replace('%',''))/100.0 + 0.5]
			#pinx,piny = utils.pin_pos([pin_name, local_pin_loc], components)
			#print(pin_name, local_pin_loc, pinx, piny)
			#pin = Point([pinx,piny])
			"""
			else:
				minx,miny,maxx,maxy = components[pin_name].bounds
				w = (maxx - minx)
				h = (maxy - miny)
				pinx = (float(l[3][1:]) - minx)/(100*w) + 0.5
				piny = (float(l[4][1:]) - miny)/(100*h) + 0.5
				local_pin_loc = [pinx, piny]
			"""
			pin = local_pin_loc
			nets[i].append([pin_name, pin])
		if pin_name in mod2net:
			mod2net[pin_name].append(i)
		else:
			mod2net[pin_name] = [i]
	return nets, mod2net

def read_nets2(fname,components,board_pins):
	"""
	Read & parse .nets (netlist) file
	:param fname: .nets filename
	"""
	nets = []
	mod2net = {} # {component: [nets]}
	i = -1
	with open(fname,'r') as f:
		lines = f.read().splitlines()
	target_ibdex = [i for i, s in enumerate(lines) if 'NetDegree' in s][0]
	lines = lines[target_ibdex:]

	t = -1
	for line in lines:
		if line[0] == '#':
			continue
		if 'NetDegree' in line:
			i += 1
			nets.append([])
			continue

		l = line.split()
		#if len(l) < 3 or l[0] in board_pins:
		if len(l) < 3:
			pin_name = l[0]
			pin = board_pins[pin_name]
			nets[i].append([pin_name,pin])
		else:
			pin_name = l[0]
			local_pin_loc = [float(l[3]), float(l[4])]
			pin = local_pin_loc
			nets[i].append([pin_name, pin])
		if pin_name in mod2net:
			mod2net[pin_name].append(i)
		else:
			mod2net[pin_name] = [i]
	return nets, mod2net

def read_blocks(fname):
	"""
	Read & parse .blocks file
	:param fname: .blocks filename
	"""
	blocks = {}
	with open(fname, 'r') as f:
		lines = f.read().splitlines()
	target_ibdex = [i for i, s in enumerate(lines) if 'NumTerminals' in s][0]
	lines = lines[target_ibdex+2:]

	components = {}
	bp = 0
	for line in lines:
		if line[0] == '#':
			continue
		if line == '':
			bp = 1
			continue
		if bp == 0:
			l = line.split()
			cname = l[0]
			vstring =  ' '.join(l[3:])
			vstring = '[' + vstring.replace(') (', '),(') + ']'
			vertices = ast.literal_eval(vstring)
			poly = Polygon(vertices)
			components[cname] = poly

	return components

def read_nodes(fname):
	"""
	Read & parse .nodes (blocks) file
	:param fname: .nodes filename
	"""
	blocks = {}
	with open(fname, 'r') as f:
		lines = f.read().splitlines()

	target_ibdex = [i for i, s in enumerate(lines) if 'NumTerminals' in s][0]
	lines = lines[target_ibdex+2:]
	components = {}
	bp = 0
	for line in lines:
		if line == '':
			bp = 1
			continue
		if bp == 0:
			l = line.split()
			cname = l[0]
			w = float(l[1])
			h = float(l[2])
			poly = Polygon([[0,0],[w,0],[w,h],[0,h]])
			#print(l)
			#print(cname, w, h)
			#print(poly.exterior)
			components[cname] = poly
	return components

def write_pl(fname,components,board_pins):
	with open(fname,'w') as f:
		f.write('UMICH blocks 1.0\n')
		f.write('\n')
		f.write('\n')
		f.write('\n')
		for cname in components:
			component = components[cname]
			f.write(cname)
			f.write('\t')
			minx,miny,maxx,maxy = component.bounds
			f.write(str(minx))
			f.write('\t')
			f.write(str(miny))
			f.write('\t')
			f.write('DIMS = (' + str(maxx - minx) + ', ' + str(maxy - miny) + ')')
			f.write(' : N\n')

		f.write('\n')

		for pname in board_pins:
			if pname in components:
				pass
			pin = board_pins[pname]
			f.write(pname)
			f.write('\t')
			f.write(str(pin.x))
			f.write('\t')
			f.write(str(pin.y))
			f.write('\t')
			f.write(' : N\n')

def write_nodes(fname,components,board_pins):
	with open(fname,'w') as f:
		f.write('UCLA blocks 1.0\n')
		f.write('\n')
		f.write('\n')
		f.write('\n')
		f.write('NumNodes\t :\t ' + str(len(set([c for c in components]))))
		f.write('\n')
		f.write('NumTerminals\t :\t ' + str(len(set([b for b in board_pins if b not in components]))))

		f.write('\n')

		for cname in components:
			component = components[cname]
			f.write(cname)
			f.write('\t')
			minx,miny,maxx,maxy = component.bounds
			f.write(str((maxx-minx)/100))
			f.write('\t')
			f.write(str((maxy-miny)/100))
			f.write('\n')
		for pname in board_pins:
			board_pin = board_pins[pname]
			f.write(pname)
			f.write('\t')
			f.write(str(1))
			f.write('\t')
			f.write(str(1))
			f.write('\t')
			f.write('terminal_NI')
			f.write('\n')

def write_newnets(fname,nets,components):
	with open(fname,'w') as f:
		f.write('UCLA nets 1.0\n')
		f.write('\n')
		f.write('\n')
		f.write('\n')
		f.write('NumNets\t :\t ' + str(len(nets)))
		f.write('\n')
		f.write('NumPins\t :\t ' + str(len(set([pin[0] for net in nets for pin in net]))))

		f.write('\n')

		for net in nets:
			f.write('NetDegree : ' + str(len(net)) + '\n')
			for pin in net:
				pinname, pinloc = pin
				if isinstance(pinloc, list): # component pin
					centroid = components[pinname].centroid
					pinx,piny = utils.pin_pos(pin, components)
					f.write(pinname + ' B : ' + str(round(pinx - centroid.x,2)/100) + ' ' +  str(round(piny - centroid.y,2)/100) + '\n')
				else: # terminal pin
					f.write(pinname + ' B \n')

def write_newpl(fname,components,board_pins):
	with open(fname,'w') as f:
		f.write('UCLA pl 1.0\n')
		f.write('\n')
		f.write('\n')
		f.write('\n')
		for cname in components:
			component = components[cname]
			f.write(cname)
			f.write('\t')
			minx,miny,maxx,maxy = component.bounds
			f.write(str(minx/100))
			f.write('\t')
			f.write(str(miny/100))
			f.write('\t')
			#f.write('DIMS = (' + str(maxx - minx) + ', ' + str(maxy - miny) + ')')
			f.write(': N\n')

		f.write('\n')

		for pname in board_pins:
			pin = board_pins[pname]
			f.write(pname)
			f.write('\t')
			f.write(str(pin.x/100))
			f.write('\t')
			f.write(str(pin.y/100))
			f.write('\t')
			f.write(' : N\t')
			f.write('terminal_NI\n')

#=blocksfile = '/Users/orange3xchicken/Downloads/merrill_place_example_1.blocks'
#blk = read_blocks(blocksfile)
#print(blk)
#print(len(blk))
#plfile = sys.argv[1]
#netsfile = sys.argv[2]
#blocksfile = sys.argv[3]
"""
components:        dictionary of placed component-polygons indexed by their name
placed_components: dictionary of placed component-polygons indexed by their name
board_pins:        dictionary of static point-pins indexed by their name
nets:              list of k nets, each net is a list of n point-pins
"""
#components = read_blocks(blocksfile)
#placed_components, board_pins = read_pl(plfile)
#nets = read_nets(netsfile, components, board_pins)

#print('pcb data loaded')
