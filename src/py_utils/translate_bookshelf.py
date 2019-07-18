import load_bookshelf
import json

import argparse

def old2new(dir,cname,outdir):
	pl_fname = dir + cname + '.pl'
	nets_fname = dir + cname + '.nets'
	blocks_fname = dir + cname + '.blocks'

	components = load_bookshelf.read_blocks(blocks_fname)
	placed_components, board_pins = load_bookshelf.read_pl(pl_fname)
	nets, mod2net = load_bookshelf.read_nets(nets_fname, components, board_pins)

	load_bookshelf.write_newpl(outdir+circuit_name+'.pl',placed_components,board_pins)
	load_bookshelf.write_nodes(outdir+circuit_name+'.nodes',components,board_pins)
	load_bookshelf.write_newnets(outdir+circuit_name+'.nets',nets,components)

def new2old(dir,cname,outdir):
	pl_fname = dir + cname + '.pl'
	nets_fname = dir + cname + '.nets'
	nodes_fname = dir + cname + '.nodes'

	components = load_bookshelf.read_nodes(nodesfile)
	components,comp2rot,_,board_pins = load_bookshelf.read_pl2(plfile,components)
	nets,mod2net = load_bookshelf.read_nets2(netsfile,components,board_pins)

	load_bookshelf.write_pl(outdir+circuit_name+'.pl',placed_components,board_pins)
	load_bookshelf.write_blocks(outdir+circuit_name+'.nodes',components,board_pins)
	load_bookshelf.write_nets(outdir+circuit_name+'.nets',nets,components)

parser = argparse.ArgumentParser(description='Bookshelf version translator')

parser.add_argument('-v', action="store", default='old2new', dest='version')
parser.add_argument('-d', action="store", dest='dir')
parser.add_argument('-c', action="store", dest='cname')
parser.add_argument('-o', action="store", default='./' dest='outdir')

if parser.version == 'old2new':
	old2new(parser.dir, parser.cname)
else:
	new2old(old2new(parser.dir, parser.cname))
