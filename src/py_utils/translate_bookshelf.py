"""
///////////////////////////////////////////////////////////////////////////////
// Authors: Chester Holtz, Devon Merrill, James (Ting-Chou) Lin, Connie (Yen-Yi) Wu
//          (respective Ph.D. advisors: Chung-Kuan Cheng, Andrew B. Kahng, Steven Swanson).
//
// BSD 3-Clause License
//
// Copyright (c) 2018, The Regents of the University of California
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
//   list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holder nor the names of its
//   contributors may be used to endorse or promote products derived from
//   this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
///////////////////////////////////////////////////////////////////////////////
"""

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
