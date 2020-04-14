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

import matplotlib
#matplotlib.use('Agg')
from matplotlib import patches
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import sys

import shapely
from shapely import geometry

import numpy as np
import load_bookshelf

def plot_circuit(circuit_name, components,radii, comp2rot, comp2layer, nets, board_dim,idx=0,figname=None, stats = None, ax=None):
        """
        board dim: [[xs],[ys]]
        """
        board_lower_left_corner = [min(board_dim[0]), min(board_dim[1])]
        board_upper_right_corner = [max(board_dim[0]), max(board_dim[1])]
        if ax is None:
                fig, ax = plt.subplots(1)
        else:
                ax.clear()
        #if stats is None:
        #       ax.set_title(circuit_name)
        #else:
                #ax.set_title("circuit: " + circuit_name+ " wirelength: " + str(round(stats[0],2)) + " overlap: " + str(round(stats[1],2)))
        ax.set_title(idx)
        boundary = patches.Rectangle((min(board_dim[0]), min(board_dim[1])), \
                                                                  max(board_dim[0]) - min(board_dim[0]), max(board_dim[1]) - min(board_dim[1]), \
                                                                  linewidth=1,edgecolor='b',facecolor='none')
        ax.add_patch(boundary)
        for name,shape in components.items():
                #print(shape.bounds)
                if name not in ['C31']:
                    continue
                if isinstance(shape, geometry.Point):
                        x = shape.x
                        y = shape.y
                        ax.scatter(x,y, c='b',s=5)
                        label = ax.annotate(name, xy=(x, y), fontsize=5, ha="right", va="top", )
                else:
                        x, y = shape.exterior.coords.xy
                        c = shape.centroid
                        x = [i for i in x]
                        y = [i for i in y]
                        points = np.array([x, y], np.float32).T
                        if comp2layer[name] == 1:
                            polygon_shape = patches.Polygon(points, linewidth=1, edgecolor='r', facecolor='r',alpha=0.5)
                        else:
                            polygon_shape = patches.Polygon(points, linewidth=1, edgecolor='g', facecolor='g',alpha=0.5)
                        ax.add_patch(polygon_shape)
                        radius = radii[name]

                        hline = Line2D([c.x-radius,c.x+radius],[c.y]*2,linewidth=0.2)
                        vline = Line2D([c.x]*2,[c.y-radius,c.y+radius], linewidth=0.2)
                        ax.add_line(hline)
                        ax.add_line(vline)
                        label = ax.annotate(name, xy=(c.x, c.y), fontsize=5, ha="right", va="top")

        # draw nets
        """
        for net in nets:
                netlist = []
                c = [np.random.rand(3)] # random colors
                for pin in net:
                        if pin[0] not in components: #or pin[0] in board_pins:
                                try:
                                        cx = pin[1].x
                                        cy = pin[1].y
                                except:
                                        cx = pin[1][0]
                                        cy = pin[1][1]
                        else:
                                cx, cy = pin_pos2(pin,components,comp2rot)
                                ##cx, cy = pin_pos(pin,components)
                                #label = ax.annotate(pin, xy=(cx, cy), fontsize=5, ha="right", va="top", )
                        netlist.append([cx,cy])
                xmax = max([p[0] for p in netlist])
                xmin = min([p[0] for p in netlist])
                ymax = max([p[1] for p in netlist])
                ymin = min([p[1] for p in netlist])
                center =  [(xmax + xmin)/2,(ymax + ymin)/2]
                # centroid - star
                #for i in range(len(netlist)):
                #       ax.plot([netlist[i][0],center[0]],[netlist[i][1],center[1]], color=tuple(map(tuple, c))[0] + (255,), linewidth=1, alpha=0.5, linestyle='dashed')
                #xs= [ x[0] for x in netlist ]
                #ys= [ x[1] for x in netlist ]
                #ax.scatter(xs,ys,marker='.',c=c)
                #ax.scatter(center[0],center[1],marker='.',c=c)
        """
        #plt.xlim(board_lower_left_corner[0] - 5,board_upper_right_corner[0] + 5)
        #plt.ylim(board_lower_left_corner[1] - 5,board_upper_right_corner[1] + 5)
        ax.set_xlim(board_lower_left_corner[0] - 5, board_upper_right_corner[0] + 5)
        ax.set_ylim(board_lower_left_corner[1] - 5, board_upper_right_corner[1] + 5)

        #plt.gca().set_aspect('equal', adjustable='box')
        print(figname)
        fig.savefig(figname)
        #plt.close()
        #plt.show()

def pin_pos2(pin_loc, modules,comp2rot):
        """
        Convert localized pin positions to position wrt
         global coordinates
        :param pin_loc: pin location of the form [pinname, [xoffset, yoffset]]
        :param modules: list of modules
        """
        module_name, local_pin_loc = pin_loc
        cx = modules[module_name].centroid.x
        cy = modules[module_name].centroid.y
        if module_name in comp2rot:
                r = comp2rot[module_name]
        else:
                r = 'N'
        if r == 'N':
                pinx = cx + local_pin_loc[0]
                piny = cy + local_pin_loc[1]
        elif r == 'S':
                pinx = cx - local_pin_loc[0]
                piny = cy - local_pin_loc[1]
        elif r == 'E':
                pinx = cx + local_pin_loc[1]
                piny = cy - local_pin_loc[0]
        elif r == 'W':
                pinx = cx - local_pin_loc[1]
                piny = cy + local_pin_loc[0]

        return pinx, piny

def pin_pos(pin_loc, modules):
        """
        Convert localized pin positions to position wrt
         global coordinates
        :param pin_loc: pin location of the form [pinname, [%x, %y]]
        :param modules: list of modules
        """
        module_name, local_pin_loc = pin_loc
        minx, miny, maxx, maxy = modules[module_name].bounds

        pinx = (maxx - minx) * local_pin_loc[0] + minx
        piny = (maxy - miny) * local_pin_loc[1] + miny
        return pinx, piny

cachedir = '../cache2/'
nodesfile = cachedir+'nodes.nodes'
figname = cachedir+'anim.gif'

boarddimmin = None
boarddim = None

board_pins = {}

from os import walk
from tqdm import tqdm

f = []
for (dirpath, dirnames, filenames) in walk(cachedir):
    if 'anim.gif' in filenames:
        filenames.remove('anim.gif')
    if 'nodes.nodes' in filenames:
        filenames.remove('nodes.nodes')
    f.extend(filenames)
f = [fname for fname in f if '.pl' in fname]
f = sorted(f, key=lambda x: int(x.split('.')[0]))


for i,ff in enumerate(tqdm(f)):
    components = load_bookshelf.read_nodes(nodesfile)
    for c in components:
        shape = components[c]
    
    if '.pl' in ff:
        components,comp2rot,board_pins,_,comp2layer = load_bookshelf.read_pl2(cachedir+ff,components)
        radii = load_bookshelf.read_rad(cachedir+ff.split('.')[0]+'.rad')
    if '.rad' in ff:
        continue
    elif '.txt' in ff:
        continue
    elif '.png' in ff:
        continue
    elif '.gif' in ff:
        continue
    elif '.nodes' in ff:
        continue
    #nets,mod2net = load_bookshelf.read_nets2(netsfile,components,board_pins)
    nets = None
    #board_dim = [[138-20,158+20],[97-20,111+20]]
    board_dim = [[95-20,200+20],[95-20,160+20]]
    #if i % 10 == 0:
    plot_circuit(ff.split('.')[0], components,radii,comp2rot,comp2layer,nets,board_dim,ff.split('.')[0],cachedir+'img/'+ff.split('.')[0]+'.png')


f = []
for (dirpath, dirnames, filenames) in walk(cachedir+'img/'):
    f.extend(filenames)
    
f = sorted(f, key=lambda x: int(x.split('.')[0]))
import imageio
with imageio.get_writer(figname, mode='I') as writer:
        for ff in tqdm(f):
                try:
                        if ff == "" or ff == "anim":
                                continue
                        ff = cachedir+'img/' + ff
                        image = imageio.imread(ff)
                        writer.append_data(image)
                except:
                        tqdm.write("exception: " + ff)
                        continue
