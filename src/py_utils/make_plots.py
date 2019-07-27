"""make_plots.
This program makes plots
Usage:
  make_plots.py -h | --help
  make_plots.py --brd <BRD> --pl <PL> --out <OUT_NAME> --reports <reports>

-h --help                      Show this message.
-i --brd BRD                   The circuit file.
-o --out OUT
-p --pl PL
-r --reports reports
-apl --anim_pl animate_reports
-arb --anim_rb animate_routability
"""

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
matplotlib.use('Agg')
from matplotlib import patches
import matplotlib.pyplot as plt
import sys

from os import walk
from tqdm import tqdm

import imageio
import shapely
from shapely import geometry

import numpy as np
import load_bookshelf
from docopt import docopt

def plot_circuit(circuit_name, components, comp2rot, nets, board_dim,figname=None, ax=None):
    """
    board dim: [[xs],[ys]]
    """
    board_lower_left_corner = [min(board_dim[0]), min(board_dim[1])]
    board_upper_right_corner = [max(board_dim[0]), max(board_dim[1])]
    if ax is None:
        fig, ax = plt.subplots(1)
    else:
        ax.clear()
    boundary = patches.Rectangle((min(board_dim[0]), min(board_dim[1])), \
                                  max(board_dim[0]) - min(board_dim[0]), max(board_dim[1]) - min(board_dim[1]), \
                                  linewidth=1,edgecolor='b',facecolor='none')
    ax.add_patch(boundary)
    for name,shape in components.items():
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
            polygon_shape = patches.Polygon(points, linewidth=1, edgecolor='r', facecolor='r',alpha=0.5)
            ax.add_patch(polygon_shape)
            label = ax.annotate(name, xy=(c.x, c.y), fontsize=5, ha="right", va="top")

    # draw nets
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
                cx, cy = pin_pos(pin,components,comp2rot)
            netlist.append([cx,cy])
        xmax = max([p[0] for p in netlist])
        xmin = min([p[0] for p in netlist])
        ymax = max([p[1] for p in netlist])
        ymin = min([p[1] for p in netlist])
        center =  [(xmax + xmin)/2,(ymax + ymin)/2]
        # centroid - star
        for i in range(len(netlist)):
            ax.plot([netlist[i][0],center[0]],[netlist[i][1],center[1]], color=tuple(map(tuple, c))[0] + (255,), linewidth=1, alpha=0.5, linestyle='dashed')
        xs= [ x[0] for x in netlist ]
        ys= [ x[1] for x in netlist ]
        ax.scatter(xs,ys,marker='.',c=c)
        ##ax.scatter(center[0],center[1],marker='.',c=c)
    plt.xlim(board_lower_left_corner[0] - 5,board_upper_right_corner[0] + 5)
    plt.ylim(board_lower_left_corner[1] - 5,board_upper_right_corner[1] + 5)
    #plt.gca().set_aspect('equal', adjustable='box')
    plt.savefig(figname)
    plt.close()

def pin_pos(pin_loc, modules,comp2rot):
    """
    Convert localized pin positions to position wrt
     global coordinates
    :param pin_loc: pin location of the form [pinname, [xoffset, yoffset]]
    :param modules: list of modules
    """
    module_name, local_pin_loc = pin_loc
    cx = modules[module_name].centroid.x
    cy = modules[module_name].centroid.y
    """
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
    """
    pinx = cx
    piny = cy
    return pinx, piny

def plot_pl(plfile, nodesfile, netsfile, figname, shapesfile=None, board_dim=None):
    board_pins = {}
    components = load_bookshelf.read_nodes(nodesfile)
    if shapesfile:
        components = load_bookshelf.read_shapes(shapesfile,components)
    components,comp2rot,board_pins,_ = load_bookshelf.read_pl2(plfile,components)
    nets,mod2net = load_bookshelf.read_nets2(netsfile,components,board_pins)
    if board_dim is None:
        if board_pins is not None and len(board_pins) > 0:
            xs = [pin[1].x for pin in board_pins.items()]
            ys = [pin[1].y for pin in board_pins.items()]
            board_dim = [xs,ys]
            pass
        else:
            #board_dim = [[-25,120],[-25,120]]
            #board_dim = [[0,25],[0,10]] #snake
            #board_dim = [[-10,120],[-8,55]] # bm1
            board_dim = [[-5,600],[-5,150]] # mt
    else:
        board_dim = board_dim#[[boarddimmin,boarddimmax],[boarddimmin, boarddimmax]]
    plot_circuit(plfile.split('.')[0], components,comp2rot,nets,board_dim,figname)

def make_heatmap_anim():
    routabilities = []
    fnames = []

    for (dirpath, dirnames, filenames) in os.walk('./cache_rudy'):
        fnames.extend(filenames)
    fnames = sorted(fnames, key=lambda x: float(x.split('.')[0]))

    for i,ff in enumerate(tqdm(fnames)):
         #mat = np.loadtxt('./cache_rudy/'+str(filename),delimiter='\t')
         mat = pd.read_table('./cache_rudy/'+str(ff),sep='\t',header=None).values[:,:-1]
         routabilities.append(mat)
    import seaborn as sns

    nx, ny = routabilities[0].shape

    fig = plt.figure()
    data = np.random.rand(nx, ny)
    sns.heatmap(data, vmax=.8, square=True)

    def init():
          sns.heatmap(np.zeros((nx, ny)), vmax=.8, square=True,cbar=False)

    def animate(i):
        data = routabilities[i]
        ax = sns.heatmap(data, vmax=.8, square=True,cbar=False)
        ax.invert_yaxis()


    plt.yticks = []
    plt.xticks = []

    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=None, repeat = True)
    plt.title('routability')

def make_placement_anim():
    f = []
    board_pins = {}
    for (dirpath, dirnames, filenames) in walk('./cache/'):
        f.extend(filenames)
    f = sorted(f, key=lambda x: float(x.split('.')[0]))

    for i,ff in enumerate(tqdm(f)):
        if i > 10000:
            break
        if i % 20 == 0:
            components = load_bookshelf.read_nodes(nodesfile)
            components,comp2rot,board_pins,_ = load_bookshelf.read_pl2('./cache/'+ff,components)
            nets,mod2net = load_bookshelf.read_nets2(netsfile,components,board_pins)
            board_dim = [[-5,100],[-5,55]]
            plot_circuit(ff.split('.')[0], components,comp2rot,nets,board_dim,'./cache/img/'+ff.split('.')[0]+'.png')

    f = []
    for (dirpath, dirnames, filenames) in walk('./cache/img/'):
        f.extend(filenames)
    f = sorted(f, key=lambda x: float(x.split('.')[0]))

    with imageio.get_writer(figname, mode='I') as writer:
        for ff in tqdm(f):
            try:
                if ff == "" or ff == "anim":
                    continue
                ff = './cache/img/' + ff
                image = imageio.imread(ff)
                writer.append_data(image)
            except:
                tqdm.write("exception: " + ff)
                continue

def make_graphs(reports_dir):
    ar_fname = reports_dir + ''
    cost_fname = reports_dir + ''
    wl_fname = reports_dir + ''
    oa_fname = reports_dir + ''

    wl = open(wl_fname).read().splitlines()
    oa = open(oa_fname).read().splitlines()
    cost = open(cost_fname).read().splitlines()
    ar = open(ar_fname).read().splitlines()

    wl = [float(w) for w in wl]
    oa = [float(o) for o in oa]
    cost = [float(c) for c in oa]
    ar = [float(a) for a in ar]

def main(brd_name, out_file, pl_file=None, pl_anim=False, heatmap_anim=False, reports_dir=False):
    circuitname = brd_name
    if pl_file:
        plfile = pl_file
    else:
        plfile = circuitname + '.pl'
    nodesfile = circuitname + '.nodes'
    netsfile = circuitname + '.nets'
    shapesfile = circuitname + '.shapes'
    figname = out_file

    bversion = 1
    boarddimmin = None
    boarddim = None

    #if reports_dir:
    #    make_graphs(reports_dir)
    if pl_file is not None:
        plot_pl(plfile, nodesfile, netsfile, figname, shapesfile)
    #if pl_anim:
    #    make_placement_anim()
    #if heatmap_anim:
    #    make_heatmap_anim()

if __name__ == "__main__":
    arguments = docopt(__doc__, version='make_plots v0.1')
    print(arguments)
    main(
        brd_name=str(arguments['--brd']),
        out_file=str(arguments['--out']),
        pl_file=str(arguments['--pl']),
        reports_dir=str(arguments['--reports'])
    )
