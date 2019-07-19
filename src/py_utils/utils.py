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

from matplotlib import patches
import matplotlib.pyplot as plt

import shapely
from shapely import geometry

import numpy as np

def plot_circuit(circuit_name, components, nets, board_dim, stats = None):
	board_lower_left_corner = [min(board_dim[0]), min(board_dim[1])]
	width_height = [max(board_dim[0]) - min(board_dim[0]), \
				    max(board_dim[1]) - min(board_dim[1])]

	fig, ax = plt.subplots(1)
	if stats is None:
		ax.set_title(circuit_name)
	else:
		ax.set_title("circuit: " + circuit_name+ " wirelength: " + str(round(stats[0],2)) + " overlap: " + str(round(stats[1],2)))
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
			points = np.array([x, y], np.int32).T
			polygon_shape = patches.Polygon(points, linewidth=1, edgecolor='r', facecolor='none')
			ax.add_patch(polygon_shape)
			label = ax.annotate(name, xy=(c.x, c.y), fontsize=5, ha="right", va="top")

	for net in nets:
		netlist = []
		for pin in net:
			if pin[0] not in components:
				cx = pin[1].x
				cy = pin[1].y
			else:
				cx, cy = utils.pin_pos(pin,components)
			netlist.append([cx,cy])
		#for i in range(len(netlist)-1):
		#	ax.plot([netlist[i][0],netlist[i+1][0]],[netlist[i][1],netlist[i+1][1]], color="blue", linewidth=1, alpha=0.25, linestyle='dashed')
		#xs= [ x[0] for x in netlist ]
		#ys= [ x[1] for x in netlist ]
		#ax.scatter(xs,ys,marker='.')
	plt.xlim(-25, max(width_height) + 25)
	plt.ylim(-25, max(width_height) + 25)
	plt.gca().set_aspect('equal', adjustable='box')
	plt.show()

def pin_pos2(pin_loc, modules,comp2rot):
	"""
	Convert localized pin positions to position wrt
	 global coordinates
	:param pin_loc: pin location of the form [pinname, [%x, %y]]
	:param modules: list of modules
	"""
	module_name, local_pin_loc = pin_loc
	cx = modules[module_name].centroid.x
	cy = modules[module_name].centroid.y
	r = comp2rot[module_name]
	if r == 'N':
		pinx = cx + local_pin_loc[0]
		piny = cy + local_pin_loc[1]
	elif r == 'S':
		pinx = cx - local_pin_loc[0]
		piny = cy - local_pin_loc[1]
	elif r == 'E':
		pinx = cx - local_pin_loc[1]
		piny = cy + local_pin_loc[0]
	elif r == 'W':
		pinx = cx + local_pin_loc[1]
		piny = cy - local_pin_loc[1]

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
