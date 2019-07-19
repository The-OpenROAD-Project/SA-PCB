"""Eagle2Bookshelf2012.

This program converts EAGLE board files (.brd) to bookshelf (DAC 2012 contest flavor) format files.

Specifically this program outputs nodes files (.nodes), netlist files (.nets), and net weight files (.wts).
These files are sutible for academic IC placement programs.
This program also outputs a placement file (.pl) from the initial placement in the Eagle file.

This version DOES account for pin placement. But, this feature is unverified.
This version does not account for net weights. All weights are set to '1'.

This program was written by Devon Merrill (devon@ucsd.edu).

Usage:
  eagle2bookshelf2012.py -h | --help
  eagle2bookshelf2012.py --brd <BRD> --output_prfx <STEM_NAME> --userid <USERID>

-h --help                      Show this message.
-i --brd BRD                   The EAGLE .brd file to convert.
-o --output_prfx STEM_NAME     The stem name for the new files (file names without suffex). Includes directory.
--userid USERID                Your name and contact.
"""

LICENCE = """
BSD 3-Clause License

Copyright (c) 2015-2018, The Regents of the University of California
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

from __future__ import print_function


import datetime

import Swoop
from docopt import docopt


class PinPercentage(object):
	"""With with percentage based offset"""
	def __init__(self, name, x_offset, y_offset, direction):
		super(PinPercentage, self).__init__()
		self.name = name
		assert isinstance(x_offset, str)
		assert isinstance(y_offset, str)

		# these offsets are with respect to the (0,0) coordinate of the component
		self.x_offset = x_offset
		self.y_offset = y_offset
		self.direction = direction



class PinAbsolute(object):
	"""Pin with absolute offset in 'units' not percents"""
	def __init__(self, name, x_offset, y_offset, direction):
		super(PinAbsolute, self).__init__()
		self.name = name
		assert isinstance(x_offset, float) or isinstance(x_offset, int)
		assert isinstance(y_offset, float) or isinstance(y_offset, int)

		# these offsets are with respect to the (0,0) coordinate of the component
		self.x_offset = x_offset
		self.y_offset = y_offset
		self.direction = direction
		

class Signal(object):
	def __init__(self, name, weight=1):
		self.name = name
		self.pins = []
		self.weight = weight

	def add_pin_percentage(self, element, pin_name, direction='B'):
		origin_offsets = element.pins[pin_name] # the pin dict stores (x,y) pairs of origin offsets
		new_pin = PinPercentage(element.name, '%0.0', '%0.0', direction)
		self.pins.append(new_pin)

		# you have bounding box for element
		# you have x, y offsets for pin
		# what we want to do is center the bounding box
		# 	then center the pin
		#	then we can get the percent offset from center
		x_min = element.x_min
		x_max = element.x_max
		y_min = element.y_min
		y_max = element.y_max

		x_center = (x_min + x_max) / 2.0
		y_center = (y_min + y_max) / 2.0

		pin_x_from_center = origin_offsets[0] - x_center
		pin_y_from_center = origin_offsets[1] - y_center

		x_percent = pin_x_from_center / (x_max - x_min)
		y_percent = pin_y_from_center / (y_max - y_min)

		assert x_percent <= 0.5 and x_percent >= -0.5, 'x_percent: ' + str(x_percent)
		assert y_percent <= 0.5 and y_percent >= -0.5, 'y_percent: ' + str(y_percent)

		new_pin.x_offset = '%' + str(x_percent * 100)
		new_pin.y_offset = '%' + str(y_percent * 100)

	def add_pin_absolute(self, element, pin_name, direction='B'):
		origin_offsets = element.pins[pin_name] # the pin dict stores (x,y) pairs of origin offsets
		new_pin = PinAbsolute(element.name, 0.0, 0.0, direction)
		self.pins.append(new_pin)

		if element.name == 'JP9':
			print('Adding pin to JP9: ' + pin_name)

		x_min = element.x_min
		x_max = element.x_max
		y_min = element.y_min
		y_max = element.y_max

		if element.name == 'JP9':
			print('\tmin and max: ' + str( ((x_min,x_max),(y_min,y_max)) ))

		x_center = (x_min + x_max) / 2.0
		y_center = (y_min + y_max) / 2.0

		if element.name == 'JP9':
			print('\tcenter: ' + str((x_center, y_center)))
			print('\torigin offset: ' + str(origin_offsets))

		pin_x_from_center = origin_offsets[0] - x_center
		pin_y_from_center = origin_offsets[1] - y_center

		if element.name == 'JP9':
			print('\tfrom center pin location: ' + str((pin_x_from_center,pin_y_from_center)))

		new_pin.x_offset = pin_x_from_center
		new_pin.y_offset = pin_y_from_center

		# new_pin.x_offset = element.pins[pin_name][0]
		# new_pin.y_offset = element.pins[pin_name][1]


class ElementEntry(object):
	def __init__(self, name, library=None, package=None):
		self.name = name
		self.library = library
		self.package = package

		# these are not height and width, these are cartisian coordiants (can be negative)
		self.x_loc = 0.0
		self.y_loc = 0.0
		self.x_min = 9e99
		self.x_max = -9e99
		self.y_min = 9e99
		self.y_max = -9e99
		self.pins = {}
		self.bounding_box_multiplier = 1.0

	def __str__(self):
		"""This is in '.blocks' format"""
		width = self.x_max - self.x_min
		height = self.y_max - self.y_min
		width_expand = width # * self.bounding_box_multiplier / 2.0
		height_expand = height # * self.bounding_box_multiplier / 2.0
		box_str = ''
		box_str += '(' + str(self.x_min-width_expand) + ', ' + str(self.y_min-height_expand) + ')' + ', '
		box_str += '(' + str(self.x_min-width_expand) + ', ' + str(self.y_max+height_expand) + ')' + ', '
		box_str += '(' + str(self.x_max+width_expand) + ', ' + str(self.y_max+height_expand) + ')' + ', '
		box_str += '(' + str(self.x_max+width_expand) + ', ' + str(self.y_min-height_expand) + ')'
		return self.name + ' ' + 'hardrectilinear' + ' ' + '4' + ' ' + box_str

	def node_str(self):
		"""This is in '.nodes' format"""
		# print self.name, self.x_loc, self.y_loc
		width = self.x_max - self.x_min
		# width = 1#self.x_max - self.x_min
		height = self.y_max - self.y_min
		# height = 1#self.y_max - self.y_min
		node_str = ''
		node_str += self.name.rjust(20) + ' ' + str(width).rjust(20) + ' ' + str(height).rjust(20)
		return node_str

	def expand_bb(self, x_min, x_max, y_min, y_max):
		self.x_min = min(x_min, self.x_min)
		self.x_max = max(x_max, self.x_max)
		self.y_min = min(y_min, self.y_min)
		self.y_max = max(y_max, self.y_max)

	def add_pin(self, pin):
		if isinstance(pin, Swoop.Smd):
			self.pins[pin.get_name()] = (pin.get_x(), pin.get_y())
		elif isinstance(pin, Swoop.Pad):
			self.pins[pin.get_name()] = (pin.get_x(), pin.get_y())
		else:
			assert False



def de_bounding_box(drawing_element):
	"""Return a bounding box for the drawing element."""
	x_min = 9e99
	x_max = -9e99
	y_min = 9e99
	y_max = -9e99
	if isinstance(drawing_element, Swoop.Wire):
		if drawing_element.get_curve():
			pass
		else:
			x_min = min(drawing_element.get_x1(), drawing_element.get_x2()) - drawing_element.get_width()
			x_max = max(drawing_element.get_x1(), drawing_element.get_x2()) + drawing_element.get_width()
			y_min = min(drawing_element.get_y1(), drawing_element.get_y2()) - drawing_element.get_width()
			y_max = max(drawing_element.get_y1(), drawing_element.get_y2()) + drawing_element.get_width()
	elif isinstance(drawing_element, Swoop.Rectangle):
		x_min = min(drawing_element.get_x1(), drawing_element.get_x2()) # - drawing_element.get_width()
		x_max = max(drawing_element.get_x1(), drawing_element.get_x2()) # + drawing_element.get_width()
		y_min = min(drawing_element.get_y1(), drawing_element.get_y2()) # - drawing_element.get_width()
		y_max = max(drawing_element.get_y1(), drawing_element.get_y2()) # + drawing_element.get_width()
	elif isinstance(drawing_element, Swoop.Hole):
		x_min = drawing_element.get_x() - (drawing_element.get_drill()/2.0)
		x_max = drawing_element.get_x() + (drawing_element.get_drill()/2.0)
		y_min = drawing_element.get_y() - (drawing_element.get_drill()/2.0)
		y_max = drawing_element.get_y() + (drawing_element.get_drill()/2.0)
	elif isinstance(drawing_element, Swoop.Circle):
		x_min = drawing_element.get_x() - (drawing_element.get_radius()/2.0) - drawing_element.get_width()
		x_max = drawing_element.get_x() + (drawing_element.get_radius()/2.0) + drawing_element.get_width()
		y_min = drawing_element.get_y() - (drawing_element.get_radius()/2.0) - drawing_element.get_width()
		y_max = drawing_element.get_y() + (drawing_element.get_radius()/2.0) + drawing_element.get_width()
	elif isinstance(drawing_element, Swoop.Polygon):
		x_min = min(map(lambda x: x.get_x(), drawing_element.get_vertices())) - drawing_element.get_width()
		x_max = max(map(lambda x: x.get_x(), drawing_element.get_vertices())) + drawing_element.get_width()
		y_min = min(map(lambda x: x.get_y(), drawing_element.get_vertices())) - drawing_element.get_width()
		y_max = max(map(lambda x: x.get_y(), drawing_element.get_vertices())) + drawing_element.get_width()
	elif isinstance(drawing_element, Swoop.Smd):
		rotation = drawing_element.get_rot()
		if not rotation:
			x_min = drawing_element.get_x() - drawing_element.get_dx() # maybe this should be get_dx() / 2.0 etc... I'm unclear on the spec
			x_max = drawing_element.get_x() + drawing_element.get_dx()
			y_min = drawing_element.get_y() - drawing_element.get_dy()
			y_max = drawing_element.get_y() + drawing_element.get_dy()
		elif rotation=='R90':
			x_min = drawing_element.get_x() - drawing_element.get_dy()
			x_max = drawing_element.get_x() + drawing_element.get_dy()
			y_min = drawing_element.get_y() - drawing_element.get_dx()
			y_max = drawing_element.get_y() + drawing_element.get_dx()
		elif rotation=='R180':
			x_min = drawing_element.get_x() - drawing_element.get_dx()
			x_max = drawing_element.get_x() + drawing_element.get_dx()
			y_min = drawing_element.get_y() - drawing_element.get_dy()
			y_max = drawing_element.get_y() + drawing_element.get_dy()
		elif rotation=='R270':
			x_min = drawing_element.get_x() - drawing_element.get_dy()
			x_max = drawing_element.get_x() + drawing_element.get_dy()
			y_min = drawing_element.get_y() - drawing_element.get_dx()
			y_max = drawing_element.get_y() + drawing_element.get_dx()

	elif isinstance(drawing_element, Swoop.Pad):
		rotation = drawing_element.get_rot()
		# this extra is wrong, it needs to account for the pad shape to be true-to-spec
		# this is just a shortcut because I don't have time to reverse engineer it
		extra = max(drawing_element.get_drill(), drawing_element.get_diameter()) / 2.0
		if not rotation:
			x_min = drawing_element.get_x() - extra
			x_max = drawing_element.get_x() + extra
			y_min = drawing_element.get_y() - extra
			y_max = drawing_element.get_y() + extra
		elif rotation: # do nothing for rotation for now. It only affect "long" pads
			x_min = drawing_element.get_x() - extra
			x_max = drawing_element.get_x() + extra
			y_min = drawing_element.get_y() - extra
			y_max = drawing_element.get_y() + extra

	# assert x_min < 9e98, str(drawing_element) + str(((x_min, x_max), (y_min, y_max)) )
	# assert x_max > -9e98, str(dr/awing_element) + str(((x_min, x_max), (y_min, y_max)) )
	# assert y_min < 9e98, str(drawing_element) + str(((x_min, x_max), (y_min, y_max)) )
	# assert y_max > -9e98, str(drawing_element) + str(((x_min, x_max), (y_min, y_max)) )

	return ((x_min, x_max), (y_min, y_max))





def run_conversion(
	user_id = 'No user ID set',
	project_name = '.',
	brd_file = 'unplaced.brd'
):

	brd = Swoop.EagleFile.from_file(brd_file)

	# get the elements (components/blocks/nodes)
	elements = {}
	for n in (Swoop.From(brd).
		get_elements()
	):
		name = n.get_name()			
		library = n.get_library()
		package = n.get_package()
		e = ElementEntry(name, library=library, package=package)
		e.x_loc = n.get_x()
		e.y_loc = n.get_y()
		e.rotation = n.get_rot()
		e.locked = n.get_locked()
		elements[name] = e

	print('total: ' + str(len(elements)) + ' elements (components/blocks/nodes)')

	# get the bounding box for the elements (from lib?)
	for n, e in elements.iteritems():
		eagle_package = Swoop.From(brd).get_library(e.library).get_package(e.package)
		pads = eagle_package.get_pads()
		smds = eagle_package.get_smds()

		for pin in pads+smds:
			e.add_pin(pin)

		drawing = eagle_package.get_drawing_elements()

		for de in drawing + pads + smds:
			allowed_types = [
				Swoop.Wire,
				Swoop.Rectangle,
				Swoop.Hole,
				Swoop.Circle,
				Swoop.Polygon,
				Swoop.Smd,
				Swoop.Pad,
			]
			allowed_layers = set(
				(None,1,16,17,18,29,30,31,32,33,34,35,36,151,39,40,41,42,44,45)
			)

			# if hasattr(de, 'get_layer') and de.get_layer() not in allowed_layers:
				# continue

			
			allowed = [isinstance(de, t) for t in allowed_types]
			if not any(allowed):
				continue

			((x_min, x_max), (y_min, y_max)) = de_bounding_box(de)
			
			e.expand_bb(
				x_min,
				x_max,
				y_min,
				y_max,
			)

			if e.name == 'K1':
				print(de)
				print((e.x_min, e.x_max), (e.y_min, e.y_max))

			# if n == 'JP9':
			# 	print de
			# 	print((x_min, x_max), (y_min, y_max))
			# 	print((e.x_min, e.x_max), (e.y_min, e.y_max))
			# 	print


	# header for nodes file
	nodes_header = ''

	nodes_header += 'UCLA nodes 1.0\n'
	nodes_header += '\n'
	nodes_header += '# Created    : ' + str(datetime.datetime.now()) + '\n'
	nodes_header += '# Created by : ' + user_id + '\n'
	nodes_header += '\n'
	nodes_header += 'NumNodes : ' + str(len(elements)) + '\n'
	nodes_header += 'NumTerminals : ' + '0' + '\n'
	nodes_header += '\n'

	nodes_str = ''
	for n, e in elements.iteritems():
		nodes_str += e.node_str() + '\n'
		if n == 'K1':
			print(e.node_str())

	# get the nets (signals/wires)
	signals = {}
	for n in (Swoop.From(brd).
		get_signals()
	):
		name = n.get_name()
		signal =  Signal(name)
		signals[name] = signal
		c_refs = n.get_contactrefs()
		for c_ref in c_refs:
			assert c_ref.get_element() in elements
			element = elements[c_ref.get_element()]
			# if element.name == 'JP9':
				# print element
				# print element.package
				# pin_names = sorted(element.pins.keys())#, key=lambda x: if isinstance(x, number): int(x) else: x)

				# last_pin = (0,0)
				# for pn in pin_names:
					# print pn, element.pins[pn]
					# print(last_pin[0] - element.pins[pn][0], last_pin[1] - element.pins[pn][1])
					# last_pin = element.pins[pn]
				# exit(-1)
			pin_name = c_ref.get_pad()
			signal.add_pin_absolute(element=element, pin_name=pin_name)


	print('Total: ' + str(len(signals)) + ' nets')

	nets_header = ''
	nets_header += 'UCLA nets 1.0\n'
	nets_header += '\n'
	nets_header += '# Created    : ' + str(datetime.datetime.now()) + '\n'
	nets_header += '# Created by : ' + user_id + '\n'
	nets_header += '\n'
	nets_header += 'NumNets : ' + str(len(signals)) + '\n'
	nets_header += 'NumPins : ' + str(sum([len(s.pins) for n, s in signals.iteritems()])) + '\n'
	nets_header += '\n'

	nets_str = ''
	for n, s in signals.iteritems():
		degree = len(s.pins)
		nets_str += 'NetDegree : ' + str(degree) + '\n'
		for p in s.pins:
			nets_str += p.name.rjust(15) + ' ' + p.direction.rjust(3) + ' : ' + str(p.x_offset).rjust(10) + ' ' + str(p.y_offset).rjust(10) + '\n'
		# nets_str += '\n'

	weights_header = ''
	weights_header += 'UCLA wts 1.0\n'
	weights_header += '\n'
	weights_header += '# Created    : ' + str(datetime.datetime.now()) + '\n'
	weights_header += '# Created by : ' + user_id + '\n'
	weights_header += '\n'

	weights_str = ''

	for i, (n, s) in enumerate(signals.iteritems()):
		name = n
		weight = s.weight
		weights_str += name + ' ' + str(weight) + '\n'


	pl_header = ''
	pl_header += 'UCLA pl 1.0\n'
	pl_header += '\n'
	pl_header += '# Created    : ' + str(datetime.datetime.now()) + '\n'
	pl_header += '# Created by : ' + user_id + '\n'
	pl_header += '\n'

	pl_str = ''
	for n, e in elements.iteritems():

		# if e.name == 'X4':
			# pl_str += e.name.rjust(15) + ' ' + str(-(e.x_max - e.x_min)/2).rjust(10) + ' ' + str(-(e.y_max-e.y_min)/2).rjust(10)
		# else:

		ll_x = e.x_loc # default to origin
		ll_y = e.y_loc # default to origin

		if n == 'K1':
			print(n, e.rotation)
		if (e.rotation is None) or (e.rotation == 'R0'): # N
			# ll_x = e.x_loc - abs(e.x_min)
			# ll_y = e.y_loc - abs(e.y_min)

			ll_x = e.x_loc + (e.x_min)
			ll_y = e.y_loc + (e.y_min)
			# print e.name, 'loc:', (e.x_loc, e.y_loc),'ll:', (ll_x, ll_y), 'x:', (e.x_min, e.x_max), 'y:', (e.y_min, e.y_max)
		elif e.rotation == 'R90':
			ll_x = e.x_loc - (e.y_max)
			ll_y = e.y_loc + (e.x_min)
			# print e.name, 'loc:', (e.x_loc, e.y_loc),'ll:', (ll_x, ll_y), 'x:', (e.x_min, e.x_max), 'y:', (e.y_min, e.y_max)
		elif e.rotation == 'R180':
			ll_x = e.x_loc - (e.x_max)
			ll_y = e.y_loc - (e.y_max)
			# print e.name, 'loc:', (e.x_loc, e.y_loc),'ll:', (ll_x, ll_y), 'x:', (e.x_min, e.x_max), 'y:', (e.y_min, e.y_max)
		elif e.rotation == 'R270':
			ll_x = e.x_loc + (e.y_min)
			ll_y = e.y_loc - (e.x_max)
			# print e.name, 'loc:', (e.x_loc, e.y_loc),'ll:', (ll_x, ll_y), 'x:', (e.x_min, e.x_max), 'y:', (e.y_min, e.y_max)
		# else: # this is wrong, but we don't handle other rotations yet
		# 	pass

		# pl_str += e.name.rjust(15) + ' ' + str(e.x_loc).rjust(10) + ' ' + str(e.y_loc).rjust(10)
		pl_str += e.name.rjust(15) + ' ' + str(ll_x).rjust(10) + ' ' + str(ll_y).rjust(10) # use ll
		if e.name == 'K1':
			print  (e.name.rjust(15) + ' ' + str(ll_x).rjust(10) + ' ' + str(ll_y).rjust(10))

		if e.rotation is None:
			pl_str += ' : N'
		elif e.rotation == 'R90':
			pl_str += ' : W' # really, EAGLE does left hand rotation for some reason
		elif e.rotation == 'R180':
			pl_str += ' : S'
		elif e.rotation == 'R270':
			pl_str += ' : E'
		else:
			pl_str += ' : N'

		if e.locked == True:
			pl_str += ' /FIXED\n'.rjust(12)
		else:
			pl_str += '\n'


	with open(project_name + '.nodes', 'w') as file:
		file.write(nodes_header + nodes_str)

	with open(project_name + '.nets', 'w') as file:
		file.write(nets_header + nets_str)

	with open(project_name + '.wts', 'w') as file:
		file.write(weights_header + weights_str)

	with open(project_name + '.pl', 'w') as file:
		file.write(pl_header + pl_str)


if __name__ == '__main__':
	arguments = docopt(__doc__, version='eagle2bookshelf v0.1')
	run_conversion(
		user_id=str(arguments['--userid']),
		project_name=str(arguments['--output_prfx']),
		brd_file=str(arguments['--brd'])
	)
