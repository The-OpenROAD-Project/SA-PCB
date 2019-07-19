"""Eagle2Bookshelf.

This program converts EAGLE board files (.brd) to bookshelf format files.
Specifically this program outputs block component files (.blocks), netlist files (.nets), and net weight files (.wts).
These files are sutible for academic IC placement programs.

This version DOES account for pin placement. But, this feature is unverified.
This version does not account for net weights. All weights are set to '1'.

This program was written by Devon Merrill (devon@ucsd.edu).

Usage:
  eagle2bookshelf.py -h | --help
  eagle2bookshelf.py --brd <BRD> --output_prfx <STEM_NAME> --userid <USERID>

-h --help                      Show this message.
-i --brd BRD                   The EAGLE .brd file to convert.
-o --output_prfx STEM_NAME     The stem name for the new files. Includes directory.
--userid USERID                Your name and contact.
"""

import datetime

import Swoop
from docopt import docopt


class Pin(object):
	def __init__(self, name, x_offset, y_offset, direction):
		super(Pin, self).__init__()
		self.name = name
		self.x_offset = x_offset
		self.y_offset = y_offset
		self.direction = direction
		

class Signal(object):
	def __init__(self, name, weight=1):
		self.name = name
		self.pins = []
		self.weight = weight

	def add_pin(self, element, pin_name, direction='B'):
		origin_offsets = element.pins[pin_name] # the pin dict stores (x,y) pairs of origin offsets
		new_pin = Pin(element.name, '%0.0', '%0.0', direction)
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


class ElementEntry(object):
	def __init__(self, name, library=None, package=None):
		self.name = name
		self.library = library
		self.package = package
		self.x_min = 0.0
		self.x_max = 0.0
		self.y_min = 0.0
		self.y_max = 0.0
		self.pins = {}
		self.bounding_box_multiplier = 1.0

	def __str__(self):
		width = self.x_max - self.x_min
		height = self.y_max - self.y_min
		width_expand = width * self.bounding_box_multiplier / 2.0
		height_expand = height * self.bounding_box_multiplier / 2.0
		box_str = ''
		box_str += '(' + str(self.x_min-width_expand) + ', ' + str(self.y_min-height_expand) + ')' + ', '
		box_str += '(' + str(self.x_min-width_expand) + ', ' + str(self.y_max+height_expand) + ')' + ', '
		box_str += '(' + str(self.x_max+width_expand) + ', ' + str(self.y_max+height_expand) + ')' + ', '
		box_str += '(' + str(self.x_max+width_expand) + ', ' + str(self.y_min-height_expand) + ')'
		return self.name + ' ' + 'hardrectilinear' + ' ' + '4' + ' ' + box_str

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

def run_conversion(
	user_id = 'No user ID set',
	project_name = '.',
	brd_file = 'unplaced.brd'
):

	brd = Swoop.EagleFile.from_file(brd_file)

	# get the elements (components/blocks)



	def de_bounding_box(drawing_element):
		"""Return a bounding box for the drawing element."""
		x_min = 9e99
		x_max = -9e99
		y_min = 9e99
		y_max = -9e99
		if isinstance(drawing_element, Swoop.Wire):
			x_min = min(drawing_element.get_x1(), drawing_element.get_x2()) - drawing_element.get_width()
			x_max = max(drawing_element.get_x1(), drawing_element.get_x2()) + drawing_element.get_width()
			y_min = min(drawing_element.get_y1(), drawing_element.get_y2()) - drawing_element.get_width()
			y_max = max(drawing_element.get_y1(), drawing_element.get_y2()) + drawing_element.get_width()
		elif isinstance(drawing_element, Swoop.Rectangle):
			x_min = min(drawing_element.get_x1(), drawing_element.get_x2()) - drawing_element.get_width()
			x_max = max(drawing_element.get_x1(), drawing_element.get_x2()) + drawing_element.get_width()
			y_min = min(drawing_element.get_y1(), drawing_element.get_y2()) - drawing_element.get_width()
			y_max = max(drawing_element.get_y1(), drawing_element.get_y2()) + drawing_element.get_width()
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
			elif rotation: # just flip x and y I think. This part could be wrong
				y_min = drawing_element.get_x() - drawing_element.get_dx() # maybe this should be get_dx() / 2.0 etc... I'm unclear on the spec
				y_max = drawing_element.get_x() + drawing_element.get_dx()
				x_min = drawing_element.get_y() - drawing_element.get_dy()
				x_max = drawing_element.get_y() + drawing_element.get_dy()
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
			elif rotation: # just flip x and y I think. This part could be wrong
				y_min = drawing_element.get_x() - extra
				y_max = drawing_element.get_x() + extra
				x_min = drawing_element.get_y() - extra
				x_max = drawing_element.get_y() + extra

		assert x_min < 9e98, str(drawing_element) + str(((x_min, x_max), (y_min, y_max)) )
		assert x_max > -9e98, str(drawing_element) + str(((x_min, x_max), (y_min, y_max)) )
		assert y_min < 9e98, str(drawing_element) + str(((x_min, x_max), (y_min, y_max)) )
		assert y_max > -9e98, str(drawing_element) + str(((x_min, x_max), (y_min, y_max)) )

		return ((x_min, x_max), (y_min, y_max))

	elements = {}
	for n in (Swoop.From(brd).
		get_elements()
	):
		name = n.get_name()
		library = n.get_library()
		package = n.get_package()
		elements[name] = ElementEntry(name, library=library, package=package)

	print 'total: ' + str(len(elements)) + ' elements'

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

	blocks_header = ''

	blocks_header += 'UCSC blocks 1.0\n'
	blocks_header += '\n'
	blocks_header += '# Created    : ' + str(datetime.datetime.now()) + '\n'
	blocks_header += '# Created by : ' + user_id + '\n'
	blocks_header += '\n'
	blocks_header += 'NumSoftRectangularBlocks : ' + '0' + '\n'
	blocks_header += 'NumHardRectilinearBlocks : ' + str(len(elements)) + '\n'
	blocks_header += 'NumTerminals : ' + '0' + '\n'
	blocks_header += '\n'

	blocks_str = ''
	for n, e in elements.iteritems():
		blocks_str += str(e) + '\n'

	# get the nets (signals)


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
			pin_name = c_ref.get_pad()
			signal.add_pin(element=element, pin_name=pin_name)

	print 'Total: ' + str(len(signals)) + ' nets'

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
			nets_str += p.name + ' ' + p.direction + ' : ' + p.x_offset + ' ' + p.y_offset + '\n'

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

	with open(project_name + '.blocks', 'w') as file:
		file.write(blocks_header + blocks_str)

	with open(project_name + '.nets', 'w') as file:
		file.write(nets_header + nets_str)

	with open(project_name + '.wts', 'w') as file:
		file.write(weights_header + weights_str)


if __name__ == '__main__':
	arguments = docopt(__doc__, version='eagle2bookshelf v0.1')
	run_conversion(
		user_id=str(arguments['--userid']),
		project_name=str(arguments['--output_prfx']),
		brd_file=str(arguments['--brd'])
	)