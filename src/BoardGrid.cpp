// BoardGrid.cpp
#include "BoardGrid.h"

bool Location::operator==(const Location &other) const
{
	return (this->x == other.x) && (this->y == other.y) && (this->z == other.z);
}

std::ostream &operator<<(std::ostream &os, Location const &l)
{
	return os << "Location(" << l.x << ", " << l.y << ", " << l.z << ")";
}

void BoardGrid::initilization(int w, int h, int l)
{
	this->w = w;
	this->h = h;
	this->l = l;
	this->size = w * h * l;

	assert(this->base_cost == nullptr);
	assert(this->working_cost == nullptr);
	assert(this->via_cost == nullptr);

	this->base_cost = new float[this->size];
	this->working_cost = new float[this->size];
	this->via_cost = new float[this->size];

	assert(this->base_cost != nullptr);
	assert(this->working_cost != nullptr);
	assert(this->via_cost != nullptr);
}

void BoardGrid::base_cost_fill(float value)
{
	for (int i = 0; i < this->size; i++)
	{
		this->base_cost[i] = value;
	}
}

void BoardGrid::working_cost_fill(float value)
{
	for (int i = 0; i < this->size; i++)
	{
		this->working_cost[i] = value;
	}
}

float BoardGrid::cost_to_occupy(const Location &l) const
{
	return this->base_cost_at(l) + this->via_cost_at(l);
}

float BoardGrid::base_cost_at(const Location &l) const
{
#ifdef BOUND_CHECKS
	assert((l.x + l.y * this->w + l.z * this->w * this->h) < this->size);
#endif
	return this->base_cost[l.x + l.y * this->w + l.z * this->w * this->h];
}

float BoardGrid::via_cost_at(const Location &l) const
{
	// std::cerr << "via_cost_at( " << l << ")" << std::endl;
	// std::cerr << this->via_cost << std::endl;

	// if (!this->via_cost) {
	// 	std::cout << "Could not allocate via_cost" << std::endl;
	//     	exit(-1);
	// }

#ifdef BOUND_CHECKS
	assert((l.x + l.y * this->w + l.z * this->w * this->h) < this->size);
#endif
	return this->via_cost[l.x + l.y * this->w + l.z * this->w * this->h];
	// std::cerr << "returning via_cost_at" << std::endl;
}

float BoardGrid::working_cost_at(const Location &l) const
{
#ifdef BOUND_CHECKS
	assert((l.x + l.y * this->w + l.z * this->w * this->h) < this->size);
#endif
	return this->working_cost[l.x + l.y * this->w + l.z * this->w * this->h];
}

void BoardGrid::base_cost_set(float value, const Location &l)
{
#ifdef BOUND_CHECKS
	assert(l.x + l.y * this->w + l.z * this->w * this->h < this->size);
#endif
	this->base_cost[l.x + l.y * this->w + l.z * this->w * this->h] = value;
}

void BoardGrid::base_cost_add(float value, const Location &l)
{
#ifdef BOUND_CHECKS
	assert(l.x + l.y * this->w + l.z * this->w * this->h < this->size);
#endif
	this->base_cost[l.x + l.y * this->w + l.z * this->w * this->h] += value;
}

void BoardGrid::working_cost_set(float value, const Location &l)
{
#ifdef BOUND_CHECKS
	assert(l.x + l.y * this->w + l.z * this->w * this->h < this->size);
#endif
	this->working_cost[l.x + l.y * this->w + l.z * this->w * this->h] = value;
}

void BoardGrid::via_cost_set(float value, const Location &l)
{
#ifdef BOUND_CHECKS
	assert(l.x + l.y * this->w + l.z * this->w * this->h < this->size);
#endif
	this->via_cost[l.x + l.y * this->w + l.z * this->w * this->h] = value;
}

void BoardGrid::via_cost_add(float value, const Location &l)
{
#ifdef BOUND_CHECKS
	assert(l.x + l.y * this->w + l.z * this->w * this->h < this->size);
#endif
	this->via_cost[l.x + l.y * this->w + l.z * this->w * this->h] += value;
}

std::unordered_map<Location, Location> BoardGrid::dijkstras_with_came_from(
	const std::vector<Location> &route,
	int via_size)
{
	std::cout << "Starting dijkstras_with_came_from ==Multipin== nets: route.features.size() = " << route.size() << std::endl;

	// For path to multiple points
	// Searches from the multiple points to every other point
	this->working_cost_fill(std::numeric_limits<float>::infinity());

	LocationQueue<Location, float> frontier;		  // search frontier
	std::unordered_map<Location, Location> came_from; // cheapest neighbor
	for (Location start : route)
	{
		this->working_cost_set(0.0, start);
		frontier.push(start, 0.0);
		came_from[start] = start;
	}

	std::cout << "came_from.size() = " << came_from.size() << ", frontier.size(): " << frontier.size() << std::endl;

	while (!frontier.empty())
	{
		Location current = frontier.front();
		frontier.pop();

		//std::cout << "Visiting " << current << ", frontierSize: "<< frontier.size() << std::endl;

		for (std::pair<float, Location> next : this->neighbors(current, via_size))
		{
			if (
				(next.second.x < 0) || (next.second.x >= this->w) ||
				(next.second.y < 0) || (next.second.y >= this->h) ||
				(next.second.z < 0) || (next.second.z >= this->l))
			{
				continue; // continue if out of bounds
			}
			// std::cerr << "next.second.x: " << next.second.x << std::endl;
			// std::cerr << "next.second.y: " << next.second.y << std::endl;
			// std::cerr << "next.second.z: " << next.second.z << std::endl;

			// std::cerr << "geting new cost" << std::endl;

			// this->via_cost_at(next.second) ??????????
			float new_cost = this->working_cost_at(current) + this->base_cost_at(next.second) + this->via_cost_at(next.second) + next.first;

			// std::cerr << "Done" << std::endl;

			if (new_cost < this->working_cost_at(next.second))
			{
				// std::cerr << "setting working cost" << std::endl;
				this->working_cost_set(new_cost, next.second);
				came_from[next.second] = current;
				frontier.push(next.second, new_cost);
				// std::cerr << "Done" << std::endl;
			}

			// std::cerr << std::endl;
		}
	}

	// std::cerr << "finished dijkstras_with_came_from" << std::endl;

	return came_from;
}

std::array<std::pair<float, Location>, 10> BoardGrid::neighbors(const Location &l, int via_size) const
{
	std::array<std::pair<float, Location>, 10> ns;

	return ns;
}

void BoardGrid::printGnuPlot()
{
	float max_val = 0.0;
	for (int i = 0; i < this->size; i += 1)
	{
		if (this->base_cost[i] > max_val)
			max_val = this->base_cost[i];
	}

	std::cout << "printGnuPlot()::Max Cost: " << max_val << std::endl;

	for (int l = 0; l < this->l; ++l)
	{
		std::string outFileName = "layer" + std::to_string(l) + "_baseCost.dat";
		outFileName = util::appendDirectory(GlobalParam::gOutputFolder, outFileName);
		std::ofstream ofs(outFileName, std::ofstream::out);
		ofs << std::fixed << std::setprecision(5);

		for (int r = 0; r < this->h; ++r)
		{
			for (int c = 0; c < this->w; ++c)
			{
				ofs << this->base_cost_at(Location(c, r, l)) << " ";
			}
			ofs << std::endl;
		}
	}
}

void BoardGrid::printMatPlot()
{
	float maxCost = std::numeric_limits<float>::min();
	float minCost = std::numeric_limits<float>::max();
	for (int i = 0; i < this->size; i += 1)
	{
		if (this->base_cost[i] > maxCost)
		{
			maxCost = this->base_cost[i];
		}
		else if (this->base_cost[i] < minCost)
		{
			minCost = this->base_cost[i];
		}
	}

	std::cout << "printGnuPlot()::Max Cost: " << maxCost << ", Min Cost: " << minCost << std::endl;

	for (int l = 0; l < this->l; ++l)
	{
		std::string outFileName = "layer" + std::to_string(l) + "_baseCost.py";
		outFileName = util::appendDirectory(GlobalParam::gOutputFolder, outFileName);
		std::ofstream ofs(outFileName, std::ofstream::out);
		std::cout << "outFileName: " << outFileName << std::endl;

		ofs << std::fixed << std::setprecision(5);
		ofs << "import numpy as np\n";
		ofs << "import matplotlib.pyplot as plt\n";
		ofs << "plt.close()\n";
		ofs << "viridis = plt.get_cmap('viridis', 12)\n";
		ofs << "data = np.array([[";

		for (int r = 0; r < this->h; ++r)
		{
			for (int c = 0; c < this->w; ++c)
			{
				ofs << this->base_cost_at(Location(c, r, l)) << " ";
				if (c < this->w - 1)
				{
					ofs << ", ";
				}
				else
				{
					ofs << "]";
				}
			}

			if (r < this->h - 1)
			{
				ofs << ", [";
			}
		}

		ofs << "])\n";
		ofs << "plt.pcolormesh(data, cmap=viridis, vmin=data.min(), vmax=data.max())\n";
		ofs << "plt.title('test123')\n";
		ofs << "plt.colorbar()\n";
		ofs << "plt.show()\n";
	}
}

void BoardGrid::pprint()
{
	char colors[11] = " .,:=+*#%@";
	// long
	// colors[36]    = {' ','.',':','░','#','▒','▓','█'};
	int num_colors = 10;

	float max_val = 0.0;
	for (int i = 0; i < this->size; i += 1)
	{
		if (this->base_cost[i] > max_val)
			max_val = this->base_cost[i];
	}

	for (int l = 0; l < this->l; l += 1)
	{
		std::cout << std::endl
				  << "Layer: " << l << std::endl;

		for (int r = 0; r < this->h; r += 1)
		{
			for (int c = 0; c < this->w; c += 1)
			{
				int current_color = this->base_cost_at(Location(c, r, l)) / max_val * num_colors;
				std::cout << colors[current_color] << " ";
				// std::cout << current_color << " ";

				// std::cout << std::setw(6) << std::setprecision(2) << std::right << this->at(c, r, l);
			}
			std::cout << std::endl;
		}
	}
}

float BoardGrid::sized_via_cost_at(const Location &l, int via_size) const
{
	int radius = via_size;
	float cost = 0.0;
	for (int z = 0; z < this->l; z += 1)
	{
		for (int y = -radius; y < radius; y += 1)
		{
			for (int x = -radius; x < radius; x += 1)
			{
				Location current_l = Location(l.x+x, l.y+y, z);
				if (!validate_location(current_l))
				{
					// std::cerr << 'Invalid location: ' << current_l << std::endl;
					cost += 1000;
					continue;
				}
				cost += this->via_cost_at(current_l); // + this->via_cost_at(l);
				cost += this->base_cost_at(current_l); // + this->via_cost_at(l);
			}
		}
	}
	// std::cerr << "sized_via_cost_at " << l << ": " << cost << std::endl;

	return cost;
}

float BoardGrid::sized_trace_cost_at(const Location &l, int traceRadius) const
{
	int radius = traceRadius;
	float cost = 0.0;
	for (int y = -radius; y < radius; y += 1)
	{
		for (int x = -radius; x < radius; x += 1)
		{
			Location current_l = Location(l.x + x, l.y + y, l.z);
			if (!validate_location(current_l))
			{
				//TODO
				cost += 100000;
				continue;
			}
			cost += this->base_cost_at(current_l);
		}
	}
	return cost;
}

void BoardGrid::print_came_from(const std::unordered_map<Location, Location> &came_from, const Location &end)
{
	// Location end(ex, ey, ez);

	std::cout << "Came froms: " << std::endl;
	std::cout << "\tsize: " << came_from.size() << std::endl;

	for (int l = 0; l < this->l; l += 1)
	{
		std::cout << std::endl
				  << "Layer: " << l << std::endl;

		for (int r = 0; r < this->h; r += 1)
		{
			for (int c = 0; c < this->w; c += 1)
			{

				Location current(c, r, l);

				if (current == end)
				{
					std::cout << "# ";
					continue; // goal
				}

				if (came_from.find(current) == came_from.end())
				{
					std::cout << ". "; // not found
					continue;		   // unexplored
				}

				Location cf = came_from.find(current)->second;

				if (current.x > cf.x)
				{
					std::cout << "> ";
				}
				else if (current.x < cf.x)
				{
					std::cout << "< ";
				}
				else if (current.y > cf.y)
				{
					std::cout << "V ";
				}
				else if (current.y < cf.y)
				{
					std::cout << "^ ";
				}
				else if (current.z > cf.z)
				{
					std::cout << "X ";
				}
				else if (current.z < cf.z)
				{
					std::cout << "O ";
				}
				else if ((current.x == cf.x) && (current.y == cf.y) && (current.z == cf.z))
				{
					std::cout << "* "; // start
				}
			}
			std::cout << std::endl;
		}
	}
}

void BoardGrid::add_via_cost(const Location &l, int layer)
{

	int radius = 10;
	float cost = 10.0;
	for (int y = -radius; y < radius; y += 1)
	{
		for (int x = -radius; x <= radius; x += 1)
		{
#ifdef BOUND_CHECKS
			assert(((l.x + x) + (l.y + y) * this->w + (layer) * this->w * this->h) < this->size);
#endif
			this->via_cost[(l.x + x) + (l.y + y) * this->w + (layer) * this->w * this->h] += cost;
		}
	}
}

void BoardGrid::remove_via_cost(const Location &l, int layer)
{
	int radius = 10;
	float cost = 10.0;
	for (int y = -radius; y < radius; y += 1)
	{
		for (int x = -radius; x <= radius; x += 1)
		{
#ifdef BOUND_CHECKS
			assert(((l.x + x) + (l.y + y) * this->w + (layer) * this->w * this->h) < this->size);
#endif
			this->via_cost[(l.x + x) + (l.y + y) * this->w + (layer) * this->w * this->h] -= cost;
		}
	}
}

void BoardGrid::add_route_to_base_cost(const MultipinRoute &route, int radius, float cost, int via_size)
{
	// std::vector<Location> features = route.features;
	Location last_location = route.features[0];
	for (Location l : route.features)
	{
		// std::cout << "setting cost for feature " << l << std::endl;
		int layer = l.z;

		// Add costs for vias
		if ((l.z != last_location.z) && (l.x == last_location.x) && (l.y == last_location.y))
		{
			for (int z = 0; z < this->l; z += 1)
			{
				this->add_via_cost(l, z);
			}
		}

		for (int current_radius = 0; current_radius <= radius; current_radius += 1)
		{
			float current_cost = cost;
			if (current_cost <= 0)
				break;

			for (int r = l.y - current_radius; r <= l.y + current_radius; r += 1)
			{
				if (r < 0)
					continue;
				if (r >= this->h)
					continue;
				for (int c = l.x - current_radius; c <= l.x + current_radius; c += 1)
				{
					if (c < 0)
						continue;
					if (c >= this->w)
						continue;
					// std::cout << "\tsetting cost at " << Location(c, r, layer) << std::endl;
					this->base_cost_set(
						this->base_cost_at(Location(c, r, layer)) + current_cost,
						Location(c, r, layer));
				}
			}
		}
	}
}

void BoardGrid::remove_route_from_base_cost(const MultipinRoute &route, int radius, float cost)
{
	// std::vector<Location> features = route.features;
	std::cout << "Starting remove_route_from_base_cost" << std::endl;
	for (Location l : route.features)
	{
		// std::cout << "setting cost for feature " << l << std::endl;
		int layer = l.z;

		if (layer > 2)
		{
			std::cout << "Bad layer: " << l << std::endl;
			exit(-1);
		}

		for (int current_radius = 0; current_radius <= radius; current_radius += 1)
		{
			float current_cost = cost - current_radius;
			if (current_cost <= 0)
				break;

			for (int r = l.y - current_radius; r <= l.y + current_radius; r += 1)
			{
				if (r < 0)
					continue;
				if (r >= this->h)
					continue;
				for (int c = l.x - current_radius; c <= l.x + current_radius; c += 1)
				{
					if (c < 0)
						continue;
					if (c >= this->w)
						continue;
					std::cout << "\tSetting cost at " << Location(c, r, layer) << std::endl;
					this->base_cost_set(this->base_cost_at(Location(c, r, layer)) - current_cost, Location(c, r, layer));
					std::cout << "\tFinised setting cost" << std::endl;
				}
			}
		}
	}
	std::cout << "Finished remove_route_from_base_cost" << std::endl;
}

void BoardGrid::came_from_to_features(
	const std::unordered_map<Location, Location> &came_from,
	const Location &end,
	std::vector<Location> &features) const
{

	// std::vector<Location> features;
	// features.insert(end);
	std::cout << "Starting came_from_to_features " << std::endl;

	if (!this->validate_location(end))
		std::cout << "Bad end for came_from_to_features" << std::endl;

	features.push_back(end);
	Location current = end;

	while (came_from.find(current) != came_from.end())
	{ // while not the start
		Location next = came_from.find(current)->second;
		// if (features.find(next) != features.end()) { // next already in features, loops or start
		// break;
		// }
		if (next == current)
		{
			break;
		}
		// features.insert(next);
		features.push_back(next);
		current = next;
	}

	for (Location l : features)
	{
		//std::cout << "\t" << l << std::endl;
		if (l.z > 2)
		{
			exit(-1);
		}
	}

	std::cout << "Finished came_from_to_features " << std::endl;

	// return features;
}

std::vector<Location> BoardGrid::came_from_to_features(
	const std::unordered_map<Location, Location> &came_from,
	const Location &end) const
{
	std::vector<Location> features;
	// features.insert(end);
	this->came_from_to_features(came_from, end, features);
	return features;
}

void BoardGrid::print_route(const std::unordered_map<Location, Location> &came_from, const Location &end)
{
	std::vector<Location> features;
	features = this->came_from_to_features(came_from, end);


	std::cout << "Printing features" << std::endl;
	for (int l = 0; l < this->l; l += 1)
	{
		std::cout << std::endl
				  << "Layer: " << l << std::endl;

		for (int r = 0; r < this->h; r += 1)
		{
			for (int c = 0; c < this->w; c += 1)
			{

				Location current(c, r, l);

				if (std::find(features.begin(), features.end(), current) != features.end())
				{
					// if (features.find(current) != features.end()) {
					std::cout << "# "; // contains feature
				}
				else
				{
					std::cout << ". ";
				}
			}
			std::cout << std::endl;
		}
	}
}

void BoardGrid::print_features(std::vector<Location> features)
{
	std::cout << "Printing features" << std::endl;
	for (int l = 0; l < this->l; l += 1)
	{
		std::cout << std::endl
				  << "Layer: " << l << std::endl;

		for (int r = 0; r < this->h; r += 1)
		{
			for (int c = 0; c < this->w; c += 1)
			{

				Location current(c, r, l);

				if (std::find(features.begin(), features.end(), current) != features.end())
				{
					// if (features.find(current) != features.end()) {
					std::cout << "# "; // contains feature
				}
				else
				{
					std::cout << ". ";
				}
			}
			std::cout << std::endl;
		}
	}
}

void BoardGrid::add_route(MultipinRoute &route)
{
	int radius = 7;
	int cost = 10;
	int via_size = 7;
	int num_pins = route.pins.size();
	current_trace_width = route.trace_width;
	current_half_trace_width = route.trace_width / 2;
	current_clearance = route.clearance;

	if (num_pins <= 0)
	{
		return;
	}
	else if (num_pins == 1)
	{
		return;
	}
	else
	{
		route.features.push_back(route.pins[0]);

		for (Location pin : route.pins)
		{
			std::unordered_map<Location, Location> came_from = this->dijkstras_with_came_from(route.features, via_size);
			//this->print_came_from(came_from, pin);
			std::vector<Location> new_features = this->came_from_to_features(came_from, pin);

			Location last_feature = new_features[0];
			for (Location f : new_features)
			{
				route.features.push_back(f);
				last_feature = f;
			}
		}
		//this->print_features(route.features);
		//TODO
		//this->add_route_to_base_cost(route, traceWidth, cost, via_size);
		this->add_route_to_base_cost(route, current_trace_width, cost, via_size);
	}
}

void BoardGrid::ripup_route(MultipinRoute &route)
{
	std::cout << "Doing ripup" << std::endl;
	for (Location l : route.features)
	{
		if (l.x > this->w || l.y > this->h || l.z > this->l)
		{
			std::cout << "Bad route to ripup: " << l << std::endl;
			exit(-1);
		}
	}
	this->remove_route_from_base_cost(route, 1, 10);
	std::cout << "Clearing features" << std::endl;

	route.features.clear();
	std::cout << "Finished ripup" << std::endl;
}

bool BoardGrid::validate_location(const Location &l) const
{
	if (l.x >= this->w || l.x < 0)
		return false;
	if (l.y >= this->h || l.y < 0)
		return false;
	if (l.z >= this->l || l.z < 0)
		return false;
	return true;
}
