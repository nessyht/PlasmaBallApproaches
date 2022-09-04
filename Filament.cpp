#include "Filament.h"
#include <random>
#include <algorithm>
#include <ctime>

Filament::Filament(Eigen::Vector2d root, double charge) : root_(root), new_point_(root), alive_(true), charge_(charge)
{
	points_.push_back(root);
}

Filament::Filament() : alive_(true)
{

}

Filament::~Filament()
{
}

/*
Iterate through all the nodes around the new point
Any valid (not a point of the Filament or already a node and inside the map) points are added to the nodes
*/
void Filament::find_nodes(MatrixXd map) // To Do: Parellise
{
	int height = map.cols();
	int width = map.rows();
	for (Vector2d point : points_)
	{
		std::vector<Vector2d> nodes = get_nodes(point, 3);
		for (auto node : nodes)
		{
			if (node(0) >= 0 && node(1) >= 0 && node(0) < height && node(1) < width)
			{
				if (map(node(0), node(1)) != -1) // inside the map, not the current point
				{
					if (!(std::find(points_.begin(), points_.end(), node) != points_.end()))  // a point
					{
						if (!(std::find(nodes_.begin(), nodes_.end(), node) != nodes_.end())) // already a node
						{
							nodes_.push_back(node);
							/* //debug
							mask(i, j) = 7;
							for (int k = 0; k < height; k++)
							{
								for (int l = 0; l < width; l++)
									std::cout << mask(k, l) << " ";
								std::cout << std::endl;
							}
							*/
						}
					}
				}
			}
		}
		
	}
	
	// sort(nodes_.begin(), nodes_.end());
	// nodes_.erase(unique(nodes_.begin(), nodes_.end()), nodes_.end()); //ensure no duplicates exist
}

std::vector<Vector2d> Filament::get_nodes(Vector2d point, int window_size)
{
	std::vector<Vector2d> nodes;
	for (int i = point(0) - int(window_size / 2); i < point(0) + int(window_size / 2) + 1; ++i) //-1 in y to +1 in y
	{
		for (int j = point(1) - int(window_size / 2); j < point(1) + int(window_size / 2) + 1; ++j) // -1 in x to +1 in x
		{
			nodes.push_back(Vector2d(i, j));
		}
	}
	return nodes;
}


/*
Using the phis_ calculate next point to be selected
Move phis into an eigen array, take each to the power of dim
Divide by the Sum of all phis to the dim
Multiply with random[0,1) and select max as new_point_ from the index in the nodes
*/
void Filament::find_next(double sum, double dim)
{
	srand(time(NULL));

	ArrayXd phis(phis_.size());
	for (int i = 0; i < phis_.size(); i++)
		phis(i) = (phis_[i] + 1)/2; // normalise between 0 <= phi < 1
	//std::cout << std::endl << phis.transpose() << std::endl << sum << std::endl; // debug
	// ArrayXd phis = ArrayXd(phis_).pow(dim) * (1/sum); // /sum
	sum = phis.sum();
	phis = phis.pow(dim) / sum;

	ArrayXd r = ArrayXd::Random(phis_.size()).abs();
	ArrayXd::Index max;
	
	double top = (phis * r).maxCoeff(&max);
	new_point_ = nodes_[max];
	/* debug
	std::cout << "phis: " << phis.transpose() << std::endl
		<< "r: " << r.transpose() << std::endl
		<< "phi*r: " << (phis * r).transpose() << std::endl
		<< "max_I: "<< std::endl << max << std::endl;
	

	std::cout << "new_point: "<< new_point_(0) << ", " << new_point_(1) << std::endl;
	// */
}

/*
add the new_point_ to the points_ of the Filament
*/
void Filament::extend()
{
	points_.push_back(new_point_);
}

