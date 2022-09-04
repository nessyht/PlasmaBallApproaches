#pragma once
#include <Eigen/Dense>
#include <vector>
#include <iostream>

using namespace Eigen;

class Filament
{
public:
	bool parent_ = false;
	Vector2d root_;
	Vector2d new_point_;
	std::vector<Vector2d> points_;
	std::vector<Vector2d> nodes_;
	std::vector<Vector2d> charges_;
	std::vector<double> phis_;
	//Field field_;
	bool alive_;
	double charge_;


	Filament(Eigen::Vector2d root, double charge = 4.);
	Filament();
	~Filament();

	void extend();

	void find_nodes(MatrixXd mask);
	std::vector<Vector2d> get_nodes(Vector2d point, int window_size);

	void find_next(double sum, double dim = 1);
};

