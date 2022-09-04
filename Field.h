#pragma once
#define CL_HPP_ENABLE_EXCEPTIONS
#define CL_HPP_TARGET_OPENCL_VERSION 200
#define CL_USE_DEPRECATED_OPENCL_2_0_APIS
#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <array>
#include <CL/cl.hpp>
#include <fstream>
#include <string>
#include "Filament.h"

using namespace Eigen;

class Field
{
private:

	MatrixXd map_;

public:


	MatrixXd phi_map_;
	MatrixXd rho_map_;
	MatrixXd base_gauss_map_;
	MatrixXd pos_res_map_;
	MatrixXd mesh_;
	MatrixXd mask_;
	MatrixXd la_map_;
	MatrixXd blank_map_;
	MatrixXd old_phi_map_;

	const int height_;
	const int width_;
	const int depth_;
	const int charge_ = -1;
	double conv_ = 0.001;
	double inner_r_;
	double outer_r_;

	const int numRows_;
	const int numCols_;
	const int count_;
	bool convergence_ = false;
	bool gpu_is_valid_ = false;

	cl::Platform platform_;
	cl::Device device_;
	cl::Context context_;
	cl::Program laplace_program_;

	//** INITS **//
	Field(bool ball, bool GPU, int gauss = 1, int height = 100, int width = 100, int depth = 0);
	~Field();
	void create_rect_field();
	void create_round_field(double w = 1);
	void reset_mask();
	void create_ring(double r, double fill, MatrixXd &map);
	void create_ring_outer(double r, double fill, MatrixXd &map);
	bool test_outside_ring(Vector2d point, double r);
	void create_mask(MatrixXd &map);

	void update_mask(const std::vector<Filament>& Filaments);
	void add_charge(Vector2d point, double charge);

	//** TESTS + GETTERS **//
	bool test_in_gas(Vector2d point);
	Vector2d get_coords(int point);
	void make_matrix();
	void clamp(MatrixXd& mat);

	//** EQUATIONS **//
	void laplace_threaded(double w);
	void laplace(double w = 1);
	void calc_laplace(double width, int slot = 1, double w = 1);
	void iter_laplace(double w = 1);
	void iter_laplace(cl::Platform &platform, cl::Device &device, cl::Context &context, int its = 100, bool reset = true, double w = 1);

	void calculate_convergence();

	void gaussian_blur(MatrixXd& map);

	//** GPU **//
	cl::Program init_GPU(std::string name);
	void laplace_GPU(double w = 1);




};

