#pragma once
#define CL_HPP_ENABLE_EXCEPTIONS
#define CL_HPP_TARGET_OPENCL_VERSION 200
#define CL_USE_DEPRECATED_OPENCL_2_0_APIS
#include <algorithm>
#include <random>
#include "cubic_spline.h"
#include "super.h"
#include <CL/cl.hpp>
#include <unordered_map>

class Ball
{
public:
	Field e_field_;
	double voltage_;
	double threshold_;
	const int radius_;
	const int pixels_;
	const int width_;
	const int height_;
	const int depth_;
	const double dim_;
	const int time_;
	const int gauss_;
	double min_charge_;
	double max_charge_;
	const double pi_ = 3.14159265358979323846;
	std::unordered_map<std::string, int> show_counters_;
	bool vtk_;
	std::vector<Filament> Filaments_;
	std::vector<Filament> fin_Filaments_;
	int iteration_;
	int leader_iter_;
	std::string dir_name_;
	bool is_ball_;
	bool alive_;
	bool debug_;
	bool GPU_;
	bool laplace_ = true;
	std::vector<Eigen::Vector2d> roots_;
	std::vector<Cubic_Spline> splines_;

	cl::Platform platform_;
	cl::Device device_;
	cl::Context context_;
	cl::Program program_;

	Ball(bool ball, std::string dirname, double dim, bool write, int radius, int pixels, int gauss = 1, bool debug = false, bool run_on_gpu = false, int time = 1, int depth = 0);
	~Ball();
	
	void make_plasma();
	void init_plasma();
	void find_roots(int n);
	void find_roots();
	void find_square_roots(int n);
	void find_round_roots(int n);
	void find_round_roots();
	bool test_outside_ring(Vector2d point, double r);
	void sim_Filaments(double w = 1);
	void find_next(MatrixXd map);
	void make_next();
	double calc_sum();
	void run_sim(int roots = 1);
	void run_field(int its, bool reset, double w = 1);
	void update_alive();
	void test_filament_join(Filament& fila);
	void update_dead_Filament(int i);
	void update_field();
	void create_pos_residual(double w = 1);
	void update_gaussian(double w = 1);
	void create_splines();
	std::vector<std::vector<char>> make_show_map();
	void show(MatrixXd &Map, std::string name = "");
	void show_raw(MatrixXd &Mat);
	void write_vtk_(std::string filename, std::string name, MatrixXd &Map);
	void write_txt_(std::vector<std::vector<char>> display, std::string name);
};
