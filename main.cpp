#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <string.h>


#include "Ball.h"

#include<iostream>
#include<vector>
#include<algorithm>
#include<cmath>

int main() 
{

	int radius = 126; // in mm
	int pixels = 100;
	// int fork_width = 16;
	double dim = 2;
	int gauss = 1;
	int filaments = 6;
	bool run_on_gpu = true;
	std::string foldername = std::to_string(radius) + "_h_with_" + std::to_string(filaments) + "_dim_" + std::to_string(dim);


	Ball globe_0(true, foldername, dim, true, radius, pixels, gauss, false, run_on_gpu, 30);
	globe_0.run_sim(50 + 4 * filaments);

	/*
	Ball globe_2(true, foldername, dim + 0.4, true, radius, pixels, gauss, true, run_on_gpu, 10);
	globe_2.run_sim(50 + 4 * filaments);
	Ball globe_3(true, foldername, dim + 0.6, true, radius, pixels, gauss, true, run_on_gpu, 10);
	globe_3.run_sim(50 + 4 * filaments);
	Ball globe_4(true, foldername, dim + 0.8, true, radius, pixels, gauss, true, run_on_gpu, 10);
	globe_4.run_sim(50 + 4 * filaments);
	Ball globe_5(true, foldername, dim + 1	, true, radius, pixels, gauss, true, run_on_gpu, 10);
	globe_5.run_sim(50 + 4 * filaments);

	filaments += 4;	
	Ball globe_10(true, foldername, dim			, true, radius, pixels, gauss, true, run_on_gpu, 10);
	globe_10.run_sim(50 + 4 * filaments);
	Ball globe_11(true, foldername, dim + 0.2, true, radius, pixels, gauss, true, run_on_gpu, 10);
	globe_11.run_sim(50 + 4);
	Ball globe_12(true, foldername, dim + 0.4, true, radius, pixels, gauss, true, run_on_gpu, 10);
	globe_12.run_sim(50 + 4 * filaments);
	Ball globe_13(true, foldername, dim + 0.6, true, radius, pixels, gauss, true, run_on_gpu, 10);
	globe_13.run_sim(50 + 4 * filaments);
	Ball globe_14(true, foldername, dim + 0.8, true, radius, pixels, gauss, true, run_on_gpu, 10);
	globe_14.run_sim(50 + 4 * filaments);
	Ball globe_15(true, foldername, dim + 1	 , true, radius, pixels, gauss, true, run_on_gpu, 10);
	globe_15.run_sim(50 + 4 * filaments);

	//Ball globe_test(true, foldername, 4, true, radius, radius, gauss, true, true);
	/*
	globe_test.e_field_.add_charge(Vector2d(12, 15), -1);
	globe_test.e_field_.add_charge(Vector2d(11, 15), -1);
	globe_test.e_field_.add_charge(Vector2d(10, 15), -1);
	globe_test.e_field_.add_charge(Vector2d(9, 15), -1);
	globe_test.e_field_.add_charge(Vector2d(8, 15), -1);
	globe_test.e_field_.pos_res_map_(10, 14) = 100;
	//globe_test.e_field_.pos_res_map_.block(22, 3, 1, 20) = MatrixXd::Ones(1, 20) * 0.5; 
	//globe_test.e_field_.pos_res_map_.block(27, 3, 1, 20) = MatrixXd::Ones(1, 20) * 0.5;
	//globe_test.e_field_.pos_res_map_.block(2, 25, 20, 1) = MatrixXd::Ones(20, 1) * 0.5;
	//globe_test.e_field_.gaussian_blur(globe_test.e_field_.pos_res_map_);
	globe_test.e_field_.rho_map_ = MatrixXd::Ones(radius * 2 + 2, radius * 2 + 2);

	//globe_test.show(globe_test.e_field_.base_gauss_map_, "test_base_gauss_map");
	//globe_test.show(globe_test.e_field_.pos_res_map_, "test_pos_res_map");
	//globe_test.e_field_.rho_map_ = globe_test.e_field_.pos_res_map_ - globe_test.e_field_.base_gauss_map_;
	globe_test.show(globe_test.e_field_.rho_map_, "test_rho_map");
	globe_test.show(globe_test.e_field_.phi_map_, "test_init_phi_map");

	// globe_test.e_field_.phi_map_.block(0, 0, 2, 50) = MatrixXd::Ones(2, 50) * 10;
	globe_test.e_field_.iter_laplace(globe_test.platform_, globe_test.device_, globe_test.context_, 10);


	globe_test.show(globe_test.e_field_.phi_map_, "test_phi_map");
	// */




	return 0;
}
