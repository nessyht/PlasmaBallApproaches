#include "Ball.h"
#include <CL/cl.hpp>

/*
#include <fileapi.h>
#include <windows.h>
#include <strsafe.h>
*/

/*
bool ball: True when a ball, false when a Rectangle
std::string dirname: directory name
int dim: dimension
bool write: true when wishing to write to vtk, false otherwise
bool debug: write each iteration in debug mode
int dimensions of grid
*/
Ball::Ball(bool ball, std::string dirname, double dim, bool write_vtk, int radius, int pixels, int gauss, bool debug, bool run_on_gpu, int time, int depth) : is_ball_(ball),
					dir_name_(dirname), radius_(radius), pixels_(pixels), height_(2 * radius), width_(2 * radius), depth_(2 * depth), e_field_(Field(ball, run_on_gpu, gauss, radius, pixels, depth)),
					alive_(true), dim_(dim), vtk_(write_vtk), GPU_(run_on_gpu), debug_(debug), time_(time), iteration_(0), leader_iter_(0), gauss_(gauss)
{
	max_charge_ = e_field_.phi_map_.maxCoeff();
	min_charge_ = e_field_.phi_map_.minCoeff();
	/*
	if (write)
	{
		if (!(CreateDirectoryA(dirname.c_str(), NULL) ||
			ERROR_ALREADY_EXISTS ==  GetLastError()))
		{
			write_vtk_ = false;
		}
	}
	*/
}

void Ball::run_sim(int roots) 
{
	voltage_ = roots;
	find_roots(voltage_); // using these roots execute DBM
	if (debug_)
	{
		show(e_field_.mask_);
		//show(e_field_.base_gauss_map_, "rho_init"); // no base map anymore
	}
	make_plasma();
}

void Ball::find_roots(int n)
{
	if (is_ball_)
		find_round_roots();
	else
		find_square_roots(n);
}

void Ball::find_roots()
{
	find_round_roots();
}

/*
int n: number of roots
roots are evenly spaced along upper edge
*/
void Ball::find_square_roots(int n)
{
	int step = width_ / n;
	for (int x = width_ / (n * 2); x < width_; x += step)
		roots_.push_back(Vector2d(height_, x - 1)); //temp single Filament
}

/*
int n: number of roots
roots are randomly placed on the outer ring of the inner sphere
*/
void Ball::find_round_roots(int n)
{
	double w = width_ / 2.;
	double h = height_ / 2.;
	std::vector<Vector2d> root_pool;
	double r = e_field_.inner_r_;

	int max_y = 0;
	int max_yi = 0;
	for (int y = 0; y < height_; ++y)
	{
		for (int x = 0; x < width_; ++x)
		{
			// if ((sqrt(double(y - h) * double(y - h) + double(x - w) * double(x - w)) < r) && // inside the circle
			//		(sqrt(double(y - h) * double(y - h) + double(x - w) * double(x - w)) >= r - 1)) // outer ring
			if (!test_outside_ring(Vector2d(x, y), r) && 
					test_outside_ring(Vector2d(x,y), r-1))
			{

				root_pool.push_back(Vector2d(y, x));
				/* 
				max_y = (y > max_y) ? y : max_y;
				max_yi = (y >= max_y) ? root_pool.size()-1 : max_yi;
				std::cout << y << ", " << max_y << std::endl;
				*/
			}
		}
	}

	// roots_.push_back(root_pool[max_yi]);
	
	// auto rd = std::random_device{};
	auto rng = std::default_random_engine{  };
	std::shuffle(std::begin(root_pool), std::end(root_pool), rng);
	for (int i = 0; i < n; ++i) // DEBUG < n -1 
	{
		roots_.push_back(root_pool.back());
		root_pool.pop_back();
		// std::cout << roots_.back() << std::endl;
	}
	/*
	for (Vector2d root : roots_)
		std::cout << root.transpose() << std::endl;
	*/
}

void Ball::find_round_roots()
{
	////////// Find root candidates
	double w = width_ / 2.;
	double h = height_ / 2.;
	//std::vector<Vector2d> root_pool;
	roots_.clear();
	double r = e_field_.inner_r_;
	Filament root_pool = Filament();
	int max_y = 0;
	int max_yi = 0;
	for (int y = 0; y < height_; ++y)
	{
		for (int x = 0; x < width_; ++x)
		{
			if (!test_outside_ring(Vector2d(x, y), r) &&
				test_outside_ring(Vector2d(x, y), r - 1))
			{
				root_pool.points_.push_back(Vector2d(y, x));
			}
		}
	}
	double remaining_voltage = voltage_ - 50;
	double filament_voltage = 4;

	root_pool.find_nodes(e_field_.mask_);

	while (remaining_voltage > threshold_)
	{

		e_field_.iter_laplace(platform_, device_, context_, 5);
		double sum = calc_sum();
		root_pool.phis_.clear();
		MatrixXd next_map = e_field_.phi_map_; // +e_field_.rho_map_;
		for (const Vector2d& node : root_pool.nodes_) // for each node for the Filament add the fitting phi value to the phis_
		{
			root_pool.phis_.push_back(next_map(node[0], node[1])); // for the nodes add the relevant phis (including offset due to buffer on phi_map)
		}
		root_pool.find_next(sum, 4);
 		roots_.push_back(root_pool.new_point_);
				
	  if(debug_)
			show(e_field_.mask_, "roots_map");
		e_field_.add_charge(root_pool.new_point_, e_field_.charge_);

		remaining_voltage -= filament_voltage;
	}
	e_field_.phi_map_ = e_field_.mask_;
}

bool Ball::test_outside_ring(Vector2d point, double r)
{
	double w = width_ / 2.;
	double h = height_ / 2.;
	return (point(0) - h) * (point(0) - h) + (point(1) - w) * (point(1) - w) >= r * r;
}

void Ball::make_plasma()
{	
	// Step leader
	std::cout << "Step_leaders" << std::endl;
	init_plasma();
	run_field(height_, false);
	// e_field_.conv_ = (std::min)(e_field_.conv_, calc_sum()/100.0); // set convergence to when changes are less than 1/100 * abs_sum of field
	sim_Filaments();
	show(e_field_.mask_, "Leaders_");
	laplace_ = false;
	// Dart leader
	std::cout << "Dart_leaders" << std::endl;
	for (int its = 1; its <= time_; ++its)
	{
		std::cout << "Epoch number: " << its << " of " << time_ << std::endl;

		alive_ = true;
		create_pos_residual();
		if (false)//is_ball_)
			find_roots();
		init_plasma();
		run_field(height_ * 2 , true);
		sim_Filaments();
		//create_splines();
		show(e_field_.mask_, "Leaders_");
	}
}

/*
Execute the running of the DBM loop, run the field, find and add the next point
Check whether a Filament is finished and update the field accordingly
*/
void Ball::sim_Filaments(double w) 
{
	int counter = 0;
	while (alive_ == true)
	{
		counter++;
		std::cout << "Iteration: " << counter << std::endl;
		run_field(height_/5 , false, counter); // calculate the E_field
		
		if (true)
			find_next(e_field_.phi_map_);										// find the new point according to DBM model
		else
		{
			MatrixXd next_map = e_field_.phi_map_; // +e_field_.rho_map_;
			show(next_map, "next_map");
			find_next(next_map);
		}

		//	update_gaussian();

		update_alive();									// check if new point is in gas
		make_next();										// if it is then add it to the filament
		update_field();									// if it is then add to field as edge

		if (debug_)
		{
			show(e_field_.phi_map_, "phi_map");
			show(e_field_.mask_);
		}
	}
}

/*
initiate Filaments with the root points
*/
void Ball::init_plasma()
{
	for (Eigen::Vector2d root : roots_)
	{
		Filaments_.push_back(Filament(root));
	}
}

/*
update the map and iterate laplace
*/
void Ball::run_field(int its, bool reset, double w) {
	e_field_.update_mask(Filaments_);
	if(GPU_)
		e_field_.iter_laplace(platform_, device_, context_, its, reset, w);
	else
		e_field_.iter_laplace(w);
	//std::cout << e_field_.phi_map_ << std::endl;
}

/*
Calculate the sum of the phi_map
For each active Filament, find the nodes, push the phi values for the nodes onto each Filaments phis_ vector
For each Filament execute find_next to find the next penetration point
*/
void Ball::find_next(MatrixXd map)
{
	double sum = calc_sum();
	for (Filament& Filament : Filaments_)
	{	
		test_filament_join(Filament);
		if (Filament.alive_)
		{
			Filament.nodes_.clear();
			Filament.find_nodes(e_field_.mask_);
			Filament.phis_.clear(); // empty out old phis

			for (const Vector2d& node : Filament.nodes_) // for each node for the Filament add the fitting phi value to the phis_
			{
				Filament.phis_.push_back(map(node[0], node[1])); // for the nodes add the relevant phis (including offset due to buffer on phi_map)
			}
			assert(Filament.phis_.size() == Filament.nodes_.size());

			Filament.find_next(sum, dim_);

		}
	}
}

/*
for each Filament execute extend to 
*/
void Ball::make_next()
{
	for (Filament& Filament : Filaments_)
	{
		if (Filament.alive_)
		{
			Filament.extend();
		}
	}
}

/*
calculate the sum of the phi_map_ after taking each point to the power of dim_
*/
double Ball::calc_sum()
{ 
	assert(ArrayXXd(e_field_.phi_map_).pow(dim_)(0) == std::pow(e_field_.phi_map_(0), dim_));
	return ArrayXXd(e_field_.phi_map_).abs().pow(dim_).sum();
}

void Ball::update_alive() 
{
	alive_ = true;
	for (int i = 0; i < Filaments_.size(); ++i)
	{
		if (is_ball_)
		{
			if (!e_field_.test_in_gas(Filaments_[i].new_point_) || !Filaments_[i].alive_)
			{
				update_dead_Filament(i);
			}

		} 
		else
		{
			if (e_field_.mask_(Filaments_[i].new_point_(0), Filaments_[i].new_point_(1)) == 1)
			{
				update_dead_Filament(i);
			}
			
		}
	}
	if (Filaments_.begin() == Filaments_.end()) // no Filaments left
		alive_ = false;
}

void Ball::test_filament_join(Filament& fila)
{

	for (Filament& filb : Filaments_)
	{
		// if a filament has a new point in the viscinity of another kill that filament
		if (fila.root_ != filb.root_ 
			&& std::find(filb.nodes_.begin(), filb.nodes_.end(), fila.new_point_) != filb.nodes_.end()
			&& filb.alive_ == true
			)
		{
			//std::cout << "a: " << fila.new_point_(0) << ", " << fila.new_point_(1) << " b: " << filb.new_point_(0) << ", " << filb.new_point_(1) << std::endl;
			fila.alive_ = false;
			filb.charge_ += fila.charge_;
		}
	}

}

/*
i is dead, remove from Filaments_ and move to fin_Filaments_
*/
void Ball::update_dead_Filament(int i)
{
	Filaments_[i].alive_ = false;
	fin_Filaments_.push_back(Filaments_[i]);
	Filaments_.erase(Filaments_.begin() + i);
}

/*
update map of field with new points which are edge values
*/
void Ball::update_field() {
	for (Filament filament : Filaments_)
	{
		if (filament.alive_)
		{
			e_field_.add_charge(filament.new_point_, e_field_.charge_);
		}
	}
}

/*
refill all Filaments with positive residual
*/
void Ball::create_pos_residual(double w)
{
	if (is_ball_)
	{
		e_field_.create_round_field(w);
		e_field_.pos_res_map_ = e_field_.blank_map_;
		e_field_.rho_map_ = e_field_.blank_map_;
		
		double pos = 0.5;
		for (const auto& Filament : fin_Filaments_)
			for (Vector2d xy : Filament.points_)
			{
				if (e_field_.test_in_gas(xy))
				{
					e_field_.pos_res_map_(xy(0), xy(1)) = pos;
					e_field_.rho_map_(xy(0), xy(1)) = pos;
					//e_field_.phi_map_(xy(0), xy(1)) = e_field_.charge_;
				}
				// e_field_.pos_res_map_(xy(0), xy(1)) = gauss_*w;
				// e_field_.mask_(xy(0), xy(1)) = 0.5;
			}
		//show(e_field_.rho_map_, "rho_pre_gauss");
	if (gauss_ < 0)
	{
		for (int i = 0; i < gauss_; i++)
		{
			e_field_.gaussian_blur(e_field_.rho_map_);
			//e_field_.gaussian_blur(e_field_.pos_res_map_);
		}
	}
		//e_field_.rho_map_ = e_field_.rho_map_ - e_field_.base_gauss_map_;

		/*
		for (int i = 0; i < e_field_.rho_map_.cols(); ++i)
		{
			for (int j = 0; j < e_field_.rho_map_.rows(); ++j)
			{
				if (e_field_.test_in_gas(Vector2d(i, j)) && e_field_.pos_res_map_(i, j) != 0)
				{
					e_field_.phi_map_(i, j) = e_field_.pos_res_map_(i, j);
					e_field_.mask_(i, j) = 1;
				}
			}
		}
		e_field_.rho_map_ = e_field_.pos_res_map_;
		*/
		//e_field_.clamp(e_field_.rho_map_);

	}
	else
	{
		e_field_.create_rect_field();
		for (const auto& Filament : fin_Filaments_)
			for (Vector2d xy : Filament.points_)
			{
				e_field_.rho_map_(xy(0), xy(1)) = 0.5;
			}
	}
	if (debug_)
		{
			show(e_field_.rho_map_, "gaussian");
			show(e_field_.pos_res_map_, "pre-gaussian");
			//show(e_field_.phi_map_, "phi_post_gauss");
			//show(e_field_.pos_res_map_, "pos_res");
		}

		fin_Filaments_.clear();

}

void Ball::update_gaussian(double w)
{
	for (const auto& Filament : Filaments_)
	{
		e_field_.pos_res_map_(Filament.root_(0), Filament.root_(1)) = -w*e_field_.charge_;
		e_field_.pos_res_map_(Filament.new_point_(0), Filament.new_point_(1)) = -w*e_field_.charge_;
	}
	e_field_.rho_map_ = e_field_.pos_res_map_;
	for (int i = 0; i < gauss_; i++)
		e_field_.gaussian_blur(e_field_.rho_map_);

	e_field_.clamp(e_field_.rho_map_);

	if (debug_)
		show(e_field_.rho_map_, "rho_post_gauss");
}
/*
calculate splines for each Filament
*/
void Ball::create_splines()
{
	for (Filament Filament : Filaments_)
		splines_.push_back(Cubic_Spline(Filament.nodes_));
}

Ball::~Ball()
{
}

std::vector<std::vector<char>> Ball::make_show_map()
{
	std::vector<std::vector<char>> display(height_, std::vector<char>(width_));
	MatrixXd blank_map = MatrixXd::Ones(height_ + 2, width_ + 2);
	e_field_.create_mask(blank_map);
	//MatrixXd blank_map = e_field_.mask_;
	for (int i = 0; i < height_ - 1; i++)
	{
		for (int j = 0; j < width_; j++)
		{
			display[i][j] = (blank_map(i, j) == -1.0) ? '-' : '.' ;
			// display[i][j] = (blank_map(i, j) >= 10.0) ? '$' : display[i][j];
			display[i][j] = (blank_map(i, j) == 1.0) ? '+' : display[i][j];
		}
	}

	if (debug_) 
	{
		for (Filament &Filament : Filaments_)
		{
			for (Vector2d point : Filament.points_)
			{
				display[point(1)][point(0)] = '#';
				// std::cout << display[node(1)][node(0)] << std::endl;
			}
			display[Filament.new_point_(1)][Filament.new_point_(0)] = 'T';
		}
		if (laplace_)
			for (Filament Filament : fin_Filaments_)
				for (Vector2d point : Filament.points_)
					display[point(1)][point(0)] = '#';
	}  
	if (!alive_)
		for (int i = fin_Filaments_.size() / 2; i < fin_Filaments_.size(); i++)
		{
			Filament Filament = fin_Filaments_[i];
			for (Vector2d point : Filament.points_)
			{
				display[point(1)][point(0)] = '#';
				// std::cout << display[node(1)][node(0)] << std::endl;
			}
		}
	return display;
}

void Ball::show_raw(MatrixXd &Mat)
{
	for (int i = 0; i < Mat.cols(); i++)
	{
		for (int j = 0; j < Mat.rows(); j++)
		{
			std::cout << Mat(i, j) << " ";
		}
		std::cout << std::endl;
	}
}

/*
name: name of process to be shown
write to text and to vtk files if required, if in debug mode show each step
*/
void Ball::show( MatrixXd &Map, std::string name)
{
	std::string filename;
	if(is_ball_)
		filename = name + "r-" + std::to_string(radius_) + "_Filaments-" + std::to_string(roots_.size()) + "_dim-" + std::to_string(dim_);
	else
		filename = name + "whm-" + std::to_string(radius_) + "_Filaments-" + std::to_string(roots_.size()) + "_dim-" + std::to_string(dim_);
	if (false)
	{	
		std::vector<std::vector<char>> display = make_show_map();
		// write_txt_(display, filename);
		for (int i = 0; i < height_; i++)
		{
			for (int j = 0; j < width_; j++)
			{
				std::cout << display[i][j] << " ";
			}
			std::cout << std::endl;
		}
	} 

	if (vtk_)
		write_vtk_(filename, name, Map);

}

void Ball::write_txt_(std::vector<std::vector<char>> display, std::string filename)
{
	std::ofstream fs("../txt_images/" + filename + ".txt");
	if (!fs)
	{
		std::cerr << "Cannot open the output file." << std::endl;
	}
	else
		std::cout << "Writing " +  filename + " as txt" << std::endl;
		fs << std::to_string(height_) << std::endl << std::endl;

		for (int i = 0; i < height_ - 1; i++)
		{
			for (int j = 0; j < width_; j++)
			{

				fs << display[i][j] << " ";
			}
			fs << std::endl;
		}
		fs.close();
}


void Ball::write_vtk_(std::string filename, std::string name, MatrixXd &Map)
{
	show_counters_[name]++;

	std::vector<char> writearray;
	std::string it;
	it = std::to_string(show_counters_[name]);
	it.insert(it.begin(), 4 - it.size(), '0');
	std::cout << "Writing " + name + " it " + it + " as vtk" << std::endl;
	for (int n = 0; n < 1; n++)
	{
		std::string dirname = (name == "Leaders_")? "../vtk_images/" : "../vtk_images/debug/";
		std::string fname = dirname + filename + "_it-" + it + ".vtk";

		int POINTS = (Map.cols() ) * (Map.rows() );
		int CELLS = (Map.cols()-1) * (Map.rows()-1);
		std::ofstream file(fname, std::ios::out | std::ios::trunc);
		if (file)
		{
			// /* Structured grid
			file << "# vtk DataFile Version 2.0\n";
			file << "Plasma Ball, radius = " << radius_ << ", eta = " << dim_ << "\n";
			file << "ASCII" << "\n";
			file << "DATASET STRUCTURED_POINTS\n";
			file << "DIMENSIONS " << Map.cols() << " " << Map.rows() << " 1" << "\n";
			//file << "POINTS " << Map.cols() * Map.rows() << "\n";
			file << "ASPECT_RATIO 1 1 1" << "\n";
			file << "ORIGIN 0 0 0" << "\n";
			file << "POINT_DATA " << POINTS << "\n";
			file << "SCALARS volume_scalars float 1\n";
			file << "LOOKUP_TABLE default";
			// */
		  /* Rectilinear grid
			file << "# vtk DataFile Version 3.0\n";
			file << "Plasma Ball, radius = " << height_ << ", eta = " << dim_ << "\n";
			file << "ASCII" << "\n";
			file << "DATASET RECTILINEAR_GRID\n";
			file << "DIMENSIONS " << Map.cols() << " " << Map.rows() << " 1" << "\n";
			file << "X_COORDINATES " << Map.cols() << " float\n";
			for (int x = 0; x < Map.cols() ; ++x)
				file << x << " ";
			file << "\n";
			file << "Y_COORDINATES " << Map.rows() << " float\n";
			for (int x = 0; x < Map.rows() ; ++x)
				file << x << " ";
			file << "\n";
			file << "Z_COORDINATES " << 1 << " float\n";
			for (int x = 0; x < 0 + 1; ++x)
				file << x << " ";
			file << "\n";
			file << "CELL_DATA " << CELLS << "\n";
			file << "POINT_DATA " << POINTS << "\n";
			file << "FIELD FieldData 1\n";
			file << "nodal 1 " << CELLS << " float";
			 */
			int rowindent = 0;
			for (int i = 0; i < Map.rows(); ++i)
			{
				for (int j = 0; j < Map.cols(); ++j)
				{
					if (rowindent % 9 == 0) // 9 points per row (convention)
						file << "\n";
					file << Map(i, j) << " ";
					rowindent++;
				}
			}
			file.close();
		}
		else {
			file.clear();
			file.close();
			std::cout << "Error writing " << iteration_ << " to vtk" << std::endl;
		}
	}

}