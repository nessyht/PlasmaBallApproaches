#include "Field.h"
#include <thread>

using namespace Eigen;

Field::Field(bool ball, bool GPU, int gauss, int radius, int pixels, int depth) : 
	width_(2 * radius), height_(2 * radius), depth_(2 * depth), numCols_((2 * radius)+2), numRows_((2 * radius)+2), count_(((2*radius)+2)*((2 * radius) + 2))
{
	if (GPU)
		laplace_program_ = init_GPU("itlaplacekernelold.cl");
	if (ball)
	{
		create_round_field();
	}
	else
	{
		create_rect_field();
	}
	//rho_map_ = MatrixXd::Zero(width_ + 2, height_ + 2);
	base_gauss_map_ = rho_map_;
	for (int i = 0; i < gauss; i++)
	{
		if (ball)
			gaussian_blur(base_gauss_map_);
	}
		


	
	mesh_ = MatrixXd::Zero(radius, radius);
	//field_ = RowVectorXcd::LinSpaced(height*width, 0, height*width - 1).replicate(height, 1);
	//std::cout << field_.rows() << ", " << field_.cols() << std::endl;
	int n = 0;

	for (int i = 0; i < radius; ++i)
		for (int j = 0; j < radius; ++j)
			mesh_(i, j) = radius * i + j;


	double conv_ = 0.01; // 0.00001 / ((width + height)*0.5);
}

void Field::create_rect_field()
{

	map_ = MatrixXd::Zero(width_ + 2, height_ + 2)*-1;
	map_.block(0, 1, width_ + 2, height_) = MatrixXd::Zero(width_ + 2, height_); // core
	map_.block(height_ + 1, 0, 1, width_ + 2) = MatrixXd::Ones(1, width_ + 2) * -1; // Bot
	map_.block(0, 0, 1, width_ + 2) = MatrixXd::Ones(1, width_ + 2) * 1;  // Top
	// map_.block(0, height_ + 2, width_ + 2, 1) = MatrixXd::Ones(width_ + 2, 1);		// Right edge
	// phi_map_ = MatrixXd::Ones(width_ + 2, height_ + 2) * 0.5; //phi_map_ has width+2, heigh+2 due to buffers
	phi_map_ = map_;
	mask_ = map_;
	rho_map_ = map_;
	//rho_map_.block(height_, 0, 1, width_ + 2) = MatrixXd::Ones(1, width_ + 2) * -1; // Bot
	//rho_map_.block(1, 0, 1, width_ + 2) = MatrixXd::Ones(1, width_ + 2);  // Top
}

void Field::create_round_field(double w)
{
	inner_r_ = fmax(fmin(width_, height_) / 15.0, 1);
	outer_r_ = (fmin(width_, height_) / 2.0) - 1;
	map_ = MatrixXd::Ones(width_ + 2, height_ + 2); // map_ has width + 2, heigh + 2 due to buffers contains a base map for valid zones and base charges
	blank_map_ = MatrixXd::Zero(width_ + 2, height_ + 2);
	//map_.block(0, 0, 3, width_ + 2) = MatrixXd::Ones(3, width_ + 2) * 3;  // Finger on top 

	create_mask(map_); // Creates valid zones that initially are the mask_ with no Filaments
	phi_map_ = map_; // Contains the potentials
	mask_ = map_; // Contains only the zones that potentials can be calculated in, Filament points are added to edges
	rho_map_ = map_; // Contains potentials
	pos_res_map_ = map_; // currently unused outside of obselete update gaussian method

}

void Field::reset_mask()
{
	mask_ = map_;
}

bool Field::test_in_gas(Vector2d point)
{
	return map_(point(0), point(1)) == 0;
}

void Field::create_ring(double r, double fill, MatrixXd &map)
{
	double w = height_ / 2.;
	double h = width_ / 2.;

	for (int y = 0; y < width_; ++y)
	{
		for (int x = 0; x < height_; ++x)
		{
			if (sqrt(double(y - h) * double(y - h) + double(x - w) * double(x - w)) < r)
			{
				map(y, x) = fill;
			}
		}
	}
}

void Field::create_ring_outer(double r, double fill, MatrixXd &map)
{
	double w = height_ / 2.;
	double h = width_ / 2.;

	for (int y = 0; y < width_; ++y)
	{
		for (int x = 0; x < height_; ++x)
		{
			if (sqrt(double(y - h) * double(y - h) + double(x - w) * double(x - w)) >= r)
			{
				map(y, x) = fill;
			}
		}
	}
}

void Field::create_mask(MatrixXd &map)
{
	create_ring(outer_r_, 0, map);
	create_ring(inner_r_, charge_, map);
}

void Field::add_charge(Vector2d point, double charge)
{
	// map_(point(0), point(1)) = charge; map is the base zones of material
	mask_(point(0), point(1)) = charge_; // add edge to mask
	phi_map_(point(0), point(1)) = charge; // add charge to edge of mask
}

/*
iterate through points of the Filaments and add them to the map
*/
void Field::update_mask(const std::vector<Filament>& Filaments)
{
	for (const Filament& filament : Filaments)
	{
		//To Do: load all points into field
		for (const Vector2d& point : filament.points_)
		{
			//map_(point(0), point(1)) = charge_;
			add_charge(point, charge_);
		}
	}
}

bool Field::test_outside_ring(Vector2d point, double r)
{
	double w = height_ / 2.;
	double h = width_ / 2.;
	return (point(0) - h) * (point(0) - h) + (point(1) - w) * (point(1) - w) >= r * r;
}

/*
while the process hasnt converged, execute the laplacian on the phi map
*/
void Field::iter_laplace(double w) {
	// phi_map_ = map_;
	// mask_ = map_;
	int counter = 0;
	while (!convergence_)
	{	
		std::cout << ++counter << ", ";
		old_phi_map_ = phi_map_;
		la_map_ = phi_map_;
		// laplace(w); // single threaded
		laplace_threaded(w); // 8 threaded
		phi_map_ = la_map_;
		calculate_convergence();

		/* debug
		std::cout << inner_r_ << "  " << outer_r_ << std::endl;
		std::cout << "mask" << std::endl;
		for (int i = 0; i < mask_.rows(); ++i)
		{
			for (int j = 0; j < mask_.cols(); ++j)
				std::cout << mask_(i, j) << " ";
			std::cout << std::endl;
		}
		// */
		/* debug 
		std::cout << "new unmasked" << std::endl;
		for (int i = 0; i < phi_map_.rows(); ++i)
		{
			for (int j = 0; j < phi_map_.cols(); ++j)
			{
				std::cout << phi_map_(i, j) << " ";
			}
			std::cout << std::endl;
		}
		// */
	}
	convergence_ = false;
}

void Field::laplace(double w)
{
	calc_laplace(height_, 1, w);
}

void Field::laplace_threaded(double w)
{	
	constexpr int threads_to_use = 8;
	std::thread t[threads_to_use];	
	double region_cols = (map_.cols()-2) / threads_to_use;

	for (size_t i = 0; i < threads_to_use; i++)
	{
		// std::cout << std::endl << region_cols*i + 1 << ", " << region_cols * (i + 1) << std::endl; // debug
		t[i] = std::thread(&Field::calc_laplace, 
					this,
					region_cols, 
					i, //start at the first non buffer element
					w);
	}

	for (int i = 0; i < threads_to_use; ++i)
	{
		t[i].join();
	}
}


void Field::calc_laplace(double width, int slot,  double w) // laplace with w = 1, poisson with w = rho
{
	for (int j = 1+slot*width; j <= (slot+1)*width; j++)
	{
		for (int i = 1; i < numCols_; i++)		
		{
			if(map_(i, j) == 0)
				// la_map_(i, j) = w * 0.25 * (phi_map_(i - 1, j) + phi_map_(i + 1, j) + phi_map_(i, j - 1) + phi_map_(i, j + 1));
				la_map_(i, j) = phi_map_(i, j) + 
																				w * 1.5 * (
																										(
																											2 * (phi_map_(i - 1, j) + phi_map_(i + 1, j) + phi_map_(i, j - 1) + phi_map_(i, j + 1)) 
																											+ phi_map_(i-1, j-1) + phi_map_(i-1, j+1) + phi_map_(i+1, j-1) + phi_map_(i+1, j+1)
																										)
																										/12 - phi_map_(i, j)
																									);
				// la_map_(i, j) = w * 0.25 * (phi_map_(i - 1, j) + phi_map_(i + 1, j) + phi_map_(i, j - 1) + phi_map_(i, j + 1)) - phi_map_(i, j);		
				// la_map_(i, j) = (phi_map_(i - 1, j - 1) + phi_map_(i + 1, j - 1) + phi_map_(i - 1, j + 1) + phi_map_(i + 1, j + 1)) / 4 - phi_map_(i, j);
		}
	}
}

/*
int its: maximum iterations allowed
bool reset: whether the previous phi_map should be used, if not: reset to base mask
double w: depreciated
*/
void Field::iter_laplace(cl::Platform &platform, cl::Device &device, cl::Context &context, int its, bool reset, double w) {
	int counter = 0;
	if (reset)
		phi_map_ = mask_;
	while (!convergence_ && counter < its)
	{
		std::cout << ++counter << ", ";
		old_phi_map_ = phi_map_; // for convergence
		// la_map_ = phi_map_; // for copy

		laplace_GPU(w);
		/* comparison debug

		laplace_threaded(w);
		auto CPU_map = la_map_;

		laplace_GPU(w);
		auto GPU_map = la_map_;
		//*/

		// phi_map_ = la_map_; // copy back
		calculate_convergence();
		/* debug
		std::cout << "phi_map" << std::endl;
		for (int i = 0; i < mask_.rows(); ++i)
		{
			for (int j = 0; j < mask_.cols(); ++j)
				std::cout << phi_map_(i, j) << " ";
			std::cout << std::endl;
		}
		// */
		/* comparison debug

		std::cout << "GPU_map" << std::endl;
		for (int i = 0; i < mask_.rows(); ++i)
		{
			for (int j = 0; j < mask_.cols(); ++j)
				std::cout << GPU_map(i, j) << " "; // phi_map_(i, j) << " ";
			std::cout << std::endl;
		}
		// */
		 /* debug
		std::cout << "CPU phi_map" << std::endl;
		for (int i = 0; i < phi_map_.rows(); ++i)
		{
			for (int j = 0; j < phi_map_.cols(); ++j)
			{
				std::cout << CPU_map(i, j) << " ";
			}
			std::cout << std::endl;
		}
		// */
		/* debug
		std::cout << "GPU - CPU phi_map" << std::endl;
		for (int i = 0; i < phi_map_.rows(); ++i)
		{
			for (int j = 0; j < phi_map_.cols(); ++j)
			{
				std::cout << GPU_map(i, j) - CPU_map(i, j) << " ";
			}
			std::cout << std::endl;
		}
		// */
	}
 	convergence_ = false;
	 /* debug
	std::cout << "Final phi_map" << std::endl;
	for (int i = 0; i < mask_.rows(); ++i)
	{
		for (int j = 0; j < mask_.cols(); ++j)
			std::cout << phi_map_(i, j) << " ";
		std::cout << std::endl;
	}
	
	std::cout << "dif\n";
	//*/
}

cl::Program Field::init_GPU(std::string name)
{
	try {
		std::vector<cl::Platform> all_platforms;
		cl::Platform::get(&all_platforms);
		if (all_platforms.size() == 0) {
			std::cout << " No platforms found. Check OpenCL installation!\n";
			exit(1);
		}
		cl::Platform platform = all_platforms[0];
		std::cout << "Using platform: " << platform.getInfo<CL_PLATFORM_NAME>() << "\n";

		//get default device of the default platform
		std::vector<cl::Device> all_devices;
		platform.getDevices(CL_DEVICE_TYPE_GPU, &all_devices);
		if (all_devices.size() == 0) {
			std::cout << " No devices found. Check OpenCL installation!\n";
			exit(1);
		}
		cl::Device device = all_devices[0];
		std::cout << "Using device: " << device.getInfo<CL_DEVICE_NAME>() << "\n";
		std::cout << "Max workgroup size : " << device.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>() << "\n";
		cl::Context context(device); // null pointer read here

		std::ifstream kernel_file(name);
		std::string src(std::istreambuf_iterator<char>(kernel_file), (std::istreambuf_iterator<char>()));

		cl::Program::Sources sources(1, std::make_pair(src.c_str(), src.length() + 1));
		std::cout << "\n errs: " << (bool)kernel_file << std::endl;
		cl::Program program(context, sources);
		cl_int err = program.build();
		std::cout << err << "\n";
		gpu_is_valid_ = true;
		return program;
	}
	catch (cl_int err)
	{
		std::cout << "Error in GPU initialisation, error code: " << err << "\n";
		gpu_is_valid_ = false;
	}
}

void Field::clamp(MatrixXd& mat)
{
	float minimum = abs(mat.minCoeff());
	float max = mat.maxCoeff() + minimum;

	for (int i = 0; i < mat.rows(); i++)
	{
		for (int j = 0; j < mat.cols(); j++)
		{
			mat(i, j) = (map_(i, j) == 0) ? (2 * (mat(i, j) + minimum) / max) - 1 : 0;
		}
	}
}

void Field::make_matrix()
{
	const int rows = mask_.rows();
	const int cols = mask_.cols();
	double dx2 = 1. / (rows*rows);
	double diag = -4 / dx2;
	double offdiag = 1 / dx2;
	MatrixXd A = MatrixXd::Zero(rows*cols, rows*cols);
	for (int i = 0; i < cols*rows; ++i)
	{
		for (int j = 0; j < cols*rows; ++j)
		{
			if (i == j)
				A(i, j) = diag;
			else if (i == j + 1)
				A(i, j) = offdiag;
			else if (i == j - 1)
				A(i, j) = offdiag;
			else if (i == j + rows)
				A(i, j) = offdiag;
			else if (i == j - rows)
				A(i, j) = offdiag;

		}
	}
}

void Field::laplace_GPU(double w)
{
	//laplace_program_ = init_GPU();
	auto context = laplace_program_.getInfo<CL_PROGRAM_CONTEXT>();
	auto devices = context.getInfo<CL_CONTEXT_DEVICES>();
	auto device = devices.front();

	if (gpu_is_valid_)
	{

		const int numCols = 256;
		const int numRows = 256;

		const int count = numCols * numRows;
		const int bufsize = sizeof(int) * count;

		float min = abs(phi_map_.minCoeff());
		float max = phi_map_.maxCoeff() + min;
		//std::cout << "min: " << min << "\n" << max << "\n";

		std::array<std::array<float, numCols>, numRows> laplace_arr;
		std::array<std::array<float, numCols>, numRows> mask;
		std::array<std::array<float, numCols>, numRows> rho_map;
		std::array<std::array<float, numCols>, numRows> outarr;

		for (int i = 0; i < numRows_; i++)
		{
			for (int j = 0; j < numCols_; j++)
			{
				//laplace_arr[i][j] = (mask_(i, j) == 1)? 1 : (phi_map_(i, j) + min) / max;
				//laplace_arr[i][j] = (phi_map_(i, j) + min) / max;
				laplace_arr[i][j] = phi_map_(i, j);
				mask[i][j] = mask_(i, j);
				//rho_map[i][j] = (rho_map_(i, j) == 0) ? 0 : (rho_map_(i, j) + min) / max;
				rho_map[i][j] = rho_map_(i, j);
				//rho_map[i][j] = (rho_map_(i, j) + min) / max;
				//std::cout << laplace_arr[i][j] << " ";
			}
			//std::cout << "\n";
		}
		cl_int err;
		//cl::Buffer b;
		//b = cl::Buffer(laplace_arr.begin(), laplace_arr.end(), true, true, &err);

		//std::cout << "\nerrs:\n";
		cl::Buffer buf(context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_COPY_HOST_PTR, bufsize, laplace_arr.data(), &err);
		//std::cout << err << std::endl;
		cl::Buffer mbuf(context, CL_MEM_READ_ONLY | CL_MEM_HOST_READ_ONLY | CL_MEM_COPY_HOST_PTR, bufsize, mask.data(), &err);
		//std::cout << err << "\n";
		cl::Buffer rhobuf(context, CL_MEM_READ_ONLY | CL_MEM_HOST_READ_ONLY | CL_MEM_COPY_HOST_PTR, bufsize, rho_map.data(), &err);
		//std::cout << err << "\n";
		cl::Buffer outbuf(context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY, bufsize, nullptr, &err);
		//std::cout << err << std::endl;
		cl::Kernel kernel(laplace_program_, "itpoisson", &err);
		//std::cout << err << std::endl;
		err = kernel.setArg(0, buf);
		//std::cout << err << std::endl;
		err = kernel.setArg(1, mbuf);
		//std::cout << err << std::endl;
	  err = kernel.setArg(2, rhobuf);
		//std::cout << err << std::endl;
		err = kernel.setArg(3, outbuf);

		cl::CommandQueue queue(context, device);
		err = queue.enqueueNDRangeKernel(kernel, cl::NullRange, cl::NDRange(numCols, numRows));
		//std::cout << err << std::endl;
		err = queue.enqueueReadBuffer(outbuf, GL_TRUE, 0, bufsize, outarr.data());
		//std::cout << err << std::endl;
		err = queue.finish();
		//std::cout << err << std::endl;
		//std::cout << "width: " << outarr.size() << "height: " << outarr[0].size() << std::endl;
		for (int i = 0; i < numRows_; i++)
		{
			for (int j = 0; j < numCols_; j++)
			{
				//phi_map_(i, j) = (mask_(i, j) == 0)? (outarr[i][j] * max) - min : mask_(i, j);
				//std::cout << phi_map_(i, j) << " ";
				phi_map_(i, j) = outarr[i][j];
			}
			//std::cout << "\n";
		}
		/*
		min = abs(phi_map_.minCoeff());
		max = phi_map_.maxCoeff() + min;

		for (int i = 0; i < numRows_; i++)
			for (int j = 0; j < numCols_; j++)
				phi_map_(i, j) = (map_(i, j) == 0)? (2 * (phi_map_(i, j) + min) / max) - 1 : map_(i, j);
		// */
	}
}

void Field::calculate_convergence() 
{
	convergence_ = (old_phi_map_ - phi_map_).squaredNorm() < conv_;
}


void Field::gaussian_blur(MatrixXd& map)
{
	//std::cout << "map: " << map_.rows() << ", " << map_.cols() << std::endl;
	//std::cout << "height, width " << numRows_ << ", " << numCols_ << std::endl;
	MatrixXd temp = map;
	for (int i = 1; i < numRows_ -2; ++i)
	{
		for (int j = 1; j < numCols_ -2; ++j)
		{
			
			if (map_(i, j) == 0) // in all positions in the gas
			{
				temp(i, j) = (1. / 256.)*((map(i - 2, j - 2) + map(i + 2, j + 2) + map(i - 2, j + 2) + map(i + 2, j - 2)) +
														  4.*(map(i - 1, j - 2) + map(i + 1, j + 2) + map(i - 1, j + 2) + map(i + 1, j - 2) +
																 map(i - 2, j - 1) + map(i + 2, j + 1) + map(i - 2, j + 1) + map(i + 2, j - 1)) +
														  6.*(map(i, j - 2) + map(i, j + 2) + map(i - 2, j) + map(i + 2, j)) +
													   16.*(map(i - 1, j - 1) + map(i + 1, j + 1) + map(i - 1, j + 1) + map(i + 1, j - 1)) +
														 24.*(map(i - 1, j) + map(i + 1, j) + map(i, j + 1) + map(i, j - 1)) +
														 36.*(map(i, j)));			
				//std::cout << temp(i, j) << " ";
			}
		}
		//std::cout << "\n";
	}
	map = temp;
	//rho_map_ = rho_map_ - base_gauss_map_;
}

Vector2d Field::get_coords(int point) 
{
	return Vector2d(int(point / height_), point % height_);
}

Field::~Field()
{
}
