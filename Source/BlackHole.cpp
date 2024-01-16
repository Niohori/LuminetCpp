#include "BlackHole.h"

BlackHole::BlackHole() {
	;
}
BlackHole::BlackHole(double mass_, double inclination_, double acc_) :
	inclination(inclination_),
	M(mass_),
	acc(acc_)
{
	t = inclination_ * M_PI / 180;
	disk_outer_edge = 50. * M;
	disk_inner_edge = 6. * M;

	//isoradials = {};
	//isoredshifts = {};
	/*initial_guesses = 10;
	midpoint_iterations = 10;
	plot_inbetween = false;*/
	min_periastron = irs_solver_params_.min_periastron * M;
	//use_ellipse = true;

	critical_b = 3.0 * std::sqrt(3) * M;
	angular_precision = 100;
}

BlackHole::~BlackHole() {
	;
}

/*
----------------------------------------------------------------------------------------------------------------
 ?
----------------------------------------------------------------------------------------------------------------
*/
Isoradial BlackHole::calc_apparent_outer_disk_edge() {
	// Implementation of calc_apparent_outer_disk_edge method
	Isoradial ir(disk_outer_edge, t, M, 0);
	std::vector<double> X_;
	std::vector<double> Y_;
	std::pair<std::vector<double>, std::vector<double>> XY_;
	XY_ = OperatorsOrder2::polar_to_cartesian_lists(ir.radii_b, ir.angles, -M_PI / 2);
	ir.X = std::get<0>(XY_);
	ir.Y = std::get<1>(XY_);
	return ir;
}

/*
----------------------------------------------------------------------------------------------------------------
 ?
----------------------------------------------------------------------------------------------------------------
*/
Isoradial BlackHole::calc_apparent_inner_disk_edge() {
	// Implementation of calc_apparent_inner_disk_edge method
	Isoradial ir(disk_inner_edge, t, M, 0);
	for (auto& b : ir.radii_b) {
		b *= 0.99; // scale slightly down
	}
	std::vector<double> X_;
	std::vector<double> Y_;
	std::pair<std::vector<double>, std::vector<double>> XY_;
	XY_ = OperatorsOrder2::polar_to_cartesian_lists(ir.radii_b, ir.angles, -M_PI / 2);
	ir.X = std::get<0>(XY_);
	ir.Y = std::get<1>(XY_);
	return ir;
}

/*
----------------------------------------------------------------------------------------------------------------
?
----------------------------------------------------------------------------------------------------------------
*/
double BlackHole::get_apparent_outer_edge_radius(Isoradial& ir, double angle, double rotation) {
	return ir.get_b_from_angle(angle + rotation);
}
/*
----------------------------------------------------------------------------------------------------------------
 ?
----------------------------------------------------------------------------------------------------------------
*/
double BlackHole::get_apparent_inner_edge_radius(Isoradial& ir, double angle, double rotation) {
	return ir.get_b_from_angle(angle + rotation);
}

/*
----------------------------------------------------------------------------------------------------------------
The apparent inner edge of the black hole(not the apparent inner disk edge).THis takes eventual ghost images
into account
----------------------------------------------------------------------------------------------------------------
*/

std::pair<std::vector<double>, std::vector<double>> BlackHole::apparent_inner_edge(Isoradial& ir, bool cartesian, double scale) {
	std::vector<double> b, a;
	std::vector<double> a__ = OperatorsOrder2::linspace(0.0, 2 * M_PI, angular_precision);
	for (double a_ : a__) {
		a.push_back(a_);

		if (M_PI / 2 < a_ && a_ < 3 * M_PI / 2) {
			b.push_back(critical_b * scale);
		}
		else {
			double b_ = std::min(critical_b, ir.get_b_from_angle(a_));
			b.push_back(b_ * scale);
		}
	}

	if (!cartesian) {
		return std::make_pair(b, a);
	}
	else {
		return OperatorsOrder2::polar_to_cartesian_lists(b, a, -M_PI / 2);
	}
}


/*
----------------------------------------------------------------------------------------------------------------

----------------------------------------------------------------------------------------------------------------
*/


std::map<double, IsoRedShift> BlackHole::calc_isoredshifts(std::vector<double> redshifts) {
	
	//std::map<double, IsoRedShift> isoredshifts;
	std::map<double, std::map<int, Isoradial>> dirty_isoradials = get_dirty_isoradials();

	for (const auto& redshift : redshifts) {
		std::cout << "Calculating redshift " << redshift << std::endl;
		std::map<double, std::map<int, Isoradial>> dirty_ir_copy = dirty_isoradials;
		IsoRedShift iz(t, redshift, M, dirty_ir_copy);
		iz.improve();//to code!!!
		isoredshifts.emplace(redshift, iz);
	}

	return isoredshifts;

}
/*
----------------------------------------------------------------------------------------------------------------

----------------------------------------------------------------------------------------------------------------
*/
std::map<double, std::map<int, Isoradial>> BlackHole::get_dirty_isoradials() {
	std::map<double, std::map<int, Isoradial>> isoradials_;
	std::vector<double> a__ = OperatorsOrder2::linspace(disk_inner_edge, disk_outer_edge, irs_solver_params_.initial_radial_precision);
	for (double radius : a__) {
		Isoradial isoradial(radius, t, M, 0);// , ir_parameters);
		isoradials_[radius] = { {0, isoradial} };
	}
	return isoradials_;
}
/*
----------------------------------------------------------------------------------------------------------------
Add isoradial to dict of isoradials. Each key is a radius corresponding to
		some set of isoradials. Each value is again a dict, with as keys the order
		of the isoradial (usually just 0 for direct and 1 for ghost image)
----------------------------------------------------------------------------------------------------------------
*/
void BlackHole::add_isoradial(Isoradial& isoradial, double radius, int order) {
	if (isoradials.find(radius) != isoradials.end()) {
		isoradials[radius][order] = isoradial;
	}
	else {
		isoradials[radius] = { {order, isoradial} };
	}
}


/*
----------------------------------------------------------------------------------------------------------------

----------------------------------------------------------------------------------------------------------------
*/
void BlackHole::calc_isoradials(const std::vector<double>& direct_r, const std::vector<double>& ghost_r) {
	std::cout << "Ghost images" << std::endl;
	plot_params_.alpha = 0.5;
	for (const auto& radius : ghost_r) {
		//Isoradial isoradial(radius, t, M, 1);
		//add_isoradial(isoradial, radius, 1);
	}

	std::cout << "Direct images" << std::endl;
	plot_params_.alpha = 1.0;
	for (const auto& radius : direct_r) {
		//Isoradial isoradial(radius, t, M, 0);
		//add_isoradial(isoradial, radius, 0);
	}
}




/*
----------------------------------------------------------------------------------------------------------------
# TODO: sample separately for direct and ghost image?
		Samples points on the accretion disk. This sampling is not done uniformly, but a bias is added towards the
		center of the accretion disk, as the observed flux is exponentially bigger here and this needs the most
		precision.
		Both the direct and ghost image for each point is calculated. It's coordinates (polar and cartesian),
		redshift and
		:param min_radius:
		:param max_radius:
		:param n_points: Amount of points to sample. 10k takes about 6 minutes and gives ok precision mostly
		:param f:
		:param f2:
		:return:
----------------------------------------------------------------------------------------------------------------
*/
void BlackHole::sample_Sources(int n_points, const std::string& f, const std::string& f2) {
	//Reading data....
	//Remark/To consider: the reading of the data could be transferred to "utilites.*" but this would imply to get acces to the struct "Source" defined in "BlackHole.h"
	std::string file1 = (f.empty()) ? "D:/dev/Cpp/LuminetCpp/Points/points_incl=85.csv" : f;
	std::string file2 = (f2.empty()) ? "D:/dev/Cpp/LuminetCpp/Points/points_secondary_incl=85.csv" : f2;
	std::cout << "File 1: " << file1 << std::endl;
	std::cout << "File 2: " << file2 << std::endl;
	std::vector<Source> points;
	std::vector<Source> points2;
	std::ifstream  file(file1);
	std::string line;
	CSVRow row;
	unsigned count = 0;

	if (file.is_open())
	{
		Source point;
		while (getline(file, line))
		{
			if (count != 0) {
				row.readNextRow(line);

				point.X = std::stod(std::string{ row[1] });
				point.Y = std::stod(std::string{ row[2] });
				point.impact_parameter = std::stod(std::string{ row[3] });
				point.angle = std::stod(std::string{ row[4] });
				point.z_factor = std::stod(std::string{ row[5] });
				point.flux_o = std::stod(std::string{ row[6] });

				points.push_back(point);
				//std::cout << point.X << "  " << point.Y << "  " << point.impact_parameter << "  " << point.angle << "  " << point.z_factor << "  " << point.flux_o << std::endl;
			}
			count++;
		}
	}

	std::ifstream  file_2(file2);
	count = 0;

	if (file_2.is_open())
	{
		Source point;
		while (getline(file_2, line))
		{
			if (count != 0) {
				row.readNextRow(line);

				point.X = std::stod(std::string{ row[1] });
				point.Y = std::stod(std::string{ row[2] });
				point.impact_parameter = std::stod(std::string{ row[3] });
				point.angle = std::stod(std::string{ row[4] });
				point.z_factor = std::stod(std::string{ row[5] });
				point.flux_o = std::stod(std::string{ row[6] });

				points2.push_back(point);
				//std::cout << point.X << "  " << point.Y << "  " << point.impact_parameter << "  " << point.angle << "  " << point.z_factor << "  " << point.flux_o << std::endl;
			}
			count++;
		}
	}
	//End of reading data

	double min_radius_ = disk_inner_edge;
	double max_radius_ = disk_outer_edge;;

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dist_theta(0.0, 2.0 * M_PI);
	std::uniform_real_distribution<double> dist_radius(min_radius_, max_radius_);

	for (int i = 0; i < n_points; ++i) {
		double theta = dist_theta(gen);
		double r = dist_radius(gen);
		double b_ = BHphysics::calc_impact_parameter(r, t, theta, M, irs_solver_params_.midpoint_iterations, irs_solver_params_.plot_inbetween, 0, irs_solver_params_.min_periastron, irs_solver_params_.initial_guesses, irs_solver_params_.use_ellipse);
		double b_2 = BHphysics::calc_impact_parameter(r, t, theta, M, irs_solver_params_.midpoint_iterations, irs_solver_params_.plot_inbetween, 1, irs_solver_params_.min_periastron, irs_solver_params_.initial_guesses, irs_solver_params_.use_ellipse);

		if (b_ != 0) {
			auto [x, y] = OperatorsOrder2::polar_to_cartesian_single(b_, theta, -M_PI / 2);
			double redshift_factor_ = BHphysics::redshift_factor(r, theta, t, M, b_);
			double f_o = BHphysics::flux_observed(r, acc, M, redshift_factor_);
			points.push_back({ x, y, b_, theta, redshift_factor_, f_o });
		}

		if (b_2 != 0) {
			auto [x, y] = OperatorsOrder2::polar_to_cartesian_single(b_2, theta, -M_PI / 2);
			double redshift_factor_2 = BHphysics::redshift_factor(r, theta, t, M, b_2);
			double f_o2 = BHphysics::flux_observed(r, acc, M, redshift_factor_2);
			points2.push_back({ x, y, b_2, theta, redshift_factor_2, f_o2 });
		}
	}

	/*
	This last part is temporary. Only ther for completeness of the Python Code.
	Data is stored in the point,points2 container and can be passed eaasely to other classes/methods.
	*/
	file1 = "D:/dev/Cpp/LuminetCpp/Points/points_incl=85_sampled.csv";
	file2 = "D:/dev/Cpp/LuminetCpp/Points/points_secondary_incl=85_sampled.csv";

	std::ofstream file1OutStream(file1);

	std::string delimitor = ",";//make your choice
	std::string head = "X,Y,impact_parameter, angle,z_factor,flux_o,";//comment out if necessary
	file1OutStream << head << std::endl;//comment out if necessary
	for (const auto& point : points) {
		file1OutStream << point.X << delimitor << point.Y << delimitor << point.impact_parameter << delimitor << point.angle << delimitor
			<< point.z_factor << delimitor << point.flux_o << delimitor << '\n';
	}
	file1OutStream.close();
	std::ofstream file2OutStream(file2);
	file2OutStream << head << std::endl;//comment out if necessary
	for (const auto& point : points2) {
		file2OutStream << point.X << delimitor << point.Y << delimitor << point.impact_parameter << delimitor << point.angle << delimitor
			<< point.z_factor << delimitor << point.flux_o << delimitor << '\n';
	}
	file2OutStream.close();
}