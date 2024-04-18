#include "BlackHolePhysics.h"

BHphysics::BHphysics() {
	;
}

// Find the minimum of the function within the given interval [start, end]
double BHphysics::findMinimum(double ir_radius, double ir_angle, double bh_mass, double incl, int order, double start, double end) {
	double y_0 = eq13(start, ir_radius, ir_angle, bh_mass, incl, order);;
	double y_1= eq13(end, ir_radius, ir_angle, bh_mass, incl, order);;
	if (y_0 * y_1 >= 0) { return NAN; }
	double step = 0.01;
	int nIterations = int((end - start) / step);
	double begin =start;
	double stop = end;
	if (order == 0) {
		stop = start;
		begin = end;
		step *= -1.0;
	}
	int sign1= int(eq13(begin, ir_radius, ir_angle, bh_mass, incl, order)* eq13(begin+step, ir_radius, ir_angle, bh_mass, incl, order)/std::abs(eq13(begin, ir_radius, ir_angle, bh_mass, incl, order) * eq13(begin + step, ir_radius, ir_angle, bh_mass, incl, order)));
	int sign2 =sign1;
	double x = begin;
	int count = 0;
	while (count < nIterations){
		count++;
		y_0 = eq13(x, ir_radius, ir_angle, bh_mass, incl, order);
		y_1 =eq13(x+step, ir_radius, ir_angle, bh_mass, incl, order);
		sign1 = int(y_0 * y_1 / std::abs(y_0*y_1));
		//std::cout << "sign1 vs. sign2: " << sign1 << " <---> " << sign2<< std::endl;
		if (sign1*sign2<0) {//change of sign =sign to stop
			return (2.0 * x + step) / 2.0;
		}
		sign2 = sign1;
		x = x + step;
	}

	return NAN;
}

/*
	  Bisection method implementation.

	  @ -- : parameters function for which the root is to be found.
	  @ start : left limit of root bracketing interval.
	  @ end : right limit of root bracketing interval.

	  @ Returns a double giving the value of the free parameter generating  the root solution.
	*/
double BHphysics::findRoot(double ir_radius, double ir_angle, double bh_mass, double incl, int n, double start, double end) {
	const size_t maxIter = 50;//maximum number of iterations allowed
	const size_t scarboroughCrit = 5; //Scarborough criterion's parameter for tracking convergence of the solution
	double  currmid, prevmid, fmid, checkLeft, checkRight;
	for (int n = 0; n < maxIter; n++)
	{
		currmid = (start + end) / 2;
		fmid = eq13(currmid, ir_radius, ir_angle, bh_mass, incl, n);
		checkLeft = eq13(start, ir_radius, ir_angle, bh_mass, incl, n) * fmid;
		checkRight = eq13(end, ir_radius, ir_angle, bh_mass, incl, n) * fmid;
		if (checkLeft < 0) end = currmid;
		else if (checkRight < 0) start = currmid;
		else if (fmid == 0)
		{
			return currmid;
		}
		else
		{
			return NAN;

		}

		if (n >= 2 && std::abs((currmid - prevmid) / currmid) < 0.5 * pow(10, 2 - static_cast<double>(scarboroughCrit)))
		{
			return currmid;
		}
		prevmid = currmid;
	}
	return currmid;
}

int BHphysics::find_index_sign_change_indices(const std::vector<double>& y) {
	/*
----------------------------------------------------------------------------------------------------------------
//=================================================================================================ed
----------------------------------------------------------------------------------------------------------------
*/
	int index = -1;
	std::vector<int> sign;
	for (unsigned i = 0; i < y.size(); i++) {
		int s = 1;
		if (y[i] < 0.0) { s = -1; }
		sign.push_back(s);
	}
	std::vector<int> diff;
	for (unsigned i = 0; i < sign.size() - 1; i++) {
		diff.push_back(sign[i + 1] - sign[i]);
	}
	for (unsigned i = 0; i < diff.size() - 1; i++) {
		if (diff[i] != 0) { index = i; break; }
	}
	return index;
}

double BHphysics::calc_q(double periastron, double bh_mass) {
	/*
----------------------------------------------------------------------------------------------------------------
 Convert Periastron distance P to the variable Q
	(easier to work with)

	//If periastron - 2. * bh_mass < tol:
	//     // limit for small values
	//     return .5 * (periastron - 2. * bh_mass) * (periastron + 6. * bh_mass)
	// If 1/periastron < tol:
	//     // limit for large values
	//     return periastron
	// If periastron <= 2 * bh_mass:
	//     throw std::invalid_argument("Non-physical periastron found (P <= 2M, aka the photon sphere)."
	//                                 "If you want to calculate non-physical values, you should implement the mpmath library");

	// Q is complex if P < 2M = r_s
----------------------------------------------------------------------------------------------------------------
*/
	return std::sqrt((periastron - 2. * bh_mass) * (periastron + 6. * bh_mass));
}

double BHphysics::calc_b_from_periastron(double periastron, double bh_mass) {
	/*
----------------------------------------------------------------------------------------------------------------
  Get impact parameter b from Periastron distance P
	limits give no substantial speed improvement
	//if (std::abs(periastron) < tol) {  // could physically never happen
	//     std::cout << "tolerance exceeded for calc_b_from_P(P_=" << periastron << ", M=" << bh_mass << ", tol=" << tol << std::endl;
	//     return std::sqrt(3 * periastron * periastron);
	// }
	WARNING: the paper most definitely has a typo here. The fraction on the right-hand side equals b², not b.
	Just fill in u_2 in equation 3, and you'll see. Only this way do the limits P -> 3M and P >> M hold true,
	as well as the value for b_c
----------------------------------------------------------------------------------------------------------------
*/
	return std::sqrt(periastron * periastron * periastron / (periastron - 2. * bh_mass));  // the impact parameter
}

double BHphysics::k(double periastron, double bh_mass) {
	/*
----------------------------------------------------------------------------------------------------------------
 Calculate modulus of elliptic integral
----------------------------------------------------------------------------------------------------------------
*/
	double q = calc_q(periastron, bh_mass);
	// adding limits does not substantially improve speed, nor stability
	// if (q < 10e-3) {  // numerical stability
	//     return std::sqrt(0.5);
	// }
	// WARNING: Paper has an error here. There should be brackets around the numerator.
	return std::sqrt((q - periastron + 6 * bh_mass) / (2 * q));  // the modulus of the elliptic integral
}

double BHphysics::k2(double periastron, double bh_mass) {
	/*
----------------------------------------------------------------------------------------------------------------
Calculate the squared modulus of elliptic integral
----------------------------------------------------------------------------------------------------------------
*/
	double q = calc_q(periastron, bh_mass);
	// adding limits does not substantially improve speed
	// if (1 / periastron <= tol) {
	//     // limit of P -> inf, Q -> P
	//     return 0.;
	// }
	// WARNING: Paper has an error here. There should be brackets around the numerator.
	return (q - periastron + 6 * bh_mass) / (2 * q);  // the modulus of the ellipitic integral
}

double BHphysics::zeta_inf(double periastron, double bh_mass) {
	/*
----------------------------------------------------------------------------------------------------------------
 Calculate Zeta_inf for elliptic integral F(Zeta_inf, k)
----------------------------------------------------------------------------------------------------------------
*/
	double q = calc_q(periastron, bh_mass);
	double arg = (q - periastron + 2 * bh_mass) / (q - periastron + 6 * bh_mass);
	double z_inf = std::asin(std::sqrt(arg));
	return z_inf;
}

double BHphysics::zeta_r(double periastron, double r, double bh_mass) {
	/*
----------------------------------------------------------------------------------------------------------------
Calculate the elliptic integral argument Zeta_r for a given value of P and r
----------------------------------------------------------------------------------------------------------------
*/

	double q = calc_q(periastron, bh_mass);
	double a = (q - periastron + 2 * bh_mass + (4 * bh_mass * periastron) / r) / (q - periastron + (6 * bh_mass));
	double s = std::asin(std::sqrt(a));
	return s;
}

double BHphysics::cos_gamma(double alpha, double incl) {
	/*
----------------------------------------------------------------------------------------------------------------
Calculate the cos of the angle gamma
----------------------------------------------------------------------------------------------------------------
*/
	if (std::abs(incl) < 1e-2) {
		return 0;
	}
	return std::cos(alpha) / std::sqrt(std::cos(alpha) * std::cos(alpha) + 1 / (std::tan(incl) * std::tan(incl)));  // real
}

double BHphysics::cos_alpha(double phi, double incl) {
	/*
----------------------------------------------------------------------------------------------------------------
Returns cos(angle) alpha in observer frame given angles phi (black hole frame) and inclination (black hole frame)
----------------------------------------------------------------------------------------------------------------
*/
	return std::cos(phi) * std::cos(incl) / std::sqrt((1 - std::sin(incl) * std::sin(incl) * std::cos(phi) * std::cos(phi)));
}

double BHphysics::alpha(double phi, double incl) {
	/*
----------------------------------------------------------------------------------------------------------------
Returns observer coordinate of photon given phi (BHF) and inclination (BHF)
----------------------------------------------------------------------------------------------------------------
*/
	return std::acos(cos_alpha(phi, incl));
}

std::vector<double> BHphysics::filter_periastrons(const std::vector<double>& periastron, double bh_mass, double tol) {
	/*
----------------------------------------------------------------------------------------------------------------
	Removes instances where P == 2*M
	returns indices where this was the case
----------------------------------------------------------------------------------------------------------------
*/
	std::vector<double> result;
	for (auto e : periastron) {
		if (std::abs(e - 2. * bh_mass) > tol) {
			result.push_back(e);
		}
	}
	return result;
}
/**
===================================================================================================================================
* @brief Mimics equation (13):Relation between radius (where photon was emitted in accretion disk), a and P.
*@brief	P can be converted to b, yielding the polar coordinates (b, a) on the photographic plate
*@brief This function get called almost everytime when you need to calculate some black hole property
*@param[in] tentative perastion
* @param[in] polar coordinate radius in the black hole equatorial plane(there wher the accrtetion disk resides) of the emmitting particle in the accration disk
* @param[in] ir_angle: the angle alpha (angle between the BH equatorial plane and the photo equatorial plane)
* @param[in] bh_mass: BH mass
* @param[in] incl : the inclination of the observer
*
*
* @return the value  Eq13(P,r) (=0, if periastron is found)
=====================================================================================================================================*/
double BHphysics::eq13(double periastron, double ir_radius, double ir_angle, double bh_mass, double incl, int n) {
	/*
----------------------------------------------------------------------------------------------------------------
	Relation between radius (where photon was emitted in accretion disk), a and P.
	P can be converted to b, yielding the polar coordinates (b, a) on the photographic plate

	This function get called almost everytime when you need to calculate some black hole property

		# Calculate the argument of sn (mod is m = k², same as the original elliptic integral)
	# WARNING: paper has an error here: \sqrt(P / Q) should be in denominator, not numerator
	# There's no way that \gamma and \sqrt(P/Q) can end up on the same side of the division
----------------------------------------------------------------------------------------------------------------
*/

	double z_inf = zeta_inf(periastron, bh_mass);
	double q = calc_q(periastron, bh_mass);
	double m_ = k2(periastron, bh_mass);
	double ell_inf = std::ellint_1(std::sqrt(m_), z_inf);  //incomplete elliptic integral
	double g = std::acos(cos_gamma(ir_angle, incl));

	double ellips_arg;
	if (n) {
		double ell_k = std::ellint_1(std::sqrt(m_), M_PI / 2);//complete elliptic integral
		ellips_arg = (g - 2. * n * M_PI) / (2. * std::sqrt(periastron / q)) - ell_inf + 2. * ell_k;
	}
	else {
		ellips_arg = g / (2. * std::sqrt(periastron / q)) + ell_inf;
	}
	double sn = boost::math::jacobi_sn(std::sqrt(m_), ellips_arg); //Jacobi elliptic function sn TO CHECK:arguments
	double term1 = -(q - periastron + 2. * bh_mass) / (4. * bh_mass * periastron);
	double term2 = ((q - periastron + 6. * bh_mass) / (4. * bh_mass * periastron)) * sn * sn;

	return 1. - ir_radius * (term1 + term2);
}

std::tuple<std::vector<double>, std::vector<double>, int>BHphysics::midpoint_method(
	/*
----------------------------------------------------------------------------------------------------------------
//
----------------------------------------------------------------------------------------------------------------
*/
const std::function<double(double, double, double, double, double, int)> func,
const std::unordered_map<std::string, double>& args,
const std::vector<double>& x,
const std::vector<double>& y,
int index_of_sign_change)
{
	std::vector<double> new_x = x;
	std::vector<double> new_y = y;

	double x0 = new_x[index_of_sign_change];
	double x1 = new_x[index_of_sign_change + 1];
	double inbetween_x = (x0 + x1) / 2.0;
	new_x.insert(new_x.begin() + index_of_sign_change + 1, inbetween_x);

	double y0 = new_y[index_of_sign_change];
	double y1 = new_y[index_of_sign_change + 1];
	double inbetween_solution = func(inbetween_x, args.at("ir_radius"), args.at("ir_angle"), args.at("bh_mass"), args.at("incl"), args.at("n"));
	new_y.insert(new_y.begin() + index_of_sign_change + 1, inbetween_solution);

	int ind_of_sign_change_ = index_of_sign_change + (y0 * inbetween_solution < 0 ? 0 : 1);
	return std::make_tuple(new_x, new_y, ind_of_sign_change_);
}

double BHphysics::improve_solutions_midpoint(
	/*
----------------------------------------------------------------------------------------------------------------
//
----------------------------------------------------------------------------------------------------------------
*/
const std::function<double(double, double, double, double, double, int)>& func,
const std::unordered_map<std::string, double>& args,
const std::vector<double>& x,
const std::vector<double>& y,
int index_of_sign_change,
int iterations
) {
	int index_of_sign_change_ = index_of_sign_change;
	std::vector<double> new_x = x;
	std::vector<double> new_y = y;

	for (int iteration = 0; iteration < iterations; ++iteration) {
		auto [updated_x, updated_y, updated_ind] = BHphysics::midpoint_method(func, args, new_x, new_y, index_of_sign_change_);
		new_x = updated_x;
		new_y = updated_y;
		index_of_sign_change_ = updated_ind;
	}

	return new_x[index_of_sign_change_];
}

double BHphysics::calc_periastron(double _r, double incl, double _alpha, double bh_mass,
	int midpoint_iterations, bool plot_inbetween, int order, double min_periastron, int initial_guesses)
{
	/*
	Given a value for r (BH frame) and alpha (BH/observer frame), calculate the corresponding periastron value
		This periastron can be converted to an impact parameter b, yielding the observer frame coordinates (b, alpha).
		Does this by generating range of periastron values, evaluating eq13 on this range and using a midpoint method
		to iteratively improve which periastron value solves equation 13.
		The considered initial periastron range must not be lower than min_periastron (i.e. the photon sphere),
		otherwise non-physical solutions will be found. These are interesting in their own right (the equation yields
		complex solutions within radii smaller than the photon sphere!), but are for now outside the scope of this project.
		Must be large enough to include solution, hence the dependency on the radius (the bigger the radius of the
		accretion disk where you want to find a solution, the bigger the periastron solution is, generally)

		Args:
			_r (float): radius on the accretion disk (BH frame)
			incl (float): inclination of the black hole
			_alpha: angle along the accretion disk (BH frame and observer frame)
			bh_mass (float): mass of the black hole
			midpoint_iterations (int): amount of midpoint iterations to do when searching a periastron value solving eq13
			plot_inbetween (bool): plot
			*/
	//auto t1 = std::chrono::high_resolution_clock::now();
	std::vector<double> periastron_range = OperatorsOrder2::linspace(min_periastron, 2.0 * _r, initial_guesses);
	std::vector<double> y_values;

	for (double periastron : periastron_range) {
		//auto a = eq13(periastron, _r, _alpha, bh_mass, incl, order, 1e-6);
		//std::cout << "Result of equation 13 = " << a << std::endl;
		y_values.push_back(eq13(periastron, _r, _alpha, bh_mass, incl, order));
	}
	auto ind = find_index_sign_change_indices(y_values);
	//std::cout << "find_index_sign_change_indices = " << ind << std::endl;
	double periastron_solution = (ind >= 0) ? periastron_range[ind] : NAN;

	if (!std::isnan(periastron_solution)) {
		std::unordered_map<std::string, double> args_eq13 = { {"ir_radius", _r},
												  {"ir_angle", _alpha},
												  {"bh_mass", bh_mass},
												  {"incl", incl},
												  {"n", static_cast<double>(order)},
												  {"tol", 1e-8}
		};

		periastron_solution = BHphysics::improve_solutions_midpoint(&BHphysics::eq13, args_eq13, periastron_range, y_values, ind, midpoint_iterations);
	}
	/*auto t2 = std::chrono::high_resolution_clock::now();
	std::cout << "Time (s) for old method  " << (t2 - t1).count() / 1e6 << " ms." << " ("<< std::endl;
	std::cout << std::endl << "=====================================================================================================================================================" << std::endl << std::endl;
	t1 = std::chrono::high_resolution_clock::now();*/
	double lower_bound =3.1 * bh_mass;
	double upper_bound = 2.0*_r;


	 double	minloc = findMinimum(_r, _alpha, bh_mass, incl, order, lower_bound, upper_bound);

	
	//return minval;
	/*t2 = std::chrono::high_resolution_clock::now();
	std::cout << "Time (s) for new method  " << (t2 - t1).count() / 1e6 << " ms." << " (" << std::endl;
	std::cout << std::endl << "=====================================================================================================================================================" << std::endl << std::endl;*/

	// Print optimized results
	//if (order == 1) {
	//std::cout <<    "ORDER " << order <<  "   Radius "<< _r << "  Periastron with old algorithm : " << periastron_solution << "  Dlib lfbgs algorithm gives : " << minloc << " Ratio : "<< minloc/ periastron_solution<<std::endl;// starting_point(0) << std::endl;
	//}
	 //periastron_solution = minloc;
	 //std::cout << "Periastron versus minloc for order "<<order << ": " << periastron_solution << " <---> " << minloc << std::endl;
	return periastron_solution;
}

double BHphysics::calc_impact_parameter(double _r, double incl, double _alpha, double bh_mass, int midpoint_iterations, bool plot_inbetween, int order, double min_periastron, int initial_guesses, bool use_ellipse) {
	/*
----------------------------------------------------------------------------------------------------------------

	Calls the main body: calc_periastron, return the impact parameter.

----------------------------------------------------------------------------------------------------------------
*/
	double periastron_solution = BHphysics::calc_periastron(_r, incl, _alpha, bh_mass, midpoint_iterations, plot_inbetween, order, min_periastron, initial_guesses);
	//std::cout << "periastron_solution = " << periastron_solution<<std::endl;
	//return BHphysics::calc_b_from_periastron(periastron_solution, bh_mass);
	if (periastron_solution == NAN || std::isnan(periastron_solution)) {
		// No periastron was found
		//std::cout << "No solution was found for the periastron." << std::endl;
		//return BHphysics::calc_b_from_periastron(periastron_solution, bh_mass);
		return BHphysics::ellipse(_r, _alpha, incl);
	}
	else if (periastron_solution <= 2. * bh_mass) {
		// Periastron found is non-physical
		if (use_ellipse) {
			//std::cout << "Ellipse used." << std::endl;
			//return BHphysics::calc_b_from_periastron(periastron_solution, bh_mass);
			return BHphysics::ellipse(_r, _alpha, incl);
		}
		else {
			//std::cout << "No physical periastron solution found." << std::endl;
		}
	}
	else {
		// Physical periastron found
		return BHphysics::calc_b_from_periastron(periastron_solution, bh_mass);
	}
}

double BHphysics::phi_inf(double periastron, double M) {
	/*
----------------------------------------------------------------------------------------------------------------
//=========================================================================================================
----------------------------------------------------------------------------------------------------------------
*/
	double q = calc_q(periastron, M);
	double ksq = (q - periastron + 6. * M) / (2. * q);
	double z_inf = zeta_inf(periastron, M);
	double phi = 2. * (std::sqrt(periastron / q)) * (std::ellint_1(std::sqrt(ksq), M_PI / 2.0) - std::ellint_1(std::sqrt(ksq), z_inf));
	return phi;
}

double BHphysics::mu(double periastron, double bh_mass) {
	/*
----------------------------------------------------------------------------------------------------------------
//=====================================================================================================
----------------------------------------------------------------------------------------------------------------
*/
	return 2 * phi_inf(periastron, bh_mass) - M_PI;
}

double BHphysics::ellipse(double r, double a, double incl) {
	/*
----------------------------------------------------------------------------------------------------------------
Equation of an ellipse, reusing the definition of cos_gamma.
	This equation can be used for calculations in the Newtonian limit (large P = b, small a)
	or to visualize the equatorial plane.
----------------------------------------------------------------------------------------------------------------
*/
	double g = std::acos(cos_gamma(a, incl));
	double b_ = r * sin(g);
	return b_;
}

double BHphysics::flux_intrinsic(double r, double acc, double bh_mass) {
	/*
----------------------------------------------------------------------------------------------------------------
//====================================================================================================
----------------------------------------------------------------------------------------------------------------
*/
	double r_ = r / bh_mass;
	double log_arg = ((std::sqrt(r_) + std::sqrt(3)) * (std::sqrt(6) - std::sqrt(3))) / ((std::sqrt(r_) - std::sqrt(3)) * (std::sqrt(6) + std::sqrt(3)));
	double f = (3. * bh_mass * acc / (8 * M_PI)) * (1 / ((r_ - 3) * std::pow(r, 2.5))) *
		(std::sqrt(r_) - std::sqrt(6) + std::pow(3, -0.5) * std::log10(log_arg));
	return f;
}

double BHphysics::flux_observed(double r, double acc, double bh_mass, double redshift_factor) {
	/*
----------------------------------------------------------------------------------------------------------------
//================================================================================
----------------------------------------------------------------------------------------------------------------
*/
	double flux_intr = BHphysics::flux_intrinsic(r, acc, bh_mass);
	return flux_intr / pow(redshift_factor, 4);
}

/**
===================================================================================================================================
* @brief The redshift_factor function calculates the gravitational redshift factor (1 + z), ignoring cosmological redshift.
* @brief Homebrewed.
*
* @param[in] radius The position in polar coordinates of the emitting point
* @param[in] angle The position in polar coordinates of the emitting
* @param[in] incl The inclination of the observer
* @param[in] bh_mass The mass of the black hole
* @param[in] b_ The parameter b in the observers frame
*
*
* @return a double as z_factor
=====================================================================================================================================*/
double BHphysics::redshift_factor(double radius, double angle, double incl, double bh_mass, double b_) {
	/*
----------------------------------------------------------------------------------------------------------------
	WARNING: the paper is absolutely incomprehensible here. Equation 18 for the redshift completely
	 leaves out important factors. It should be:
	1 + z = (1 - Ω*b*cos(η)) * (-g_tt -2Ω*g_tϕ - Ω²*g_ϕϕ)^(-1/2)
	The expressions for the metric components, Ω and the final result of Equation 19 are correct though
	TODO perhaps implement other metrics? e.g. Kerr, where g_tϕ != 0
	gff = (radius * np.sin(incl) * np.sin(angle)) ** 2
	gtt = - (1 - (2. * M) / radius)
----------------------------------------------------------------------------------------------------------------
*/

	double z_factor = (1. + sqrt(bh_mass / pow(radius, 3)) * b_ * sin(incl) * sin(angle)) *
		pow((1 - 3. * bh_mass / radius), -0.5);
	return z_factor;
}

/**
===================================================================================================================================
* @brief The wavelengthToRGB function implements an algorithm to convert a given temperature in Kelvin to corresponding RGB values.
* @brief Homebrewed.
*
* @param[in] temperature The black body temperature in kelvin
* @param[in] brightness The brightness (standard = 100)
*
* @return a rgb vector in the (0.0,1.0) range
=====================================================================================================================================*/
std::vector<double> BHphysics::wavelengthToRGB(const double& temperature, const double& brightness) {//wavelength in Kelvin
	// Normalize temperature to the range [0, 1]
	double t_normalized = static_cast<double>(temperature / 100.0 - 20) / (70 - 20);

	// Linear interpolation between red and blue
	double red = 1.0 - t_normalized;
	double blue = t_normalized;

	// Intensity adjustment
	const double intensityMax = 1.0;
	const double gamma = 0.8;
	red = intensityMax * std::pow(red, gamma);
	blue = intensityMax * std::pow(blue, gamma);

	// Create the RGB color
	std::vector<double> rgb = { red, 0, blue }; // No green component
	//if (temperature >= 4500.0 && temperature <= 5500.0)rgb = { 1.0,1.0,1.0 };
	return rgb;
}

/**
===================================================================================================================================
* @brief The convert_TH function implements an algorithm to convert a given temperature in Kelvin to corresponding RGB values.
* @brief Follows the Tanner Helland formulas.
*
* @param[in] temperature The black body temperature in kelvin
* @param[in] brightness The brightness (standard = 100)
*
* @return a rgb vector in the (0.0,1.0) range
=====================================================================================================================================*/
std::vector<double> BHphysics::convert_TH(const double& temperature, const double& brightness)
{
	std::vector<double> rgb = std::vector < double>(3, 0.0);
	double red = 0.0;
	double green = 0.0;
	double blue = 0.0;
	double t = temperature * 0.01;
	if (t <= 66)
	{
		red = 255;
		green = (99.4708025861 * log(t)) - 161.1195681661;
		if (t > 19)
		{
			blue = (138.5177312231 * log(t - 10)) - 305.0447927307;
		}
		else blue = 0;
	}
	else
	{
		red = 329.698727466 * pow(t - 60, -0.1332047592);
		green = 288.1221695283 * pow(t - 60, -0.0755148492);
		blue = 255;
	}
	rgb = { red,green,blue };
	normalizeRGB(rgb, brightness);
	return rgb;
}

/**
===================================================================================================================================
* @brief The convert_NB function implements an algorithm to convert a given temperature in Kelvin to corresponding RGB values.
* @brief Follows the Neil Bartlett formulas.
*
* @param[in] temperature The black body temperature in kelvin
* @param[in] brightness The brightness (standard = 100)
*
* @return a rgb vector in the (0.0,1.0) range
=====================================================================================================================================*/
std::vector<double> BHphysics::convert_NB(const double& temperature, const double& brightness)
{
	std::vector<double> rgb = std::vector < double>(3, 0.0);
	double red = 0.0;
	double green = 0.0;
	double blue = 0.0;
	double t = temperature * 0.01;
	if (t <= 66)
	{
		red = 255;
		green = t - 2;
		green = -155.25485562709179 - 0.44596950469579133 * green + 104.49216199393888 * log(green);
		blue = 0;
		if (t > 20)
		{
			blue = t - 10;
			blue = -254.76935184120902 + 0.8274096064007395 * blue + 115.67994401066147 * log(blue);
		}
	}
	else
	{
		red = t - 55.0;
		red = 351.97690566805693 + 0.114206453784165 * red - 40.25366309332127 * log(red);
		green = t - 50.0;
		green = 325.4494125711974 + 0.07943456536662342 * green - 28.0852963507957 * log(green);
		blue = 255;
	}
	rgb = { red,green,blue };
	normalizeRGB(rgb, brightness);
	return rgb;
}

void BHphysics::normalizeRGB(std::vector<double>& rgb, const double& brightness)
{
	double f = 0.01 * brightness;
	if (f == 0.0)f = 0.0001;
	//  divide by 255 to get factors between 0..1
	rgb[0] /= (255 / f);
	rgb[1] /= (255 / f);
	rgb[2] /= (255 / f);
}

/**
===================================================================================================================================
* @brief The kelvinToRGB function implements an algorithm to convert a given temperature in Kelvin to corresponding RGB values.
* @brief The algorithm used here is based on the Planckian locus approximation, which maps temperatures to colors in the CIE 1931 color space.
*
* @param[in] temp_kelvin The black body temperature in kelvin
*
* @return a rgb vector in the (0.0,1.0) range
=====================================================================================================================================*/
std::vector<double> BHphysics::kelvinToRGB(double temp_kelvin) {
	std::vector<double> rgb = { 0.0,0.0,0.0 };
	double temperature = temp_kelvin / 100;

	double red, green, blue;

	// Red calculation
	if (temperature <= 66) {
		red = 255;
	}
	else {
		red = temperature - 60;
		red = 329.698727446 * pow(red, -0.1332047592);
		if (red < 0) red = 0;
		if (red > 255) red = 255;
	}

	// Green calculation
	if (temperature <= 66) {
		green = temperature;
		green = 99.4708025861 * log(green) - 161.1195681661;
	}
	else {
		green = temperature - 60;
		green = 288.1221695283 * pow(green, -0.0755148492);
	}
	if (green < 0) green = 0;
	if (green > 255) green = 255;

	// Blue calculation
	if (temperature >= 66) {
		blue = 255;
	}
	else if (temperature <= 19) {
		blue = 0;
	}
	else {
		blue = temperature - 10;
		blue = 138.5177312231 * log(blue) - 305.0447927307;
		if (blue < 0) blue = 0;
		if (blue > 255) blue = 255;
	}

	rgb[0] = red / 255;
	rgb[1] = green / 255;
	rgb[2] = blue / 255;

	return rgb;
}