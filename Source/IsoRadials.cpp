#include "IsoRadials.h"

Isoradial::Isoradial() {
	;
}
Isoradial::Isoradial(double radius_, double incl_, double bh_mass_, int order_) :
	/*const std::unordered_map<std::string, double>& params_,
	const std::unordered_map<std::string, bool>& plot_params_,
	const std::unordered_map<std::string, double>& angular_properties_)*/
	M(bh_mass_),
	theta_0(incl_),
	radius(radius_),
	order(order_)
	// params(params_),
	 //plot_params(plot_params_),
	 //angular_properties(angular_properties_)
{
	_radii_b = {};
	_angles = {};
	X = {};
	Y = {};
	cartesian_co = { X, Y };
	redshift_factors = {};
};

Isoradial::Isoradial(double radius_, double incl_, double bh_mass_, int order_, angular_properties  _angular_properties) :
	/*const std::unordered_map<std::string, double>& params_,
	const std::unordered_map<std::string, bool>& plot_params_,
	const std::unordered_map<std::string, double>& angular_properties_)*/
	M(bh_mass_),
	theta_0(incl_),
	radius(radius_),
	order(order_),
	angular_properties_(_angular_properties)
	// params(params_),
	 //plot_params(plot_params_),
	 //angular_properties(angular_properties_)
{
	_radii_b = {};
	_angles = {};
	X = {};
	Y = {};
	cartesian_co = { X, Y };
	redshift_factors = {};
}

std::pair<std::vector<double>, std::vector<double>> Isoradial::get_bare_isoradials() {
	//convert polar to xy

	return OperatorsOrder2::polar_to_cartesian_lists(std::get<1>(bare_isoradials), std::get<0>(bare_isoradials),-M_PI/2);
}
std::vector<double>  Isoradial::get_redshift_factors() {

	return redshift_factors;

}

/*
* **************Necessary?
Isoradial::Isoradial(const std::vector<double>& angles, const std::vector<double>& radius_b)
	: angles(angles), radii_b(radius_b) {
	// Initialize other members if needed
}*/
/*
----------------------------------------------------------------------------------------------------------------
//TO IMPLEMENT ?: for the moment being, the parameters are hardcoded
----------------------------------------------------------------------------------------------------------------
*/
/*std::unordered_map<std::string, double> Isoradial::__read_default_solver_params() {
	// Implement the default solver parameters if needed
	// Return the default parameters as an unordered_map
	std::unordered_map<std::string, double> a;
	return a ;
}*/

/*
----------------------------------------------------------------------------------------------------------------
//TO IMPLEMENT ?: for the moment being, the parameters are hardcoded
----------------------------------------------------------------------------------------------------------------
*/
/*std::unordered_map<std::string, bool> Isoradial::default_plot_params() {
	// Implement the default plot parameters if needed
	// Return the default parameters as an unordered_map
	std::unordered_map<std::string, bool> a;
	return a;
}*/

std::pair<std::vector<double>, std::vector<double>> Isoradial::calculate_coordinates() {
	/*
----------------------------------------------------------------------------------------------------------------
 Calculates the angles (alpha) and radii (b) of the photons emitted at radius self.radius as they would appear
		on the observer's photographic plate. Also saves the corresponding values for the impact parameters (P).

		Args:

		Returns:
			tuple: Tuple containing the angles (alpha) and radii (b) for the image on the observer's photographic plate
----------------------------------------------------------------------------------------------------------------
*/
	double start_angle = angular_properties_.start_angle;
	double end_angle = angular_properties_.end_angle;
	unsigned angular_precision = angular_properties_.angular_precision;
	//angular_precision =10;
	std::vector<double > angles = OperatorsOrder2::linspace(start_angle, end_angle, angular_precision);
	
	std::vector<double> impact_parameters;
	for (auto alpha_ : angles) {
		double b_ = BHphysics::calc_impact_parameter(radius, theta_0, alpha_, M, solver_params_.midpoint_iterations,
			solver_params_.plot_inbetween, order, M * solver_params_.min_periastron, solver_params_.initial_guesses, solver_params_.use_ellipse);
		if (b_ < 1e100) {
			_angles.push_back(alpha_);
			impact_parameters.push_back(b_);
		}
	}

	if (order > 0) {
		// TODO: fix dirty manual flip for ghost images
		std::transform(angles.begin(), angles.end(), angles.begin(), [](double a) { return a + M_PI; });
	}

	// flip image if necessary
	if (theta_0 > M_PI / 2) {
		std::cout <<"Theta gt than 90 degrees" << std::endl;
		std::transform(angles.begin(), angles.end(), angles.begin(), [](double a) { return std::fmod((a + M_PI), (2 * M_PI)); });
	}

	if (angular_properties_.mirror) {
		// by default True. Halves computation time for calculating the full isoradial
		// add second half of image (left half if 0° is set at South)
		size_t s = angles.size();
		for (size_t i = 0; i < s;i++) {
			angles.push_back(2*M_PI-angles[s - i - 1]);
			impact_parameters.push_back(impact_parameters[s-i-1]);
		}
	}

	return std::make_pair(angles, impact_parameters);
}

std::vector<double> Isoradial::calc_redshift_factors() {
	/*
----------------------------------------------------------------------------------------------------------------
Calculates the redshift factor (1 + z) over the line of the isoradial
 ----------------------------------------------------------------------------------------------------------------
 */
 std::vector<double> redshift_factors_;
 std::cout << "redshiftfactors " << std::endl;
 std::cout << "********************************* " << std::endl;
 for (int i = 0; i < _radii_b.size(); ++i) {
	 double redshift = BHphysics::redshift_factor(radius, _angles[i], theta_0, M, _radii_b[i]);
	 //std::cout << i << "i): " << redshift << std::endl;
	 redshift_factors_.push_back(redshift);
 }
 redshift_factors = redshift_factors_;
 return redshift_factors_;

}
/*
----------------------------------------------------------------------------------------------------------------
 ?
----------------------------------------------------------------------------------------------------------------
*/
void Isoradial::calculate() {
	auto result = calculate_coordinates();
	bare_isoradials = { std::get<0>(result) ,std::get<1>(result) };//angle in [0,pi/2]
	_angles = std::get<0>(result);
	_radii_b = std::get<1>(result);
	/*for (int i = 0; i < angles.size(); i++) {
		std::cout << i << "):  angle and radius = (" << angles[i] << ", " << radii_b[i] << ")" << std::endl;
	}*/
	std::vector<double> redshift_factors_=calc_redshift_factors();
}
/*
----------------------------------------------------------------------------------------------------------------
Returns angle at which the isoradial redshift equals some value z
		Args:
			z: The redshift value z. Do not confuse with redshift factor 1 + z
----------------------------------------------------------------------------------------------------------------
*/
std::vector<double> Isoradial::find_angle(double z) {
	/*std::vector<int> indices;
	for (int i = 0; i < redshift_factors.size() - 1; ++i) {
		if ((redshift_factors[i] - z - 1) * (redshift_factors[i + 1] - z - 1) < 0) {
			indices.push_back(i + 1);
		}
	}

	std::vector<double> result;
	for (int index : indices) {
		result.push_back(angles[index]);
	}
	return result;*/
	std::vector<double> a;
	return a;
}

/*
----------------------------------------------------------------------------------------------------------------
# TODO: this method only works if angles augment from index 0 to end
		# if image is flipped, then the mod operator makes it so they jump back to 0 about halfway
		# yielding a fake intersection
----------------------------------------------------------------------------------------------------------------
*/
double Isoradial::get_b_from_angle(double angle) {
	// TODO: this method only works if angles augment from index 0 to end
	// if the image is flipped, then the mod operator makes it so they jump back to 0 about halfway
	// yielding a fake intersection
	/*std::transform(angles.begin(), angles.end(), std::back_inserter(d),
		[angle](double a) { return std::abs(fmod(a, 2 * M_PI) - fmod(angle, 2 * M_PI)); });

	auto minIt = std::min_element(d.begin(), d.end());
	if (minIt != d.end()) {
		int resIndex = std::distance(d.begin(), minIt);
		return radii_b[resIndex];
	}
	else {
		return -1e100;
	}*/
	return 0.;

	/*
	----------------------------------------------------------------------------------------------------------------
	Calculates the impact parameter and redshift factor at the
			isoradial angle between place ind and ind + 1

			Args:
				ind: the index denoting the location at which the middle point should be calculated. The impact parameter,
				redshift factor, b (observer plane) and alpha (observer/BH coordinate system) will be calculated on the
				isoradial between location ind and ind + 1

			Returns:
				None: Nothing. Updates the isoradial.
	----------------------------------------------------------------------------------------------------------------
	*/
	//void Isoradial::calc_between(int ind) {
		/*double mid_angle = 0.5 * (angles[ind] + angles[ind + 1]);
		double b_ = BHphysics::calc_impact_parameter(radius, t, mid_angle, M, solver_params_.midpoint_iterations, solver_params_.plot_inbetween, order, M*solver_params_.min_periastron, solver_params_.initial_guesses, solver_params_.use_ellipse);
		double z_ = BHphysics::redshift_factor(radius, mid_angle, t, M, b_);

		radii_b.insert(radii_b.begin() + ind + 1, b_);
		angles.insert(angles.begin() + ind + 1, mid_angle);
		redshift_factors.insert(redshift_factors.begin() + ind + 1, z_);
		*/
}//

/*
----------------------------------------------------------------------------------------------------------------
# TODO: improve this method, currently does not seem to work

		If you know a redshift should exist on the isoradial, use this function to calculate the isoradial until
		it finds it. Useful for when the redshift you're looking for equals (or is close to) the maximum
		redshift along some isoradial line.

		Only works if the redshift can be found within the isoradial begin and end angle.
----------------------------------------------------------------------------------------------------------------
*/

std::vector<double> Isoradial::force_intersection(double redshift) {
	/*if (angles.size() == 2) {
		calc_between(0);
	}

	std::vector<double> diff;
	for (double z_ : redshift_factors) {
		diff.push_back(redshift + 1 - z_);
	}

	std::vector<int> cross;
	for (int i = 0; i < diff.size() - 1; ++i) {
		if (diff[i] * diff[i + 1] < 0) {
			cross.push_back(i);
		}
	}

	if (!cross.empty()) {
		return diff;  // intersection is found
	}

	int it = 0;
	while (cross.empty() && it < find_redshift_params_.max_force_iter) {
		// calc derivatives
		std::vector<double> delta;
		for (int i = 0; i < redshift_factors.size() - 1; ++i) {
			delta.push_back(redshift_factors[i + 1] - redshift_factors[i]);
		}

		// where does the redshift go back up/down before it reaches the redshift we want to find
		std::vector<int> initial_guess_indices;
		for (int i = 0; i < delta.size() - 1; ++i) {
			if (delta[i] * delta[i + 1] < 0) {
				initial_guess_indices.push_back(i);
			}
		}

		int new_ind = initial_guess_indices[0];  // initialize the initial guess.
		calc_between(new_ind);  // insert a more accurate solution

		// calc new interval
		diff.clear();
		for (double z_ : redshift_factors) {
			diff.push_back(redshift + 1 - z_);
		}

		// where does the redshift go back up/down before it reaches the redshift we want to find
		cross.clear();
		for (int i = 0; i < diff.size() - 1; ++i) {
			if (diff[i] * diff[i + 1] < 0) {
				cross.push_back(i);
			}
		}

		it += 1;
	}

	return diff;*/
	std::vector<double> a;
	return a;
}

/*
----------------------------------------------------------------------------------------------------------------
Calculates which location on the isoradial has some redshift value (not redshift factor)
		Does this by means of a midpoint method, with midpoint_steps steps (defined in parameters.ini).
		The (b, alpha, z) coordinates of the isoradial are calculated closer and closer to the desired z.
		It does not matter all that much how high the isoradial resolution is, since midpoint_steps is
		much more important to find an accurate location.
----------------------------------------------------------------------------------------------------------------
*/
std::pair<std::vector<double>, std::vector<double>> Isoradial::calc_redshift_location_on_ir(double redshift, bool cartesian) {
	/*/std::vector<double> diff;
	for (double z_ : redshift_factors) {
		diff.push_back(redshift + 1 - z_);
	}

	// if (find_redshift_params["force_redshift_solution"]) {
	//     // TODO: force_intersection does not always seem to work
	//     diff = force_intersection(redshift);
	// }

	std::vector<int> initial_guess_indices;
	for (int i = 0; i < diff.size() - 1; ++i) {
		if (diff[i] * diff[i + 1] < 0) {
			initial_guess_indices.push_back(i);
		}
	}

	std::vector<double> angle_solutions;
	std::vector<double> b_solutions;

	if (!initial_guess_indices.empty()) {
		for (int s = 0; s < initial_guess_indices.size(); ++s) {
			int new_ind = initial_guess_indices[s];

			for (int i = 0; i < solver_params_.midpoint_iterations; ++i) {
				calc_between(new_ind);  // insert a more accurate solution
				std::vector<double> diff_;
				for (int i = new_ind; i < new_ind + 3; ++i) {
					diff_.push_back(redshift + 1 - redshift_factors[i]);
				}

				// calc new interval
				std::vector<int> start;
				for (int i = 0; i < diff_.size() - 1; ++i) {
					if (diff_[i] * diff_[i + 1] < 0) {
						start.push_back(i);
					}
				}

				new_ind += start[0];  // index of new redshift solution in refined isoradial
			}

			// append average values of the final interval
			angle_solutions.push_back(0.5 * (angles[new_ind] + angles[new_ind + 1]));
			b_solutions.push_back(0.5 * (radii_b[new_ind] + radii_b[new_ind + 1]));

			// update the initial guess indices, as the indexing has changed due to inserted solutions
			for (int& e : initial_guess_indices) {
				e += solver_params_.midpoint_iterations;
			}
		}
	}

	if (cartesian) {
		return OperatorsOrder2::polar_to_cartesian_lists(b_solutions, angle_solutions,0.0);
	}

	return std::make_pair(angle_solutions, b_solutions);*/
	std::pair<std::vector<double>, std::vector<double>> a;
	return a;
}