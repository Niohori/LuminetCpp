#include "IsoRedShift.h"

IsoRedShift::IsoRedShift() {
	//constructor;
};

IsoRedShift::~IsoRedShift() {
	//destructor;
};

IsoRedShift::IsoRedShift(const double& angle, const double& redshift_, const double& bh_mass, const std::map<double, std::map<int, Isoradial>>) {
	;
}

void IsoRedShift::improve() {
	auto [r_w_s, r_wo_s] = split_co_on_solutions();

	if (!r_w_s.empty()) {  // at least one solution is found
		recalc_isoradials_wo_redshift_solutions(false);
		improve_tip(irs_solver_params_.retry_tip);

		for (int n = 0; n < irs_solver_params_.times_inbetween; ++n) {
			improve_between_all_solutions_once();
			order_coordinates("calculating inbetween", irs_solver_params_.plot_inbetween);
		}
	}
}
void IsoRedShift::recalc_isoradials_wo_redshift_solutions(bool plot_inbetween) {
	auto [r_w_s, r_wo_s] = split_co_on_solutions();

	if (r_wo_s.size() > 0) {
		std::vector<double> a, b;
		for (const auto& r_wo_s_item : r_wo_s) {
			std::tie(a, b) = recalc_redshift_on_closest_isoradial_wo_z(r_wo_s_item.first);
			order_coordinates("improving tip angular", plot_inbetween);
		}
	}
}

void IsoRedShift::improve_tip(int iterations) {
	auto [r_w_so, r_wo_s] = split_co_on_solutions();
	if (r_wo_s.size() > 0) {
		for (int it = 0; it < iterations; ++it) {
			calc_ir_before_closest_ir_wo_z(irs_solver_params_.angular_margin);
			order_coordinates("Improving tip iteration " + std::to_string(it), irs_solver_params_.plot_inbetween);
		}
	}
}



void IsoRedShift::improve_between_all_solutions_once() {

	order_coordinates();
	std::pair<std::map<double, std::vector<std::vector<double>>>, std::map<double, std::vector<std::vector<double>>>>  a_ = split_co_on_solutions();
	std::map<double, std::vector<std::vector<double>>> co_w_s = std::get<0>(a_);//contains for a certain redshift, the polar coordinates with this value
	std::vector<double> redshft = {};
	std::vector<double> angles_ = {};
	std::vector<double> radii_ = {};
	for (std::map<double, std::vector<std::vector<double>>>::iterator it = co_w_s.begin(); it != co_w_s.end(); ++it) {
		redshft.push_back(it->first);
		auto b_ = it->second;
		for (auto polcoor : b_) {
			angles_.push_back(polcoor[0]);
			radii_.push_back(polcoor[0]);
		}
	}

	//std::map<double, std::vector<std::vector<double>>> co_wo_s = std::get<1>(a_);

	for (int i = 0; i < redshft.size() - 1; i++) {
		double co1 = redshft[i];
		double co2 = redshft[i + 1];
		for (int j = 0; j < angles_.size() - 1; j++) {
			double r_inbetw = 0.5 * (radii_[j] + radii_[j + 1]);//midpoint of radius
			double begin_angle = angles_[j] - 0.1;
			double end_angle = angles_[j + 1] + 0.1;
			auto [a, r] = calc_redshift_on_ir_between_angles(r_inbetw, begin_angle, end_angle,
				irs_solver_params_.retry_angular_precision, false, false, "", true);

			if (!a.empty()) {
				add_solutions(a, r, r_inbetw);
			}


		}
	}
}
///////////////////////////////////////////////////////////////////////////////
std::pair<std::vector<double>, std::vector<double>> IsoRedShift::calc_redshift_on_ir_between_angles(
	double radius, double begin_angle, double end_angle, int angular_precision, bool mirror,
	bool plot_inbetween, const std::string& title, bool force_solution) {
	angular_properties _angular_properties_;
	_angular_properties_.start_angle = begin_angle;
	_angular_properties_.end_angle = end_angle;
	_angular_properties_.angular_precision = angular_precision;
	_angular_properties_.mirror = mirror;

	Isoradial ir(radius, t, M, 0, _angular_properties_);

	ir.find_redshift_params_.force_redshift_solution = force_solution;

	auto [a, r] = ir.calc_redshift_location_on_ir(redshift, false);

	if (plot_inbetween) {
		// Assuming you have a plot function in your C++ class
		// You need to implement the plotting functionality in your class
		// This is just a placeholder, and you should adjust it based on your actual implementation
		//*********************************    TO DO ? *******************************************
		//ir.plot_redshift(title);
	}

	return { a, r };

}

void IsoRedShift::calc_ir_before_closest_ir_wo_z(double angular_margin) {
	auto [r_w_s, r_wo_s] = split_co_on_solutions();

	if (r_wo_s.size() > 0 && r_w_s.size() > 0) {
		double first_r_wo_s = r_wo_s.begin()->first;
		double last_r_w_s = r_w_s.rbegin()->first;
		double inbetween_r = 0.5 * (first_r_wo_s + last_r_w_s);
		std::vector<double> angle_interval = radii_w_coordinates_dict[last_r_w_s][0];
		std::vector<double> last_radii = radii_w_coordinates_dict[last_r_w_s][1];

		if (angle_interval.size() > 1) {
			double begin_angle = angle_interval[0];
			double end_angle = angle_interval[1];

			if (end_angle - begin_angle > M_PI) {
				std::swap(begin_angle, end_angle);
			}

			auto [a, r] = calc_redshift_on_ir_between_angles(inbetween_r, begin_angle - angular_margin,
				end_angle + angular_margin,
				irs_solver_params_.retry_angular_precision, false, false, "", false);

			if (!a.empty()) {
				add_solutions(a, r, inbetween_r);
			}
			else {
				radii_w_coordinates_dict[inbetween_r] = { {}, {} };
			}
		}
	}
}

std::pair<std::vector<double>, std::vector<double>> IsoRedShift::recalc_redshift_on_closest_isoradial_wo_z(double angular_margin) {
	
	auto [r_w_s, r_wo_s] = split_co_on_solutions();
	ir_params _ir_params_;
	if (r_wo_s.size() > 0 && r_w_s.size() > 0) {
		double first_r_wo_s = r_wo_s.begin()->first;
		double last_r_w_s = r_w_s.rbegin()->first;
		double inbetween_r = 0.5 * (first_r_wo_s + last_r_w_s);

		std::vector<double> angle_interval = radii_w_coordinates_dict[last_r_w_s][0];
		std::vector<double> last_radii = radii_w_coordinates_dict[last_r_w_s][1];

		if (angle_interval.size() > 1) {
			double begin_angle = angle_interval[0];
			double end_angle = angle_interval[1];

			if (end_angle - begin_angle > M_PI) {
				std::swap(begin_angle, end_angle);
			}

			auto [a, b] = calc_redshift_on_ir_between_angles(inbetween_r, begin_angle - angular_margin,
				end_angle + angular_margin,
				irs_solver_params_.retry_angular_precision, false, false, "", false);

		

			if (!a.empty()) {
				add_solutions(a, b, inbetween_r);
			}
		}
	}
	
	// Placeholder return values
	return { {}, {} };;
}

void IsoRedShift::update() {
	// Implement the update function
}

void IsoRedShift::add_solutions(const std::vector<double>& angles, const std::vector<double>& impact_parameters, double radius_ir) {
	// Implement the add_solutions function
}
/*
std::unordered_map<std::pair<double, double>, double> IsoRedShift::init_co_to_radii_dict() {
	// Implement the init_co_to_radii_dict function
	std::unordered_map<std::pair<double, double>, double> a;
	return a;
}*

std::pair<std::vector<double>, std::vector<double>> IsoRedShift::extract_co_from_solutions_dict() {
	// Implement the extract_co_from_solutions_dict function
	std::pair<std::vector<double>, std::vector<double>> a;
	return a;
}

void IsoRedShift::calc_from_isoradials(const std::vector<Isoradial>& isoradials, bool cartesian) {
	// Implement the calc_from_isoradials function
}

/*
Iterates the dictionary of coordinates that looks like {r_0: [[angle1, angle2], [b_1, b_2]],
		r_1: [[...], [...]]}
		Checks if each key (radius corresponding to an isoradial) has solutions for the isoredshift or not.
		Splits the original dict in two: one with solutions and one without solutions

		:returns: two dictionaries: one with solutions and one without.
*/
std::pair<std::map<double, std::vector<std::vector<double>>>, std::map<double, std::vector<std::vector<double>>>> IsoRedShift::split_co_on_solutions() {
	std::map<double, std::vector<std::vector<double>>> dict_w_s;
	std::map<double, std::vector<std::vector<double>>> dict_wo_s;
		for (const auto& pair : radii_w_coordinates_dict) {
		double key = pair.first;
		const auto& coordinates = pair.second;

		if (coordinates[0].empty()) {
			dict_wo_s[key] = coordinates;
		}
		else {
			dict_w_s[key] = coordinates;
		}
	}
	
	return std::make_pair(dict_w_s, dict_wo_s);
}



std::pair<std::vector<double>, std::vector<double>> IsoRedShift::calc_core_coordinates() {
	// Assuming you have defined Isoradial class and its constructor in your C++ code.
	
Isoradial ir(6. * M, t, M, 0);
return ir.calc_redshift_location_on_ir(redshift);
std::pair<std::vector<double>, std::vector<double>>a;
return a;
}

void IsoRedShift::order_coordinates(const std::string& plot_title , bool plot_inbetween ) {
	/*std::vector<std::pair<double, double>> co;
	for (size_t i = 0; i < angles.size(); ++i) {
		co.emplace_back(angles[i], radii[i]);
	}
	std::pair<std::vector<double>, std::vector<double>> a = OperatorsOrder2::polar_to_cartesian_lists(radii, angles, 0.0);
	x = std::get<0>(a);
	y = std::get<1>(a);
	double cx = std::accumulate(x.begin(), x.end(), 0.0) / x.size();
	double cy = std::accumulate(y.begin(), y.end(), 0.0) / y.size();
	std::vector<double> order_around = { 0.3 * cx, 0.8 * cy };

	std::sort(co.begin(), co.end(), [&](const auto& p1, const auto& p2) {
		return OperatorsOrder2::get_angle_around(order_around, OperatorsOrder2::polar_to_cartesian_single(p1.first, p1.second,0.0)) <
			OperatorsOrder2::get_angle_around(order_around, OperatorsOrder2::polar_to_cartesian_single(p2.first, p2.second,0.0));
		});

	if (plot_inbetween) {
		// Replace this with your actual plotting code in C++
		// or use a suitable third-party library for plotting.
		//plot_coordinates(co, order_around, plot_title);
	}

	angles.clear();
	radii.clear();
	for (const auto& p : co) {
		angles.push_back(p.first);
		radii.push_back(p.second);
	}
	a = OperatorsOrder2::polar_to_cartesian_lists(radii, angles, 0.0);
	x = std::get<0>(a);
	y = std::get<1>(a);*/
}