#include "IsoRedShift.h"

IsoRedShift::IsoRedShift() {
	//constructor;
};

IsoRedShift::~IsoRedShift() {
	//destructor;
};

IsoRedShift::IsoRedShift(const double& angle, const double& redshift_, const double& bh_mass, const std::map<double, std::pair<int, Isoradial*> > & from_isoradials)
	: theta_0(angle), M(bh_mass), redshift(redshift_) {
	if (from_isoradials.empty()) {
		// TODO: initialize from photon sphere isoradial?
	}
	else {
		calc_from_isoradials(from_isoradials,false);
	}

	coordinates_with_radii_dict = init_co_to_radii_dict();
	 get_ir_radii_with_co();
	std::tie(angles, radii) = extract_co_from_solutions_dict();
	max_radius = radii.empty() ? 0 : *std::max_element(radii.begin(), radii.end());
	std::tie(x, y) = OperatorsOrder2::polar_to_cartesian_lists(radii, angles, 0); // Assuming polar_to_cartesian_lists is a defined function
	order_coordinates();
}

void IsoRedShift::get_ir_radii_with_co() {
	ir_radii_w_co.clear();
	for (const auto& entry : radii_w_coordinates_dict) {
		if (!entry.second.first.empty()) {  // if radius has solution
			ir_radii_w_co.push_back(entry.first);
		}
	}
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
	std::pair<std::map<double, std::pair<std::vector<double>, std::vector<double>>>, std::map<double, std::pair<std::vector<double>, std::vector<double>>> >  a_ = split_co_on_solutions();
	std::map<double, std::pair<std::vector<double>, std::vector<double> > > co_w_s = std::get<0>(a_);//contains for a certain redshift, the polar coordinates with this value
	std::vector<double> redshft = {};
	std::vector<double> angles_ = {};
	std::vector<double> radii_ = {};
	for (std::map<double, std::pair<std::vector<double>, std::vector<double>> >::iterator it = co_w_s.begin(); it != co_w_s.end(); ++it) {
		redshft.push_back(it->first);
		auto b_ = it->second;
		for (size_t i = 0; i < b_.first.size(); i++) {
			angles_.push_back(b_.first[i]);
			radii_.push_back(b_.second[i]);
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

	Isoradial ir(radius, theta_0, M, 0, _angular_properties_);

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
		std::vector<double> angle_interval = radii_w_coordinates_dict[last_r_w_s].first;
		std::vector<double> last_radii = radii_w_coordinates_dict[last_r_w_s].second;

		if (angle_interval.size() > 1) {
			double begin_angle = angle_interval[0];
			double end_angle = angle_interval[1];

			if (end_angle - begin_angle > M_PI) {
				std::swap(begin_angle, end_angle);
			}

			auto [a, r] = calc_redshift_on_ir_between_angles(inbetween_r, begin_angle - angular_margin,
				end_angle + angular_margin,
				irs_solver_params_.retry_angular_precision, false, false, "", false);
			if (a.size() == 1) {//safegard if no redshift was found for this isoradial
				if (a[0] == 0.0 && r[0] == 0.0) {
					a.clear();
					r.clear();
				}
			}

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

		std::vector<double> angle_interval = radii_w_coordinates_dict[last_r_w_s].first;
		std::vector<double> last_radii = radii_w_coordinates_dict[last_r_w_s].second;

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

/*std::unordered_map<std::pair<double, double>, double> IsoRedShift::init_co_to_radii_dict() {
	// Implement the init_co_to_radii_dict function
	std::unordered_map<std::pair<double, double>, double> a;
	return a;

	*************************  ???? USELESS ***********************************************
}*/

std::unordered_map< double, std::pair<std::vector<double>, std::vector<double> > > IsoRedShift::init_co_to_radii_dict() {
	std::unordered_map< double, std::pair<std::vector<double>, std::vector<double> > > to_return;

	/*for (const auto& [radius, co] : radii_w_coordinates_dict) {
		if (!co.first.empty()) {  // if radius has solution
			std::vector<std::pair<double, double>> coordinates;
			for (size_t i = 0; i < co.first.size(); ++i) {
				coordinates.emplace_back(co.first[i], co.second[i]);
			}
			to_return[radius] = coordinates;
			for (const auto& co_ : coordinates) {  // either one or two solutions
				//to_return[co_] = radius;
				to_return[radius] = co_;
			}
		}
	}*/
	to_return = radii_w_coordinates_dict;
	return to_return;
}

std::vector<double> IsoRedShift::get_ir_radii_w_co() {
	std::vector<double> result;
	for (const auto& entry : radii_w_coordinates_dict) {
		if (entry.second.first.size() > 0) {
			result.push_back(entry.first);
		}
	}
	return result;
}

/*
Iterates the dictionary of coordinates that looks like {r_0: [[angle1, angle2], [b_1, b_2]],
		r_1: [[...], [...]]}
		Checks if each key (radius corresponding to an isoradial) has solutions for the isoredshift or not.
		Splits the original dict in two: one with solutions and one without solutions

		:returns: two dictionaries: one with solutions and one without.
*/
std::pair<std::map<double, std::pair<std::vector<double>, std::vector<double>>>, std::map<double, std::pair<std::vector<double>, std::vector<double>>> > IsoRedShift::split_co_on_solutions() {
	std::map<double, std::pair<std::vector<double>, std::vector<double>>> dict_w_s;
	std::map<double, std::pair<std::vector<double>, std::vector<double>>> dict_wo_s;
	for (const auto& pair : radii_w_coordinates_dict) {
		double key = pair.first;
		const auto& coordinates = pair.second;

		if (coordinates.first.empty()) {
			dict_wo_s[key] = coordinates;
		}
		else {
			dict_w_s[key] = coordinates;
		}
	}

	return std::make_pair(dict_w_s, dict_wo_s);
}

std::pair<std::vector<double>, std::vector<double>> IsoRedShift::calc_core_coordinates() {
	//

	Isoradial ir(6. * M, theta_0, M, 0);
	return ir.calc_redshift_location_on_ir(redshift);
	std::pair<std::vector<double>, std::vector<double>>a;
	return a;
}

void IsoRedShift::order_coordinates(const std::string& plot_title, bool plot_inbetween) {
	std::vector<std::pair<double, double>> co;
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
		std::vector<double> theta_1 = { OperatorsOrder2::polar_to_cartesian_single(p1.first, p1.second, 0.0).first,OperatorsOrder2::polar_to_cartesian_single(p1.first, p1.second, 0.0).second };
		std::vector<double> theta_2 = { OperatorsOrder2::polar_to_cartesian_single(p2.first, p2.second, 0.0).first,OperatorsOrder2::polar_to_cartesian_single(p2.first, p2.second, 0.0).second };
		return OperatorsOrder2::get_angle_around(order_around, theta_1) <
			OperatorsOrder2::get_angle_around(order_around, theta_2);
		});

	if (plot_inbetween) {
		// Use Dislin?
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
	y = std::get<1>(a);
}

std::pair<std::vector<double>, std::vector<double> > IsoRedShift::extract_co_from_solutions_dict() {
	std::vector<double> ang;
	std::vector<double> rad;
	//std::unordered_map<double, std::pair<std::vector<double>,std::vector<double>>> radii_w_coordinates_dict_;
	for (const auto& entry : radii_w_coordinates_dict) {
		const auto key = entry.first;
		const auto val = entry.second;
		//const auto& val = entry.first;
		if (val.first.size() > 0) {  // at least one solution was found
			const auto& angles = val.first;
			const auto& radii = val.second;

			for (const auto& angle : angles) {
				ang.push_back(angle);
			}

			for (const auto& radius : radii) {
				rad.push_back(radius);
			}
		}
	}
	angles = ang;
	radii = rad;
	return std::make_pair(ang, rad);
}

void IsoRedShift::calc_from_isoradials(std::map<double, std::pair<int, Isoradial*> > isoradials, bool cartesian) {
	std::unordered_map<double, std::pair<std::vector<double>, std::vector<double>> > solutions;
	double _max_radius = 0.0;
	for (auto ir_ : isoradials) {
		std::pair<std::vector<double>, std::vector<double>> coor = ir_.second.second->calc_redshift_location_on_ir(redshift, cartesian);
		solutions[ir_.second.second->get_radius()] = coor;

		// Update _max_radius if needed
		if (ir_.second.second->get_radius() > _max_radius) {
			_max_radius = ir_.second.second->get_radius();
		}
	}
	radii_w_coordinates_dict = solutions;
	update();  
}

void IsoRedShift::update() {
	std::vector<double> ir_radii_w_co = {};
	for (const auto& [key, val] : radii_w_coordinates_dict) {
		if (!val.first.empty()) {
			ir_radii_w_co.push_back(key);
		}
	}
	std::pair<std::vector<double>, std::vector<double>> co = extract_co_from_solutions_dict();
	std::tie(x, y) = OperatorsOrder2::polar_to_cartesian_lists(radii, angles, 0); 
	order_coordinates();
}

void IsoRedShift::add_solution(double angle, double radius_b, double radius_ir) {
	auto it = radii_w_coordinates_dict.find(radius_ir);

	if (it != radii_w_coordinates_dict.end()) {
		if (!it->second.first.empty()) {
			it->second.first.push_back(angle);
			it->second.second.push_back(radius_b);
		}
		else {
			it->second = { std::vector<double>{angle}, std::vector<double>{radius_b} };
		}
	}
	else {
		radii_w_coordinates_dict[radius_ir] = { std::vector<double>{angle}, std::vector<double>{radius_b} };
	}

	//coordinates_with_radii_dict.emplace_back(angle, radius_b);
	update();
}

void IsoRedShift::add_solutions(const std::vector<double>& angles, const std::vector<double>& impact_parameters, double radius_ir) {
	for (size_t i = 0; i < angles.size() && i < impact_parameters.size(); ++i) {
		add_solution(angles[i], impact_parameters[i], radius_ir);
	}
}