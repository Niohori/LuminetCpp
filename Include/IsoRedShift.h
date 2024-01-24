#pragma once
#include <iostream>
#include <cmath>
#include <fstream>
#include <random>
#include <algorithm>
#include <unordered_map>
#include <iostream>
#include <cmath>
#include <vector>
#include <cmath>
#include <functional>
#include <map>
#include <numeric>

#include "TensorCalculus.h"
#include "BlackHolePhysics.h"
#include "IsoRadials.h"
#include "utilities.h"

class IsoRedShift {
private://variables
	double theta_0;  // Inclination
	double redshift;
	double M;  // Black hole mass
	std::unordered_map<double, std::vector<std::vector<double>>> radii_w_coordinates_dict;
	//*****************************************************************************************
	//                                             v
	//std::unordered_map<std::pair<double, double>, double> coordinates_with_radii_dict;//This give a problem
	// I have to rewrite this for the methods
	// 	   std::unordered_map<std::pair<double, double>, double> init_co_to_radii_dict();
	// and std::pair<std::vector<double>, std::vector<double>> extract_co_from_solutions_dict();
	//                                             o
	//*****************************************************************************************

	std::vector<double> angles;
	std::vector<double> radii;
	std::vector<double> x;
	std::vector<double> y;

	//std::map<double, std::vector<double>> radii_w_coordinates_dict;
	//std::vector<double> ir_radii_w_co;
	//std::pair<std::vector<double>, std::vector<double>> co;
	double max_radius;
	irs_solver_params irs_solver_params_;
	find_redshift_params find_redshift_params_;
	angular_properties angular_properties_;
	solver_params solver_params_;
	plot_params plot_params_;
	ir_params ir_params_;

private://methods
	//void update();
	/*void add_solutions(const std::vector<double>& angles, const std::vector<double>& impact_parameters, double radius_ir);
	std::unordered_map<std::pair<double, double>, double> init_co_to_radii_dict();
	std::pair<std::vector<double>, std::vector<double>> extract_co_from_solutions_dict();
	void calc_from_isoradials(const std::vector<Isoradial>& isoradials, bool cartesian = false);
	std::pair<std::map<double, std::vector<std::vector<double>>>, std::map<double, std::vector<std::vector<double>>>> split_co_on_solutions();
	std::pair<std::vector<double>, std::vector<double>> calc_core_coordinates();
	void order_coordinates(const std::string& plot_title = "", bool plot_inbetween = false);
	void calc_ir_before_closest_ir_wo_z(double angular_margin);
	std::pair<std::vector<double>, std::vector<double>> calc_redshift_on_ir_between_angles(
		double radius, double begin_angle, double end_angle, int angular_precision, bool mirror,
		bool plot_inbetween, const std::string& title, bool force_solution);*/

public://methods
	IsoRedShift();
	IsoRedShift(const double&, const double&, const double&, const std::map<double, std::map<int, Isoradial>>);
	~IsoRedShift();
	void improve();
	void update();

	void add_solutions(const std::vector<double>& angles, const std::vector<double>& impact_parameters, double radius_ir);
	//std::unordered_map<std::pair<double, double>, double> init_co_to_radii_dict();
	std::pair<std::vector<double>, std::vector<double>> extract_co_from_solutions_dict();
	void calc_from_isoradials(const std::vector<Isoradial>& isoradials, bool cartesian = false);

	std::pair<std::vector<double>, std::vector<double>> calc_core_coordinates();
	void order_coordinates(const std::string& plot_title = "", bool plot_inbetween = false);
	void calc_ir_before_closest_ir_wo_z(double angular_margin);
	void recalc_isoradials_wo_redshift_solutions(bool plot_inbetween);
	void improve_tip(int iterations);
	void improve_between_all_solutions_once();
	std::pair<std::vector<double>, std::vector<double>> calc_redshift_on_ir_between_angles(
		double radius, double begin_angle, double end_angle, int angular_precision, bool mirror,
		bool plot_inbetween, const std::string& title, bool force_solution);
	std::pair<std::vector<double>, std::vector<double>> recalc_redshift_on_closest_isoradial_wo_z(double angular_margin);
	std::pair<std::map<double, std::vector<std::vector<double>>>, std::map<double, std::vector<std::vector<double>>>> split_co_on_solutions();

	std::vector<double> get_ir_radii_w_co();
public://variables ? make get method?
};