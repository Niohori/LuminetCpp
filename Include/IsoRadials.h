#pragma once
/*#include <iostream>
#include <cmath>
#include <fstream>
#include <random>
#include <algorithm>
#include <unordered_map>*/

#include "TensorCalculus.h"
#include "BlackHolePhysics.h"
#include "utilities.h"

class Isoradial {
public:
	Isoradial();
	Isoradial(double radius, double incl, double bh_mass, int order = 0);
	Isoradial(double radius, double incl, double bh_mass, int order, angular_properties);
	//Isoradial(const std::vector<double>& angles, const std::vector<double>& radius_b);

private://variables
	double M;  // mass of the black hole containing this isoradial
	double t;  // inclination of the observer's plane
	double radius;
	int order;
	struct params {
		std::string param = "isoradial_solver_parameters";
	};

	
	std::vector<double> redshift_factors;
	std::tuple<std::vector<double>, std::vector<double>> cartesian_co;

private://methods

public://methods

	std::pair<std::vector<double>, std::vector<double>> calculate_coordinates();
	std::vector<double> calc_redshift_factors();
	void calculate();
	std::vector<double> find_angle(double z);
	double get_b_from_angle(double angle);
	void calc_between(int ind);
	std::vector<double> force_intersection(double redshift);
	std::pair<std::vector<double>, std::vector<double>> calc_redshift_location_on_ir(double redshift, bool cartesian = false);

public://variables ? make get method?

	find_redshift_params find_redshift_params_;
	angular_properties angular_properties_;
	solver_params solver_params_;
	plot_params plot_params_;
	ir_params ir_params_;

	std::vector<double> X;
	std::vector<double> Y;
	std::vector<double> radii_b;
	std::vector<double> angles;
};