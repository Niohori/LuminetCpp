#pragma once
#ifndef ISORADIALS_H
#define ISORADIALS_H
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
	std::pair<std::vector<double>, std::vector<double>> get_bare_isoradials();
	void calculate_ISCO_region();
	std::pair<std::vector<double>, std::vector<double>> get_ISCO_region();
	std::pair<std::vector<double>, std::vector<double>> get_ISCO_curve();
	std::vector<double> get_redshift_factors();
	void calculate(const bool& isc=false);
	double get_radius();
	double  getMaxX();
	
private://variables
	double M;  // mass of the black hole containing this isoradial
	double rS;//Schwarzschild radius
	double rIsco;// innermost stable circular orbit
	double theta_0;  // inclination of the observer's plane
	double radius;
	int order;
	struct params {
		std::string param = "isoradial_solver_parameters";
	};

	std::vector<double> redshift_factors;//TO DO: pack isoradials and redshift
	std::tuple<std::vector<double>, std::vector<double>> cartesian_co;
	std::pair<std::vector<double>, std::vector<double> > bare_isoradials;//TEMPORARY for debugging: holds the polar coordinates (angles, radii) of the projected isoradial
	std::pair<std::vector<double>, std::vector<double> > ISCO_region;//TEMPORARY for debugging: holds the polar coordinates (angles, radii) of the projected isoradial
	std::pair<std::vector<double>, std::vector<double> > ISCO_boundary;

private://methods
	

public://methods

	std::pair<std::vector<double>, std::vector<double>> calculate_coordinates(const bool& isc=false);
	void calc_redshift_factors();
	bool isPointInsidePolygon(const std::pair<double, double>&, const std::pair<std::vector<double>, std::vector<double> >&);
	std::vector<double> find_angle(double z);
	double get_b_from_angle(double angle);
	void calc_between(int ind);
	
public://variables ? make get method?

	find_redshift_params find_redshift_params_;
	angular_properties angular_properties_;
	solver_params solver_params_;
	plot_params plot_params_;
	ir_params ir_params_;
	double maxCoordinate;
	std::vector<double> X;
	std::vector<double> Y;
	std::vector<double> _radii_b;
	std::vector<double> _angles;
};
#endif