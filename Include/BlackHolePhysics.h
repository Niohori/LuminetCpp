#pragma once
#ifndef BLACKHOLEPHYSICS_H
#define BLACKHOLEPHYSICS_H
#include <iostream>
#include <cmath>
#include <vector>
#include <cmath>
#include <functional>
#include <map>
#include <boost/math/special_functions/jacobi_elliptic.hpp>
#include "TensorCalculus.h"
#include "utilities.h"

//using namespace std;
//-----------------------------------------------------------BLACK-HOLE-MATH.PY

class BHphysics
{
public:
	BHphysics();
	static double calc_q(double periastron, double bh_mass, double);
	static double calc_b_from_periastron(double periastron, double bh_mass, double);

	static double k(double periastron, double bh_mass);

	static double k2(double periastron, double bh_mass, double);
	static double zeta_inf(double periastron, double bh_mass, double tol);
	double zeta_r(double periastron, double r, double bh_mass);
	static double cos_gamma(double _a, double incl, double tol);
	void get_plot(const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, const double&);

	double cos_alpha(double phi, double incl);

	double alpha(double phi, double incl);
	std::vector<double> filter_periastrons(const std::vector<double>& periastron, double bh_mass, double tol);

	static double eq13(double periastron, double ir_radius, double ir_angle, double bh_mass, double incl, int n, double tol);
	static std::tuple<std::vector<double>, std::vector<double>, int> midpoint_method(
		const std::function<double(double, double, double, double, double, int, double)> func,
		const std::unordered_map<std::string, double>& args,
		const std::vector<double>& x,
		const std::vector<double>& y,
		int index_of_sign_change);

	static double improve_solutions_midpoint(
		const std::function<double(double, double, double, double, double, int, double)>& func,
		const std::unordered_map<std::string, double>& args,
		const std::vector<double>& x,
		const std::vector<double>& y,
		int index_of_sign_change,
		int iterations
	);

	static double calc_periastron(double _r, double incl, double _alpha, double bh_mass, int midpoint_iterations, bool plot_inbetween, int n, double min_periastron, int initial_guesses);

	static double calc_impact_parameter(double _r, double incl, double _alpha, double bh_mass, int midpoint_iterations, bool plot_inbetween, int n, double min_periastron, int initial_guesses, bool use_ellipse);

	double phi_inf(double periastron, double M);

	double mu(double periastron, double bh_mass);
	static double ellipse(double r, double a, double incl);

	static double flux_intrinsic(double r, double acc, double bh_mass);

	static double flux_observed(double r, double acc, double bh_mass, double redshift_factor);
	static double redshift_factor(double radius, double angle, double incl, double bh_mass, double b_);
	static int find_index_sign_change_indices(const std::vector<double>&);
	//Black body temperature to RGB conversion

static std::vector<double> wavelengthToRGB(const double& , const double& );//artist's impression
////////////////////////////////////////////////////////////////
//
//  Tanner Helland formulas
//
static std::vector<double> convert_TH(const double& temperature, const double& brightness);



////////////////////////////////////////////////////////////////
//
//  Neil Bartlett formulas
//
static std::vector<double> convert_NB(const double& temperature, const double& brightness);
private:
static  void normalizeRGB(std::vector<double>&, const double&);

};
#endif