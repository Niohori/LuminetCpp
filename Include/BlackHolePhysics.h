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
#include <dlib/optimization.h>
#include <dlib/global_optimization.h>

#include "TensorCalculus.h"
#include "utilities.h"
#include "Brent.h"
using namespace dlib;
//-----------------------------------------------------------BLACK-HOLE-MATH.PY
class Eq13Brent;
class BHphysics
{
public:
	BHphysics();
	static double calc_q(double periastron, double bh_mass);
	static double calc_b_from_periastron(double periastron, double bh_mass);

	static double k(double periastron, double bh_mass);

	static double k2(double periastron, double bh_mass);
	static double zeta_inf(double periastron, double bh_mass);
	
	static double zeta_r(double periastron, double r, double bh_mass);
	static double cos_gamma(double _a, double incl);
	//void get_plot(const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, const double&);

	static double cos_alpha(double phi, double incl);

	static double alpha(double phi, double incl);
	static std::vector<double> filter_periastrons(const std::vector<double>& periastron, double bh_mass, double tol);

	static double eq13(double periastron, double ir_radius, double ir_angle, double bh_mass, double incl, int n);

	//solve_g(const double& ir_radius, const double& ir_angle, const double& bh_mass, const double& incl, const int& n, const double& tol = 1e - 6, double initial_guess, int max_iterations);
	static std::tuple<std::vector<double>, std::vector<double>, int> midpoint_method(
		const std::function<double(double, double, double, double, double, int)> func,
		const std::unordered_map<std::string, double>& args,
		const std::vector<double>& x,
		const std::vector<double>& y,
		int index_of_sign_change);

	static double improve_solutions_midpoint(
		const std::function<double(double, double, double, double, double, int)>& func,
		const std::unordered_map<std::string, double>& args,
		const std::vector<double>& x,
		const std::vector<double>& y,
		int index_of_sign_change,
		int iterations
	);

	static double calc_periastron(double _r, double incl, double _alpha, double bh_mass, int midpoint_iterations, bool plot_inbetween, int n, double min_periastron, int initial_guesses);

	static double calc_impact_parameter(double _r, double incl, double _alpha, double bh_mass, int midpoint_iterations, bool plot_inbetween, int n, double min_periastron, int initial_guesses, bool use_ellipse);

	static double phi_inf(double periastron, double M);

	static double mu(double periastron, double bh_mass);
	static double ellipse(double r, double a, double incl);

	static double flux_intrinsic(double r, double acc, double bh_mass);

	static double flux_observed(double r, double acc, double bh_mass, double redshift_factor);
	static double redshift_factor(double radius, double angle, double incl, double bh_mass, double b_);
	static int find_index_sign_change_indices(const std::vector<double>&);
	static double findMinimum(double ir_radius, double ir_angle, double bh_mass, double incl, int n, double start, double end);
	static double findRoot(double ir_radius, double ir_angle, double bh_mass, double incl, int n, double start, double end);



	//Black body temperature to RGB conversion

	static std::vector<double> wavelengthToRGB(const double&, const double&);//artist's impression
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

	static std::vector<double> kelvinToRGB(double temp_kelvin);
	
private:
	static  void normalizeRGB(std::vector<double>&, const double&);
};

#endif
// Define your objective function functor
class Eq13Objective {
public:
	Eq13Objective(double ir_radius, double ir_angle, double bh_mass, double incl, int n) :
		a(ir_radius),
		b(ir_angle),
		c(bh_mass),
		d(incl),
		e(n){}

	double operator()(const matrix<double>& x) const {
		// Define your function y = F(x, a, b, c)
		double y = BHphysics::eq13(x(0), a, b, c, d, e);
		return y*y;
	}

private:
	double a, b, c, d;
	int e;
};

// Define your objective function functor
class Eq13Brent {
public:
	Eq13Brent(double ir_radius, double ir_angle, double bh_mass, double incl, int n) :
		a(ir_radius),
		b(ir_angle),
		c(bh_mass),
		d(incl),
		e(n) {
		;
	}

	double operator()( double x) const {
		// Define your function y = F(x, a, b, c)
		double y = BHphysics::eq13(x, a, b, c, d, e);
		return y * y;
	}

private:
	double a, b, c, d;
	int e;
};