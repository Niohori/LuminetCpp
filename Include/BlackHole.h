#pragma once
#ifndef BLACKHOLE_H
#define BLACKHOLE_H
#include <iostream>
#include <cmath>
#include <fstream>
#include <random>
#include <algorithm>
#include <unordered_map>

#include "TensorCalculus.h"
#include "BlackHolePhysics.h"
#include "IsoRadials.h"
#include "IsoRedShift.h"
#include "utilities.h"

struct Source {
	double X, Y, impact_parameter, angle, z_factor, flux_o;
};

class BlackHole {
public:
	// Constructor and other methods go here...
	BlackHole();
	BlackHole(double mass = 1.0, double inclination = 80, double acc = 1e-8);
	~BlackHole();
	void sample_Sources(int n_points = 1000, const std::string& f = "", const std::string& f2 = "");
	void calc_isoradials(const std::vector<double>& direct_r, const std::vector<double>& ghost_r);
	void add_isoradial(Isoradial& isoradial, double radius, int order);
	std::map<double, IsoRedShift> calc_isoredshifts(std::vector<double> redshifts = { -0.15, 0.0, 0.1, 0.2, 0.5 });

private://variables
	double inclination;
	double t;
	double M;
	double acc;
	double disk_outer_edge;
	double disk_inner_edge;

	//std::vector<double>  isoradials;
	//std::vector<double> isoredshifts;
	std::map<double, std::map<int, Isoradial>> isoradials;
	std::map<double, IsoRedShift> isoredshifts;

	//make struct from this variables
	//std::unordered_map<std::string, double> settings;
	plot_params plot_params_;
	ir_params ir_parameters_;
	angular_properties angular_properties_;
	irs_solver_params irs_solver_params_;
	solver_params solver_params_;

	/*int initial_guesses;
	int midpoint_iterations;
	bool plot_inbetween;
	bool use_ellipse;*/
	double min_periastron;

	double critical_b;
	int angular_precision;
private://methods
	Isoradial calc_apparent_outer_disk_edge();
	Isoradial calc_apparent_inner_disk_edge();
	double get_apparent_outer_edge_radius(Isoradial&, double angle, double rotation);
	double get_apparent_inner_edge_radius(Isoradial&, double angle, double rotation);
	std::pair<std::vector<double>, std::vector<double>> apparent_inner_edge(Isoradial&, bool cartesian = true, double scale = 0.99);
	std::map<double, std::map<int, Isoradial>> get_dirty_isoradials();
};
#endif