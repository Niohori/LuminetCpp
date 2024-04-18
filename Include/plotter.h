#pragma once
#ifndef PLOTTER_H
#define PLOTTER_H
//#include "winuser.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <type_traits>
#include <sstream>
#include <array>
#include <exception>
#include <stdexcept>
#include <utility>
#include <unordered_map>
#include <algorithm>
#include <map>

#include "IsoRedShift.h"

#include <discpp.h>
#include <Windows.h>
// Undefine macros that may cause conflicts
#undef min
#undef max
const double a_PI = 3.14159265358979323846;

class Plotter {
public:
	Plotter();
	~Plotter();
	void plot_isoradials(std::vector<double>&, std::vector<double>&, std::vector<double>&, std::vector<double>&);//bare isoradials
	void plot_isoradials(double, std::vector<double>&, std::vector<double>&, std::vector<double>&, std::vector<double>&, std::vector<double>&, std::vector<double>&, const double&, bool = false);//isoradials with redshift
	void plot_iso_redshifts(const double&, const double&, const std::multimap<double, std::vector<meshes::Point> >&, const double&, const double&, const double&, const double&, const std::pair<std::vector<double>, std::vector<double>>&, const std::pair<std::vector<double>, std::vector<double>>&, const bool&);
	void plot_BlackHole(double inclination, std::vector<double>& xx_, std::vector<double>&, std::vector<double>&, std::vector<double>&, std::vector<double>&, std::vector<double>&, const bool&);
	void plot_Eq13(const std::vector<double> & peri, const std::vector<double> & val);
private:// Functions to convert a value to RGB format
	std::vector<std::tuple<double, double, double> >convertToRGB(const std::vector<double>&);
	std::vector<std::tuple<double, double, double> >convertToRGBbis(const std::vector<double>&, const double&);
	std::vector<double> normalize_vector(std::vector<double>);
	std::vector<double> flatten_matrix(const std::vector<std::vector<double> >&);
private://variables
	int screenWidth;
	int screenHeight;
	Dislin* g;
	int ic;//dislin
	double x_max;
	double x_min;
	double y_max;
	double y_min;
	int Npoints;
	// Vectors to store RGB color values
	std::vector<std::tuple<double, double, double>> rgbVector;
	std::vector<std::tuple<double, double, double>> rgbVector_g;
};
#endif