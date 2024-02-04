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
	void plot_isoradials(double, std::vector<double>&, std::vector<double>&, std::vector<double>&, std::vector<double>&, std::vector<double>&, std::vector<double>&,bool =false);//isoradials with redshift
	void plot_iso_redshifts(const double&, const std::multimap<double, std::vector<delaunay::Segment> >&, const double&, const double&, const double&, const double&, const bool&);
	void plot_iso_redshifts(const double&, const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, const double&, const double&, const double&, const double&);//
	//void plot_iso_redshifts(IsoRedShift& Irs, const std::vector<std::vector<double>> irsgrid, const double& inclination, std::multimap<double, std::vector<std::pair<std::vector<double>, std::vector<double> > > >& isolines, const double& x_max_, const double& x_min_, const double& max_rs, const double& min_rs, const bool& loop);
	//void plot_iso_redshifts(const double&, const std::vector<std::vector<double>>&, const std::vector<double>& X, const std::vector<double>& Y, const double&, const double&, const double&, const double&, const bool&);
	//void plot_iso_redshifts(const double&, const std::vector<double>&, const std::vector<double>& , const std::vector<double>& , const double&, const double&, const double&, const double&, const bool&);

private:// Functions to convert a value to RGB format
	std::vector<std::tuple<double, double, double> >convertToRGB(std::vector<double>);
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