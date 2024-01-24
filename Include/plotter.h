#pragma once
#ifndef PLOTTER_H
#define PLOTTER_H

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

#include <discpp.h>
const double a_PI = 3.14159265358979323846;

class Plotter {
public:
	Plotter();
	~Plotter();
	void plot(std::vector<double>&, std::vector<double>&, std::vector<double>&, std::vector<double>&);//bare isoradials
	void plot(double, std::vector<double>&, std::vector<double>&, std::vector<double>&, std::vector<double>&, std::vector<double>&, std::vector<double>&);//isoradials with redshift

private:// Functions to convert a value to RGB format
	std::vector<std::tuple<double, double, double> >convertToRGB(std::vector<double>);
	std::vector<double> normalize_vector(std::vector<double>);
private://variables
	Dislin* g;
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