#pragma once
#ifndef TENSORCALCULUS_H
#define TENSORCALCULUS_H

/*
General library in progress with common operators
To extend with tensor algebra operations.
*/
#include <iostream>
#include <vector>
#include <tuple>
#include <functional>
#include <cmath>

#include "utilities.h"
//#include <Eigen/Dense>

//using namespace Eigen;

//const double M_PI = 3.14159265358979323846;
class OperatorsOrder2 {
public:
	OperatorsOrder2(int nGrid, double delta);
	~OperatorsOrder2();

	std::vector<std::vector<std::vector<double>>> laplace(const std::vector<std::vector<std::vector<double>>>& fct);
	std::tuple<std::vector<std::vector<std::vector<double>>>, std::vector<std::vector<std::vector<double>>>, std::vector<std::vector<std::vector<double>>>> gradient(const std::vector<std::vector<std::vector<double>>>& fct);
	std::vector<std::vector<std::vector<double>>> divergence(const std::vector<std::vector<std::vector<double>>>& v_x,
		const std::vector<std::vector<std::vector<double>>>& v_y,
		const std::vector<std::vector<std::vector<double>>>& v_z);
	std::tuple<std::vector<std::vector<std::vector<double>>>, std::vector<std::vector<std::vector<double>>>, std::vector<std::vector<std::vector<double>>>> gradDiv(const std::vector<std::vector<std::vector<double>>>& A_x,
		const std::vector<std::vector<std::vector<double>>>& A_y,
		const std::vector<std::vector<std::vector<double>>>& A_z);
	std::tuple<double, double, double> partialDerivs(const std::vector<std::vector<std::vector<double>>>& fct, int i, int j, int k);
	std::vector<std::vector<std::vector<double>>> Add(const std::vector<std::vector<std::vector<double>>>&, const double&, const std::vector<std::vector<std::vector<double>>>&, const double&);
	static std::vector<double>  linspace(const double&, const double&, const int&);
	static std::vector<double> nonLinSpace(const bool&, const double&, const int&);
	static double phi(const double&, const double&);
	static std::vector<double>  logspace(const double&, const double&, const int&);
	//static std::vector<double> ellipticspace(const bool&, const double& , const int& );

	static std::pair<std::vector<double>, std::vector<double>> polar_to_cartesian_lists(const std::vector<double>& radii, const std::vector<double>& angles, const double& rotation);
	static std::pair<double, double> polar_to_cartesian_single(double th, double radius, double rotation);
	//static std::vector<double> polar_to_cartesian_single_as_vector(double th, double radius, double rotation);
	static std::pair<double, double> cartesian_to_polar(double x, double y);
	static double get_angle_around(const std::vector<double>& p1, const std::vector<double>& p2);

private:
	static std::vector<double> matrixMultiply(const std::vector<std::vector<double>>& matrix, const std::vector<double>& vector);
	int nGrid;
	double delta;
};
#endif
