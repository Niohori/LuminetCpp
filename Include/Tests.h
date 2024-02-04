#pragma once
#ifndef TESTS_H
#define TESTS_H

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
#include "BlackHole.h"
#include "BlackHolePhysics.h"
#include "plotter.h"
#include <conio.h>
#include <iostream>
#include <memory>
#include <chrono>

//const double M_PI = 3.14159265358979323846;
class Tests{
	public:
static void test_functions(const double&, const double&);
static void test_BH_rendering(const double&, const double&);
static void test_iso_radials(const double&, const double&, const double&, const unsigned&);
static void test_iso_redshifts(const double&, const double&, const double&, const int &,const unsigned&);
};
#endif