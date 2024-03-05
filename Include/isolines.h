#pragma once
#ifndef ISOLINES_H
#define ISOLINES_H
#include <algorithm>
#include <cmath>
#include <exception>
#include <iostream>
#include <limits>
#include <memory>
#include <utility>
#include <vector>
#include <map>
#include <list>
#include < iomanip >
#include <unordered_map>
#include <map>
#include <stack>
#include <queue>
#include <deque>
#include <chrono>

#include <algorithm>
#include <iostream>
#include <limits>
#include <set>
#include <stdexcept>
#include <vector>
#include<random>//TEMPORARY
#include <dlib/threads.h>

#include <memory>
#include "mesh.h"
#include "utilities.h"
using std::make_shared;


using namespace dlib;
using namespace meshes;
// Define the maximum value for infinity
//const double INF = std::numeric_limits<double>::max();

	class Isolines {
	public:

		Isolines(const std::shared_ptr<meshes::Mesh>&, const double&  );
		
		std::vector<Point>  get_iso_lines();

	public://variables


	private://methods
		std::vector<Segment> generateContourLines();
		std::vector<Point> createPointsList(const std::vector<Segment>&);
		std::vector<Point> reorderPoints(const std::vector<Point>&);
		static double distance(const Point& p1, const Point& p2);


	private://variables
	const double contourValue;
	const std::shared_ptr<Mesh> meshPointer;

	};

#endif