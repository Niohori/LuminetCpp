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
#include <dlib/threads.h>

#include <memory>
#include "mesh.h"
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
		bool isNoise(Segment);
		std::vector<std::vector<Point> > reduceToCurves(std::vector<Segment>&);
		void mergeSegments(std::vector<Segment>&);
		//void buildAdjacencyList(const std::vector<Segment>&);
		std::vector<Point> createPointsList(const std::vector<Segment>&);
		std::vector<Point> reorderPoints(const std::vector<Point>&);
		//std::pair<std::vector<double>, std::vector<double> > smoothCurve(const std::vector<Point>&);
		double distance(const Point&, const Point&);
		//std::vector<std::vector<Point> > constructCurveDFS(std::vector<Segment>& );
		std::vector<std::vector<double> > makeAdjacencyMatrix(const std::vector< Segment >&);
		std::vector<int> dijkstra(const std::vector<std::vector<double>>&, int);
	private://variables
	const double contourValue;
	const std::shared_ptr<Mesh> meshPointer;

	};

#endif