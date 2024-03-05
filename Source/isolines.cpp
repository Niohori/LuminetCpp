#include "isolines.h"

/**
===================================================================================================================================
* @brief Class constructor (1 object per desired isovalue) .

*
* @param[in] std::shared_ptr<meshes::Mesh> : a pointer to the Delauney mesh generated in the Mesh class
* @param[in] double  : the desired isovalue
*
* @return None: instantiates the necessary class members
=====================================================================================================================================*/

Isolines::Isolines(const std::shared_ptr<meshes::Mesh>& meshptr, const double& isovalue) :
	meshPointer(meshptr),
	contourValue(isovalue)
{
	;
}

/**
===================================================================================================================================
* @brief Public function called to generate isolines.

*
* @param[in] None (uses class members)
*
* @return a vector of "Point" struct representing the isolines
=====================================================================================================================================*/

std::vector<meshes::Point> Isolines::get_iso_lines() {
	//auto t1 = std::chrono::high_resolution_clock::now();
	std::vector<meshes::Segment>all_isolines;
	size_t nt = 0;
	double average_area = 0.0;
	double max_area = -10000000.0;
	double min_area = 10000000.0;
	all_isolines = generateContourLines();
	std::vector<Point> points = createPointsList(all_isolines);
	return points;
}

/**
===================================================================================================================================
* @brief Private core function to generate isolines.

*
* @param[in] None (usesclass members)
*
* @return a vector of "segments" struct i.e. unordered segments representing the isolines
=====================================================================================================================================*/
std::vector<meshes::Segment> Isolines::generateContourLines() {
	//auto t1 = std::chrono::high_resolution_clock::now();
	//std::vector<Point> Isolines::generateContourLines(const std::vector<driehoek>&driehoeken, double contourValue) {
	std::vector<meshes::Segment> contourSegments;

	for (const driehoek& driehk : meshPointer->driehoeken) {
		double z1 = driehk.Rs1;
		double z2 = driehk.Rs2;
		double z3 = driehk.Rs3;
		meshes::Point p1 = driehk.p1;
		meshes::Point p2 = driehk.p2;
		meshes::Point p3 = driehk.p3;
		auto ID = driehk.ID;


		std::vector<Point> intersections;
		// Check each edge of the triangle
		if ((z1 - contourValue) * (z2 - contourValue) < 0) {
			double t = (contourValue - z1) / (z2 - z1);
			meshes::Point intersectionPoint(p1.x + t * (p2.x - p1.x), p1.y + t * (p2.y - p1.y));
			
		}
		if ((z2 - contourValue) * (z3 - contourValue) < 0) {
			double t = (contourValue - z2) / (z3 - z2);
			meshes::Point intersectionPoint(p2.x + t * (p3.x - p2.x), p2.y + t * (p3.y - p2.y));
			intersections.push_back(intersectionPoint);
		}
		if ((z3 - contourValue) * (z1 - contourValue) < 0) {
			double t = (contourValue - z3) / (z1 - z3);
			meshes::Point intersectionPoint(p3.x + t * (p1.x - p3.x), p3.y + t * (p1.y - p3.y));
			intersections.push_back(intersectionPoint);
		}
		if (intersections.size() == 2) {
			/*if (isPointInsideCurve(intersections[0]) || isPointInsideCurve(intersections[1])) { std::cout << "inside 1" << std::endl; continue; };
			if (isInside(intersections[0]) || isInside(intersections[1])) { std::cout << "inside 2" << std::endl; continue; };*/
			Segment segment(ID, intersections[0], intersections[1]);
			contourSegments.push_back(segment);
			//std::cout << "Segment added" << std::endl;
		}
	}
	//
	//auto t2 = std::chrono::high_resolution_clock::now();
	//std::cout << "Elapsed time (s) for the " << contourSegments.size() << " contour Segments: " << (t2 - t1).count() / 1e6 << "ms"<<std::endl;
	//std::cout << "==================================================================" << std::endl << std::endl;
	return contourSegments;
}

/**
===================================================================================================================================
* @brief Private  function to simlify the isolines (from segments to points).

*
* @param[in] std::vector<Segment>& the vector containing segments representing the isolines
*
* @return a vector of ordered points epresenting the isolines
=====================================================================================================================================*/
std::vector<Point> Isolines::createPointsList(const std::vector<Segment>& segments) {
	//auto t1 = std::chrono::high_resolution_clock::now();
	std::vector<Point> zePoints = {};
	if (segments.size() == 0) {
		return zePoints;
	}
	for (auto s : segments) {
		zePoints.push_back(s.p1);
		zePoints.push_back(s.p2);
	}
	//return zePoints;
	return reorderPoints(zePoints);
}


/**
===================================================================================================================================
* @brief Private  function to reorder points using a nearest neighbor algorithm to sort the points.

*
* @param[in] std::vector<Point>& the vector containing unordered points representing the isolines
*
* @return a vector of ordered points representing the isolines
=====================================================================================================================================*/
std::vector<Point> Isolines::reorderPoints(const std::vector<Point>& points) {
	//auto t1 = std::chrono::high_resolution_clock::now();
	std::vector<Point> reorderedPoints;
	std::vector<Point> remainingPoints = points;

	// Start with the first point in the input vector
	reorderedPoints.push_back(remainingPoints[0]);
	remainingPoints.erase(remainingPoints.begin());

	// Repeat until all points are reordered
	while (!remainingPoints.empty()) {
		// Find the nearest point to the last point in the reordered list
		double minDistance = std::numeric_limits<double>::max();
		size_t nearestIndex = 0;

		for (size_t i = 0; i < remainingPoints.size(); ++i) {
			double d = distance(reorderedPoints.back(), remainingPoints[i]);
			if (d < minDistance) {
				minDistance = d;
				nearestIndex = i;
			}
		}

		// Add the nearest point to the reordered list and remove it from the remaining points
		reorderedPoints.push_back(remainingPoints[nearestIndex]);
		remainingPoints.erase(remainingPoints.begin() + nearestIndex);
	}
	//auto t2 = std::chrono::high_resolution_clock::now();

	//std::cout << "Elapsed time (s) for reordering the points : " << (t2 - t1).count() / 1e6 << " ms."<< std::endl;
	return reorderedPoints;
}

/**
===================================================================================================================================
* @brief Private  function calculating the 2D Euclidean distance.

*
* @param[in] Two "Point"'s
*
* @return a real valued distance
=====================================================================================================================================*/
double Isolines::distance(const Point& p1, const Point& p2) {
	double dx = p1.x - p2.x;
	double dy = p1.y - p2.y;
	return std::sqrt(dx * dx + dy * dy);
}

