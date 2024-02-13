#include "isolines.h"

Isolines::Isolines(const std::shared_ptr<meshes::Mesh>& meshptr, const double& isovalue) :
	meshPointer(meshptr),
	contourValue(isovalue)
{
	;
}

std::vector<meshes::Point> Isolines::get_iso_lines() {
	//auto t1 = std::chrono::high_resolution_clock::now();
	std::vector<meshes::Segment>all_isolines;

	size_t nt = 0;
	double average_area = 0.0;
	double max_area = -10000000.0;
	double min_area = 10000000.0;

	//auto t2 = std::chrono::high_resolution_clock::now();
	//std::cout << "Elapsed time (s) in get_iso_lines to populate the -- driehoeken --: " << (t2 - t1).count() / 1e9 << std::endl;
	all_isolines = generateContourLines();
	//mergeSegments(all_isolines);
	std::vector<Point> points = createPointsList(all_isolines);
	//std::cout << "Number of segments : " << all_isolines.size() << " --> expected number of points : " << 2 * all_isolines.size() << " and got " << points.size() << " points." << std::endl;

	return points;
}

// Function to generate contour lines
std::vector<meshes::Segment> Isolines::generateContourLines() {
	//auto t1 = std::chrono::high_resolution_clock::now();
	//std::vector<Point> Isolines::generateContourLines(const std::vector<driehoek>&driehoeken, double contourValue) {
	std::vector<meshes::Segment> contourSegments;

	for (const driehoek& driehk : meshPointer->driehoeken) {//search for neighbors for each driehoek
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
			intersections.push_back(intersectionPoint);
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
			Segment segment(ID, intersections[0], intersections[1]);
			contourSegments.push_back(segment);
			//std::cout << "Segment added" << std::endl;
		}
	}
	
	//auto t2 = std::chrono::high_resolution_clock::now();
	//std::cout << "Elapsed time (s) for the " << contourSegments.size() << " contour Segments: " << (t2 - t1).count() / 1e9 << std::endl;
	return contourSegments;
}
bool Isolines::isNoise(Segment seg) {
	double diff = std::abs(seg.p1.x - seg.p2.x);
	if (diff > 0.05)return false;
	return true;
}

void  Isolines::mergeSegments(std::vector<Segment>& segments) {
	//auto t1 = std::chrono::high_resolution_clock::now();
	std::vector<Segment>::iterator it;
	std::vector<Segment>::iterator jt;
	if (segments.size() < 2) return;
	it = segments.begin();
	/*
	using two iterators we walk through the entire vector testing to
	see if some combination of the start and end points match. If we
	find matching points the two vectorsrs are merged. since when we go
	thru the vector once we are garanteed that all vectors that can
	connect to that oone have been merged we only have to merge the
	vector less the processed nodes at the begining. every merge does
	force jt back to the beginning of the search tho since a merge will
	change either the start or the end of the vector
	*/
	while (it != segments.end())
	{
		jt = it + 1;
		while (jt != segments.end())
		{
			/*
			if the end of *it matches the start ot *jt we can copy
			*jt to the end of *it and remove jt, the erase funtion
			does us a favour and increments the iterator to the
			next element so we continue to test the next element
			*/
			if ((*it).p2 == (*jt).p1 || (*it).p1 == (*jt).p2)
			{
				segments.erase(jt);
				jt = it + 1;
			}
			else
			{
				jt++;
			}
		}
		it++;
	}
	//auto t2 = std::chrono::high_resolution_clock::now();
	//std::cout << "Elapsed time (s) for the merge of segments: " << (t2 - t1).count() / 1e9 << std::endl;
	return;
}

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
	return reorderPoints(zePoints);
	Point startingPoint(0.0, 0.0);
	int startIndex = 0;
	int closestSegment = 0;
	double minDist = 100000000000000000000000.0;
	int id = 0;
	for (auto s : segments) {
		for (auto x : meshPointer->ISCO) {
			double dist1 = s.p1.dist(x);
			double dist2 = s.p2.dist(x);
			if (dist1 < minDist) {
				minDist = dist1;
				closestSegment = startIndex;
				id = 1;
			}
			if (dist2 < minDist) {
				minDist = dist2;
				closestSegment = startIndex;
				id = 2;
			}
		}
		startIndex++;
	}
	//std::cout << "closestSegment index = " << closestSegment << "to compare with size of segment of " << segments.size() << std::endl;
	Segment candidate = segments[closestSegment];
	if (id == 1) {
		zePoints.push_back(candidate.p1);
		zePoints.push_back(candidate.p2);
		//std::cout << "Found a point closest to the ISO with coordinates (" << candidate.p1.x << ", " << candidate.p1.y << ")" << std::endl;
	}
	if (id == 2) {
		zePoints.push_back(candidate.p2);
		zePoints.push_back(candidate.p1);
		//std::cout << "Found a point closest to the ISO with coordinates (" << candidate.p2.x << ", " << candidate.p2.y << ")" << std::endl;
	}
	int count = 0;
	for (auto s : segments) {
		if (count == closestSegment)continue;
		zePoints.push_back(s.p1);
		zePoints.push_back(s.p2);
		count++;
	}
	//auto t2 = std::chrono::high_resolution_clock::now();

	//std::cout << "Elapsed time (s) for createPointList(std::vector<Segment>& segments) " << (t2 - t1).count() / 1e9 << std::endl;
	return zePoints;
	//return reorderPoints(zePoints);
}
std::vector<Point> Isolines::reorderPoints(const std::vector<Point>& points) {
	//using a nearest neighbor algorithm to sort the points
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

	//std::cout << "Elapsed time (s) for reordering the points" << (t2 - t1).count() / 1e9 << std::endl;
	return reorderedPoints;
}
double Isolines::distance(const Point& p1, const Point& p2) {
	double dx = p1.x - p2.x;
	double dy = p1.y - p2.y;
	return std::sqrt(dx * dx + dy * dy);
}

//std::pair<std::vector<double>, std::vector<double> > Isolines::smoothCurve(const std::vector<Point>& curve) {
//	std::vector<double> fitted(curve.size()), resids(curve.size());
//	std::vector x = std::vector<double>(curve.size(), 0.0);
//	std::vector y = std::vector<double>(curve.size(), 0.0);
//	for (int i = 0; i < curve.size(); i++) {
//		x[i] = curve[i].x;
//		y[i] = curve[i].y;
//	}
//	WeightedLowess::WeightedLowess smoother;
//	smoother.run(curve.size(), x.data(), y.data(), NULL, fitted.data(), resids.data());
//	//for (int i = 0; i < curve.size(); i++) {
//	//	std::cout << "(" << fitted[i] << ", " << resids[i] << std::endl;
//	//}
//
//	return std::make_pair(x, fitted);
//}

std::vector<std::vector<double> > Isolines::makeAdjacencyMatrix(const std::vector< Segment >& segments) {
	int n = segments.size();
	std::vector<std::vector<double> > adj(n, std::vector<double>(n, 0.0));// Initialize the adjacency matrix with zeros
	auto t1 = std::chrono::high_resolution_clock::now();
	//// Populate the adjacency matrix
	parallel_for(0, n, [&](int i) {
		adj[i][i] = 1000.0;
		parallel_for(i + 1, n, [&](int j) {
			Point p1m = Point(0.5 * (segments[i].p1.x + segments[i].p2.x), 0.5 * (segments[i].p1.y + segments[i].p2.y));
			Point p2m = Point(0.5 * (segments[j].p1.x + segments[j].p2.x), 0.5 * (segments[j].p1.y + segments[j].p2.y));
			double d = distance(p1m, p2m);
			adj[i][j] = d;
			adj[j][i] = d;
			});
		}

	);
	// Ensure graph is connected
	for (int i = 0; i < n; ++i) {
		adj[i][i] = 1000000000000000000000000000.0; // Assuming a large value for non-connected vertices
	}
	auto t2 = std::chrono::high_resolution_clock::now();

	std::cout << "Elapsed time (s) for make distance adjacency matrix " << (t2 - t1).count() / 1e9 << " seconds " << std::endl;
	// Output the adjacency matrix
	//for (int i = 0; i < n; ++i) {
	//	for (int j = 0; j < n; ++j) {
	//		std::cout << adj[i][j] << " ";
	//	}
	//	std::cout << std::endl;
	//	if (i > 1)break;
	//}
	// Find shortest path starting from the first point
	std::vector<int> shortestPathIndices = dijkstra(adj, 0);

	//// Output the shortest path indices
	/*std::cout << "Shortest Path Indices:";
	for (int i : shortestPathIndices) {
		std::cout << " " << i << std::endl;;
	}
	std::cout << std::endl;*/
	return adj;
}

std::vector<int>  Isolines::dijkstra(const std::vector<std::vector<double>>& adjMatrix, int source) {
	auto t1 = std::chrono::high_resolution_clock::now();
	std::vector<int>  shortestPath;
	int n = adjMatrix.size();
	std::vector<double> dist(n, INF); // Initialize distances to infinity
	std::vector<int> prev(n, -1); // Initialize previous vertices

	dist[source] = 0; // Distance from source to itself is 0

	// Priority queue to store vertices with their distances
	std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, std::greater<std::pair<double, int>>> pq;
	pq.push({ 0, source });

	while (!pq.empty()) {
		int u = pq.top().second;
		pq.pop();

		for (int v = 0; v < n; ++v) {
			if (adjMatrix[u][v] > 0) { // If there's an edge from u to v
				double newDist = dist[u] + adjMatrix[u][v];
				if (newDist < dist[v]) {
					dist[v] = newDist;
					prev[v] = u;
					pq.push({ dist[v], v });
				}
			}
		}
	}

	// Print the shortest distances and paths
	for (int i = 0; i < n; ++i) {
		//std::cout << "Distance from vertex " << source << " to " << i << ": " << dist[i] << std::endl;
		// To get the path, backtrack from destination to source
		int curr = i;
		while (curr != -1) {
			if (curr != source) shortestPath.push_back(curr);
			//std::cout << curr << " ";
			curr = prev[curr];
		}
		//std::cout << std::endl;
	}
	auto t2 = std::chrono::high_resolution_clock::now();

	std::cout << "Elapsed time (s) for Dijkstra algo " << (t2 - t1).count() / 1e9 << " seconds " << std::endl;
	return shortestPath;
}