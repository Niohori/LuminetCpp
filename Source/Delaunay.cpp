#include"delaunay.h"

using namespace delaunay;

bool operator ==(Segment s1, Segment s2)
{
	return((EQ(s1.p1.x, s2.p1.x) && EQ(s1.p1.y, s2.p1.y) && EQ(s1.p2.x, s2.p2.x) && EQ(s1.p2.y, s2.p2.y)) ||
		(EQ(s1.p1.x, s2.p2.x) && EQ(s1.p1.y, s2.p2.y) && EQ(s1.p2.x, s2.p1.x) && EQ(s1.p2.y, s2.p1.y)));
}

bool operator !=(Segment s1, Segment s2)
{
	return(!(EQ(s1.p1.x, s2.p1.x) && EQ(s1.p1.y, s2.p1.y) && EQ(s1.p2.x, s2.p2.x) && EQ(s1.p2.y, s2.p2.y)) &&
		!(EQ(s1.p1.x, s2.p2.x) && EQ(s1.p1.y, s2.p2.y) && EQ(s1.p2.x, s2.p1.x) && EQ(s1.p2.y, s2.p1.y)));
}

Delaunay::Delaunay(std::vector<double> const& x_coords_, std::vector<double> const& y_coords_, std::vector<double> const& redshift_field_) :
	coords(x_coords_),
	x_coords(x_coords_),
	y_coords(y_coords_),
	redshift_field(redshift_field_),
	triangles(),
	halfedges(),
	hull_prev(),
	hull_next(),
	hull_tri(),
	hull_start(),
	m_hash(),
	m_center_x(),
	m_center_y(),
	m_hash_size(),
	m_edge_stack()
{
	coords.clear();
	for (int i = 0; i < x_coords_.size(); ++i) {
		coords.push_back(x_coords_[i]);
		coords.push_back(y_coords_[i]);
	}
	triangulate();
}
void Delaunay::triangulate() {
	std::size_t n = coords.size() >> 1;

	double max_x = std::numeric_limits<double>::min();
	double max_y = std::numeric_limits<double>::min();
	double min_x = std::numeric_limits<double>::max();
	double min_y = std::numeric_limits<double>::max();
	std::vector<std::size_t> ids;
	ids.reserve(n);

	for (std::size_t i = 0; i < n; i++) {
		const double x = coords[2 * i];
		const double y = coords[2 * i + 1];

		if (x < min_x) min_x = x;
		if (y < min_y) min_y = y;
		if (x > max_x) max_x = x;
		if (y > max_y) max_y = y;

		ids.push_back(i);
	}
	const double cx = (min_x + max_x) / 2;
	const double cy = (min_y + max_y) / 2;
	double min_dist = std::numeric_limits<double>::max();

	std::size_t i0 = INVALID_INDEX;
	std::size_t i1 = INVALID_INDEX;
	std::size_t i2 = INVALID_INDEX;

	// pick a seed point close to the centroid
	for (std::size_t i = 0; i < n; i++) {
		const double d = dist(cx, cy, coords[2 * i], coords[2 * i + 1]);
		if (d < min_dist) {
			i0 = i;
			min_dist = d;
		}
	}

	const double i0x = coords[2 * i0];
	const double i0y = coords[2 * i0 + 1];

	min_dist = std::numeric_limits<double>::max();

	// find the point closest to the seed
	for (std::size_t i = 0; i < n; i++) {
		if (i == i0) continue;
		const double d = dist(i0x, i0y, coords[2 * i], coords[2 * i + 1]);
		if (d < min_dist && d > 0.0) {
			i1 = i;
			min_dist = d;
		}
	}

	double i1x = coords[2 * i1];
	double i1y = coords[2 * i1 + 1];

	double min_radius = std::numeric_limits<double>::max();

	// find the third point which forms the smallest circumcircle with the first two
	for (std::size_t i = 0; i < n; i++) {
		if (i == i0 || i == i1) continue;

		const double r = circumradius(
			i0x, i0y, i1x, i1y, coords[2 * i], coords[2 * i + 1]);

		if (r < min_radius) {
			i2 = i;
			min_radius = r;
		}
	}

	if (!(min_radius < std::numeric_limits<double>::max())) {
		throw std::runtime_error("not triangulation");
	}

	double i2x = coords[2 * i2];
	double i2y = coords[2 * i2 + 1];

	if (orient(i0x, i0y, i1x, i1y, i2x, i2y)) {
		std::swap(i1, i2);
		std::swap(i1x, i2x);
		std::swap(i1y, i2y);
	}

	std::tie(m_center_x, m_center_y) = circumcenter(i0x, i0y, i1x, i1y, i2x, i2y);

	// sort the points by distance from the seed triangle circumcenter
	std::sort(ids.begin(), ids.end(), compare{ coords, m_center_x, m_center_y });

	// initialize a hash table for storing edges of the advancing convex hull
	m_hash_size = static_cast<std::size_t>(std::llround(std::ceil(std::sqrt(n))));
	m_hash.resize(m_hash_size);
	std::fill(m_hash.begin(), m_hash.end(), INVALID_INDEX);

	// initialize arrays for tracking the edges of the advancing convex hull
	hull_prev.resize(n);
	hull_next.resize(n);
	hull_tri.resize(n);

	hull_start = i0;

	size_t hull_size = 3;

	hull_next[i0] = hull_prev[i2] = i1;
	hull_next[i1] = hull_prev[i0] = i2;
	hull_next[i2] = hull_prev[i1] = i0;

	hull_tri[i0] = 0;
	hull_tri[i1] = 1;
	hull_tri[i2] = 2;

	m_hash[hash_key(i0x, i0y)] = i0;
	m_hash[hash_key(i1x, i1y)] = i1;
	m_hash[hash_key(i2x, i2y)] = i2;

	std::size_t max_triangles = n < 3 ? 1 : 2 * n - 5;
	triangles.reserve(max_triangles * 3);
	halfedges.reserve(max_triangles * 3);
	add_triangle(i0, i1, i2, INVALID_INDEX, INVALID_INDEX, INVALID_INDEX);
	double xp = std::numeric_limits<double>::quiet_NaN();
	double yp = std::numeric_limits<double>::quiet_NaN();
	for (std::size_t k = 0; k < n; k++) {
		const std::size_t i = ids[k];
		const double x = coords[2 * i];
		const double y = coords[2 * i + 1];

		// skip near-duplicate points
		if (k > 0 && check_pts_equal(x, y, xp, yp)) continue;
		xp = x;
		yp = y;

		// skip seed triangle points
		if (
			check_pts_equal(x, y, i0x, i0y) ||
			check_pts_equal(x, y, i1x, i1y) ||
			check_pts_equal(x, y, i2x, i2y)) continue;

		// find a visible edge on the convex hull using edge hash
		std::size_t start = 0;

		size_t key = hash_key(x, y);
		for (size_t j = 0; j < m_hash_size; j++) {
			start = m_hash[fast_mod(key + j, m_hash_size)];
			if (start != INVALID_INDEX && start != hull_next[start]) break;
		}

		start = hull_prev[start];
		size_t e = start;
		size_t q;

		while (q = hull_next[e], !orient(x, y, coords[2 * e], coords[2 * e + 1], coords[2 * q], coords[2 * q + 1])) { //TODO: does it works in a same way as in JS
			e = q;
			if (e == start) {
				e = INVALID_INDEX;
				break;
			}
		}

		if (e == INVALID_INDEX) continue; // likely a near-duplicate point; skip it

		// add the first triangle from the point
		std::size_t t = add_triangle(
			e,
			i,
			hull_next[e],
			INVALID_INDEX,
			INVALID_INDEX,
			hull_tri[e]);

		hull_tri[i] = legalize(t + 2);
		hull_tri[e] = t;
		hull_size++;

		// walk forward through the hull, adding more triangles and flipping recursively
		std::size_t next = hull_next[e];
		while (
			q = hull_next[next],
			orient(x, y, coords[2 * next], coords[2 * next + 1], coords[2 * q], coords[2 * q + 1])) {
			t = add_triangle(next, i, q, hull_tri[i], INVALID_INDEX, hull_tri[next]);
			hull_tri[i] = legalize(t + 2);
			hull_next[next] = next; // mark as removed
			hull_size--;
			next = q;
		}

		// walk backward from the other side, adding more triangles and flipping
		if (e == start) {
			while (
				q = hull_prev[e],
				orient(x, y, coords[2 * q], coords[2 * q + 1], coords[2 * e], coords[2 * e + 1])) {
				t = add_triangle(q, i, e, INVALID_INDEX, hull_tri[e], hull_tri[q]);
				legalize(t + 2);
				hull_tri[q] = t;
				hull_next[e] = e; // mark as removed
				hull_size--;
				e = q;
			}
		}

		// update the hull indices
		hull_prev[i] = e;
		hull_start = e;
		hull_prev[next] = i;
		hull_next[e] = i;
		hull_next[i] = next;

		m_hash[hash_key(x, y)] = i;
		m_hash[hash_key(coords[2 * e], coords[2 * e + 1])] = e;
	}
	auto a = concavehull(0.1);
	a.clear();

}

double Delaunay::get_hull_area() {
	std::vector<double> hull_area;
	size_t e = hull_start;
	do {
		hull_area.push_back((coords[2 * e] - coords[2 * hull_prev[e]]) * (coords[2 * e + 1] + coords[2 * hull_prev[e] + 1]));
		e = hull_next[e];
	} while (e != hull_start);
	return sum(hull_area);
}

std::size_t Delaunay::legalize(std::size_t a) {
	std::size_t i = 0;
	std::size_t ar = 0;
	m_edge_stack.clear();

	// recursion eliminated with a fixed-size stack
	while (true) {
		const size_t b = halfedges[a];

		/* if the pair of triangles doesn't satisfy the Delaunay condition
		 * (p1 is inside the circumcircle of [p0, pl, pr]), flip them,
		 * then do the same check/flip recursively for the new pair of triangles
		 *
		 *           pl                    pl
		 *          /||\                  /  \
		 *       al/ || \bl            al/    \a
		 *        /  ||  \              /      \
		 *       /  a||b  \    flip    /___ar___\
		 *     p0\   ||   /p1   =>   p0\---bl---/p1
		 *        \  ||  /              \      /
		 *       ar\ || /br             b\    /br
		 *          \||/                  \  /
		 *           pr                    pr
		 */
		const size_t a0 = 3 * (a / 3);
		ar = a0 + (a + 2) % 3;

		if (b == INVALID_INDEX) {
			if (i > 0) {
				i--;
				a = m_edge_stack[i];
				continue;
			}
			else {
				//i = INVALID_INDEX;
				break;
			}
		}

		const size_t b0 = 3 * (b / 3);
		const size_t al = a0 + (a + 1) % 3;
		const size_t bl = b0 + (b + 2) % 3;

		const std::size_t p0 = triangles[ar];
		const std::size_t pr = triangles[a];
		const std::size_t pl = triangles[al];
		const std::size_t p1 = triangles[bl];

		const bool illegal = in_circle(
			coords[2 * p0],
			coords[2 * p0 + 1],
			coords[2 * pr],
			coords[2 * pr + 1],
			coords[2 * pl],
			coords[2 * pl + 1],
			coords[2 * p1],
			coords[2 * p1 + 1]);

		if (illegal) {
			triangles[a] = p1;
			triangles[b] = p0;

			auto hbl = halfedges[bl];

			// edge swapped on the other side of the hull (rare); fix the halfedge reference
			if (hbl == INVALID_INDEX) {
				std::size_t e = hull_start;
				do {
					if (hull_tri[e] == bl) {
						hull_tri[e] = a;
						break;
					}
					e = hull_next[e];
				} while (e != hull_start);
			}
			link(a, hbl);
			link(b, halfedges[ar]);
			link(ar, bl);
			std::size_t br = b0 + (b + 1) % 3;

			if (i < m_edge_stack.size()) {
				m_edge_stack[i] = br;
			}
			else {
				m_edge_stack.push_back(br);
			}
			i++;

		}
		else {
			if (i > 0) {
				i--;
				a = m_edge_stack[i];
				continue;
			}
			else {
				break;
			}
		}
	}
	return ar;
}

inline std::size_t Delaunay::hash_key(const double x, const double y) const {
	const double dx = x - m_center_x;
	const double dy = y - m_center_y;
	return fast_mod(
		static_cast<std::size_t>(std::llround(std::floor(pseudo_angle(dx, dy) * static_cast<double>(m_hash_size)))),
		m_hash_size);
}

std::size_t Delaunay::add_triangle(
	std::size_t i0,
	std::size_t i1,
	std::size_t i2,
	std::size_t a,
	std::size_t b,
	std::size_t c) {
	std::size_t t = triangles.size();
	triangles.push_back(i0);
	triangles.push_back(i1);
	triangles.push_back(i2);
	link(t, a);
	link(t + 1, b);
	link(t + 2, c);
	return t;
}

void Delaunay::link(const std::size_t a, const std::size_t b) {
	std::size_t s = halfedges.size();
	if (a == s) {
		halfedges.push_back(b);
	}
	else if (a < s) {
		halfedges[a] = b;
	}
	else {
		throw std::runtime_error("Cannot link edge");
	}
	if (b != INVALID_INDEX) {
		std::size_t s2 = halfedges.size();
		if (b == s2) {
			halfedges.push_back(a);
		}
		else if (b < s2) {
			halfedges[b] = a;
		}
		else {
			throw std::runtime_error("Cannot link edge");
		}
	}
}

std::vector<size_t> Delaunay::get_hull_points() {
	std::vector<size_t> hull_pts;

	size_t point = hull_start;
	do {
		hull_pts.push_back(point);
		point = hull_next[point];
	} while (point != hull_start);

	// Wrap back around
	hull_pts.push_back(hull_start);

	return hull_pts;
}

std::vector<double> Delaunay::get_hull_coords() {
	const std::vector<size_t> hull_pts = get_hull_points();

	std::vector<double> hull_coords;
	hull_coords.reserve(2 * hull_pts.size());

	for (size_t point : hull_pts) {
		double x = coords[2 * point];
		double y = coords[2 * point + 1];
		hull_coords.push_back(x);
		hull_coords.push_back(y);
	}

	return hull_coords;
}

double Delaunay::edge_length(size_t e_a) {
	size_t e_b = next_halfedge(e_a);

	double x_a = coords[2 * triangles[e_a]];
	double y_a = coords[2 * triangles[e_a] + 1];

	double x_b = coords[2 * triangles[e_b]];
	double y_b = coords[2 * triangles[e_b] + 1];

	return sqrt(pow(x_a - x_b, 2) + pow(y_a - y_b, 2));
}

size_t Delaunay::get_interior_point(size_t e) {
	return triangles[next_halfedge(next_halfedge(e))];
}


std::vector<double> Delaunay::concavehull(double chi_factor) {
	if (chi_factor < 0 || chi_factor > 1) {
		throw std::invalid_argument("Chi factor must be between 0 and 1 inclusive");
	}

	//delaunay::Delaunay d(coords);

	// Determine initial points on outside hull
	std::vector<size_t> bpoints = get_hull_points();
	std::set<size_t> bset(bpoints.begin(), bpoints.end());

	// Make max heap of boundary edges with lengths
	typedef std::pair<size_t, double> hpair;

	auto cmp = [](hpair left, hpair right) {
		return left.second < right.second;
		};

	std::vector<hpair> bheap(bpoints.size());

	double max_len = std::numeric_limits<double>::min();
	double min_len = std::numeric_limits<double>::max();

	for (auto point : bpoints) {
		size_t e = hull_tri[point];
		double len = edge_length(e);

		bheap.push_back({ e, len });
		std::push_heap(bheap.begin(), bheap.end(), cmp);

		min_len = std::min(len, min_len);
		max_len = std::max(len, max_len);
	}

	// Determine length parameter
	double length_param = chi_factor * max_len + (1 - chi_factor) * min_len;

	// Iteratively add points to boundary by iterating over the triangles on the hull
	while (!bheap.empty()) {

		// Get edge with the largest length
		std::pop_heap(bheap.begin(), bheap.end(), cmp);
		const auto [e, len] = bheap.back();
		bheap.pop_back();

		// Length of edge too small for our chi factor
		if (len <= length_param) {
			break;
		}

		// Find interior point given edge e (a -> b)
		//       e
		//  b <----- a
		//     \   /
		//  e_b \ / e_a
		//       c
		size_t c = get_interior_point(e);

		// Point already belongs to boundary
		if (bset.count(c)) {
			continue;
		}

		// Get two edges connected to interior point
		//  c -> b
		size_t e_b = halfedges[next_halfedge(e)];
		//  a -> c
		size_t e_a = halfedges[next_halfedge(next_halfedge(e))];

		// Add edges to heap
		double len_a = edge_length(e_a);
		double len_b = edge_length(e_b);

		bheap.push_back({ e_a, len_a });
		std::push_heap(bheap.begin(), bheap.end(), cmp);
		bheap.push_back({ e_b, len_b });
		std::push_heap(bheap.begin(), bheap.end(), cmp);

		// Update outer hull and connect new edges
		size_t a = triangles[e];
		size_t b = triangles[next_halfedge(e)];

		hull_next[c] = b;
		hull_prev[c] = a;
		hull_next[a] = hull_prev[b] = c;

		bset.insert(c);
	}
	auto conchul = get_hull_coords();
	std::vector<Point> hull;

	for (std::size_t i = 0; i < conchul.size(); i += 2)
	{
		const double x = conchul[i];
		const double y = conchul[i + 1];
		hull.push_back(Point(x, y));
	}
	concaveHull = hull;
	nConcaveHull = concaveHull.size();
	hull.clear();
	return conchul;
}

size_t  Delaunay::next_halfedge(size_t e) {
	return (e % 3 == 2) ? e - 2 : e + 1;
}

size_t  Delaunay::prev_halfedge(size_t e) {
	return (e % 3 == 0) ? e + 2 : e - 1;
}

std::vector<Segment> Delaunay::get_iso_lines(const double& requested_level, const std::vector<double>& xisco, const std::vector<double>& yisco) {
	double IscoMax = *std::max_element(xisco.begin(), xisco.end());
	double IscoMin = *std::min_element(xisco.begin(), xisco.end());
	double EPSILON = 0;
	driehoeken.clear();
	//std::pair<std::vector<Point>, std::vector<Point>> isolines;
	std::vector<Segment>all_isolines;

	size_t nt = 0;
	double average_area = 0.0;
	double max_area = -10000000.0;
	double min_area = 10000000.0;

	//create ISCO container;
	std::vector<Point> ISCO;
	//std::cout << "isco sizes : " << xisco.size() << ",   " << yisco.size() << std::endl;
	for (int i = 0; i < xisco.size(); i++) {
		Point p(xisco[i], yisco[i]);
		ISCO.push_back(p);
	}
	for (std::size_t i = 0; i < triangles.size(); i += 3) {
		auto p1 = Point(coords[2 * triangles[i]], coords[2 * triangles[i] + 1]);
		auto p2 = Point(coords[2 * triangles[i + 1]], coords[2 * triangles[i + 1] + 1]);
		auto p3 = Point(coords[2 * triangles[i + 2]], coords[2 * triangles[i + 2] + 1]);

		auto p_middle = Point((p1.x + p2.x + p3.x) / 3.0, (p1.y + p2.y + p3.y) / 3.0);

		if (isOutsideConcaveHull(p_middle)) { continue; }//the (noise)  triangles generated by the convex hull are because they lie beyond the largest radius.
		if (isOnPolygon(EPSILON, p1, p2, p3, xisco, yisco, IscoMax, IscoMin)) { continue; }//this (noise)  triangle is skipped because it lies beyond the black hole ISCO.
		driehoeken.push_back(driehoek(nt, p1, p2, p3, redshift_field[triangles[i]], redshift_field[triangles[i + 1]], redshift_field[triangles[i + 2]]));
		nt++;
	}
	all_isolines = generateContourLines(driehoeken, requested_level);

	return all_isolines;
}

// Function to generate contour lines
std::vector<Segment> Delaunay::generateContourLines(const std::vector<driehoek>& driehoeken, double contourValue) {
	//std::vector<Point> Delaunay::generateContourLines(const std::vector<driehoek>&driehoeken, double contourValue) {
	std::vector<Point> contourLines;
	std::vector<Segment> contourSegments;

	for (const driehoek& driehk : driehoeken) {
		double z1 = driehk.Rs1;
		double z2 = driehk.Rs2;
		double z3 = driehk.Rs3;
		Point p1 = driehk.p1;
		Point p2 = driehk.p2;
		Point p3 = driehk.p3;
		auto ID = driehk.ID;
		std::vector<Point> intersections;
		// Check each edge of the triangle
		if ((z1 - contourValue) * (z2 - contourValue) < 0) {
			double t = (contourValue - z1) / (z2 - z1);
			Point intersectionPoint(p1.x + t * (p2.x - p1.x), p1.y + t * (p2.y - p1.y));
			//contourLines.push_back(intersectionPoint);
			intersections.push_back(intersectionPoint);
		}
		if ((z2 - contourValue) * (z3 - contourValue) < 0) {
			double t = (contourValue - z2) / (z3 - z2);
			Point intersectionPoint(p2.x + t * (p3.x - p2.x), p2.y + t * (p3.y - p2.y));
			//contourLines.push_back(intersectionPoint);
			intersections.push_back(intersectionPoint);
		}
		if ((z3 - contourValue) * (z1 - contourValue) < 0) {
			double t = (contourValue - z3) / (z1 - z3);
			Point intersectionPoint(p3.x + t * (p1.x - p3.x), p3.y + t * (p1.y - p3.y));
			//contourLines.push_back(intersectionPoint);
			intersections.push_back(intersectionPoint);
		}
		if (intersections.size() == 2) {

			Segment segment(ID, intersections[0], intersections[1]);
			//if (isNoise(segment) ){
			//	//continue;
			//	//std::cout << "identical x's " << std::endl;
			//}
			contourSegments.push_back(segment);
			//std::cout << "Segment added" << std::endl;
		}
	}

	//return contourLines;
	return contourSegments;
}
bool Delaunay::isNoise(Segment seg) {
	double diff = std::abs(seg.p1.x - seg.p2.x);
	if (diff > 0.05)return false;
	return true;

}
/*=============================================================================================================

Removes the triangle inside the black Hole silhouet.

==============================================================================================================*/
bool Delaunay::isOnPolygon(const double& epsilon, const Point& p1, const Point& p2, const Point& p3, const std::vector<double>& X, const  std::vector<double>& Y, const double& xIscoMax, const double& xIscoMin) {
	int nCross = 0;
	double d1 = 0;
	double d2 = 0;
	double d3 = 0;
	for (int i = 0; i < X.size(); i++) {
		d1 = dist(p1.x, p1.y, X[i], Y[i]);
		d2 = dist(p1.x, p1.y, X[i], Y[i]);
		d3 = dist(p1.x, p1.y, X[i], Y[i]);
		if (d1 <= epsilon)nCross++;
		if (d2 <= epsilon)nCross++;
		if (d3 <= epsilon)nCross++;
	}
	if (nCross >= 2) {//edges of the triangle cross the ISCO line
		return true;
	}
	if ((xIscoMin <= p1.x && p1.x <= xIscoMax) || (xIscoMin <= p2.x && p2.x <= xIscoMax) || (xIscoMin <= p3.x && p3.x <= xIscoMax)) {
		if ((p1.y * p2.y < 0) || (p1.y * p3.y < 0) || (p3.y * p2.y < 0)) {
			return true;
		}

	}

	return false;
}

std::vector<std::size_t> Delaunay::get_triangles() {

	return triangles;
}
std::vector<double> Delaunay::get_tri_coordinates() {

	return coords;
}

bool Delaunay::is_left(const Point& a, const Point& b, const Point& c) {
	return ((b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y)) > 0;
}

bool Delaunay::isOutsideConcaveHull(const Point& point) {
	//bool Delaunay::isOutsideConcaveHull(const Point & point, const std::vector<Point>&concave_hull) 
	bool inside = true;

	for (int i = 0, j = nConcaveHull - 1; i < nConcaveHull; j = i++) {
		if (((concaveHull[i].y <= point.y && point.y < concaveHull[j].y) ||
			(concaveHull[j].y <= point.y && point.y < concaveHull[i].y)) &&
			(point.x < (concaveHull[j].x - concaveHull[i].x) * (point.y - concaveHull[i].y) /
				(concaveHull[j].y - concaveHull[i].y) +
				concaveHull[i].x)) {
			inside = !inside;
		}
	}

	return inside;
}