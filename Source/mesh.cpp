#include"mesh.h"

using namespace meshes;
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

/**
===================================================================================================================================
* @brief The Mesh class creates a Delaunay triangulation ( a grid) which is use to calculate the redshift isolines in the Isolines class
* @brief Everything happens in the observers frame
*
* @param[in] x_coords_ The cartesian x-component coordinates
* @param[in] y_coords_ The cartesian y-component coordinates
* @param[in] redshift_field_ The redshift factor at that point
* @param[in] bh_mass The mass of the black hole
* @param[in] xisco The cartesian x-component coordinates of the Innermost Stable Circular Orbit
* @param[in] yisco The cartesian y-component coordinates of the Innermost Stable Circular Orbit
*
*
* @return creates a vector of triangles enriched with necessary information for the calculation of the isolines
=====================================================================================================================================*/
Mesh::Mesh(std::vector<double> const& x_coords_, std::vector<double> const& y_coords_, std::vector<double> const& redshift_field_, const std::vector<double>& xisco, const std::vector<double>& yisco) :
	coords(x_coords_),
	x_coords(x_coords_),
	y_coords(y_coords_),
	redshift_field(redshift_field_),
	xISCO(xisco),
	yISCO(yisco),
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
	makeISCO();
	make_driehoeken();
}

/**
===================================================================================================================================
* @brief a copy constructor of the Mesh class creates a Delaunay triangulation ( a grid) which is used to calculate the redshift isolines in the Isolines class
* @brief Everything happens in the observers frame
*
* @param[in] x_coords_ The cartesian x-component coordinates
* @param[in] y_coords_ The cartesian y-component coordinates
* @param[in] redshift_field_ The redshift factor at that point
* @param[in] bh_mass The mass of the black hole
* @param[in] xisco The cartesian x-component coordinates of the Innermost Stable Circular Orbit
* @param[in] yisco The cartesian y-component coordinates of the Innermost Stable Circular Orbit
*
*
* @return creates a vector of triangles enriched with necessary information for the calculation of the isolines
=====================================================================================================================================*/
Mesh::Mesh(const Mesh& other) :
	coords(other.coords),
	triangles(other.triangles),
	driehoeken(other.driehoeken),
	x_coords(other.x_coords),
	y_coords(other.y_coords),
	redshift_field(other.redshift_field),
	xISCO(other.xISCO),
	yISCO(other.yISCO),
	xIscoMax(other.xIscoMax),
	xIscoMin(other.xIscoMin),
	yIscoMax(other.yIscoMax),
	yIscoMin(other.yIscoMin),
	ISCO(other.ISCO),
	IscoEpsilon(other.IscoEpsilon),
	concaveHull(other.concaveHull),
	nConcaveHull(other.nConcaveHull),
	halfedges(other.halfedges),
	hull_prev(other.hull_prev),
	hull_next(other.hull_next),
	hull_tri(other.hull_tri),
	hull_start(other.hull_start)
{
}



std::vector<std::size_t> Mesh::getMesh() {
	return triangles;
}

std::vector<double> Mesh::getCoords() {
	return coords;
}
/**
===================================================================================================================================
* @brief The core function to create a Delaunay mesh.
*
*@param[in] None, uses private class members
* 
*
* @return creates a vector of bare triangles 
=====================================================================================================================================*/
void Mesh::triangulate() {
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

double Mesh::get_hull_area() {
	std::vector<double> hull_area;
	size_t e = hull_start;
	do {
		hull_area.push_back((coords[2 * e] - coords[2 * hull_prev[e]]) * (coords[2 * e + 1] + coords[2 * hull_prev[e] + 1]));
		e = hull_next[e];
	} while (e != hull_start);
	return sum(hull_area);
}

std::size_t Mesh::legalize(std::size_t a) {
	std::size_t i = 0;
	std::size_t ar = 0;
	m_edge_stack.clear();

	// recursion eliminated with a fixed-size stack
	while (true) {
		const size_t b = halfedges[a];

		/* if the pair of triangles doesn't satisfy the Mesh condition
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

inline std::size_t Mesh::hash_key(const double x, const double y) const {
	const double dx = x - m_center_x;
	const double dy = y - m_center_y;
	return fast_mod(
		static_cast<std::size_t>(std::llround(std::floor(pseudo_angle(dx, dy) * static_cast<double>(m_hash_size)))),
		m_hash_size);
}

std::size_t Mesh::add_triangle(
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

void Mesh::link(const std::size_t a, const std::size_t b) {
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

std::vector<size_t> Mesh::get_hull_points() {
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

std::vector<double> Mesh::get_hull_coords() {
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

double Mesh::edge_length(size_t e_a) {
	size_t e_b = next_halfedge(e_a);

	double x_a = coords[2 * triangles[e_a]];
	double y_a = coords[2 * triangles[e_a] + 1];

	double x_b = coords[2 * triangles[e_b]];
	double y_b = coords[2 * triangles[e_b] + 1];

	return sqrt(pow(x_a - x_b, 2) + pow(y_a - y_b, 2));
}

size_t Mesh::get_interior_point(size_t e) {
	return triangles[next_halfedge(next_halfedge(e))];
}
/**
===================================================================================================================================
* @brief Calculates the concave hull of the given points. (as a Delaunay triangulation is convex)
*  @brief As a Delaunay triangulation is convex, "noise" triangles are generated. Knowing the concave hull allows ius to filter out this noise.
*
*@param[in] chi_factor A factor controlling the concavity of the hull
*
*
* @return creates a vector of coordinates used for removing "noise" triangles generated by the Dealaunay triangulation
=====================================================================================================================================*/
std::vector<double> Mesh::concavehull(double chi_factor) {
	//auto t1 = std::chrono::high_resolution_clock::now();
	if (chi_factor < 0 || chi_factor > 1) {
		throw std::invalid_argument("Chi factor must be between 0 and 1 inclusive");
	}

	//mesh::Mesh d(coords);

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
	//auto t2 = std::chrono::high_resolution_clock::now();
	//std::cout << "Elapsed time (s) for the concave hull : " << (t2 - t1).count() / 1e9 << std::endl;
	return conchul;
}

size_t  Mesh::next_halfedge(size_t e) {
	return (e % 3 == 2) ? e - 2 : e + 1;
}

size_t  Mesh::prev_halfedge(size_t e) {
	return (e % 3 == 0) ? e + 2 : e - 1;
}

/**
===================================================================================================================================
* @brief creates a vector of points representing the Inner Stable Circular Orbit @6M
*
*@param[in] None, uses private class members
*
* @return None, feeds the "ISCO" vector class member
=====================================================================================================================================*/
void  Mesh::makeISCO() {
	xIscoMax = *std::max_element(xISCO.begin(), xISCO.end());
	xIscoMin = *std::min_element(xISCO.begin(), xISCO.end());
	yIscoMax = *std::max_element(yISCO.begin(), yISCO.end());
	yIscoMin = *std::min_element(yISCO.begin(), yISCO.end());
	//IscoEpsilon = 0.1;
	for (int i = 0; i < xISCO.size(); i++) {
		Point p(xISCO[i], yISCO[i]);
		ISCO.push_back(p);
	}
}

/**
===================================================================================================================================
* @brief creates a vector of triangles enriched with necessary information for the calculation of the isolines 
*@brief (as opposed to the bare Delaunay triangles conatining only spacial info)
*
*@param[in] None, uses private class members
*
* @return None: creates a class memebr vector of triangles "driehoeken" enriched with necessary information for the calculation of the isolines
=====================================================================================================================================*/
void Mesh::make_driehoeken() {
	//auto t1 = std::chrono::high_resolution_clock::now();
	driehoeken.clear();
	size_t nt = 0;
	for (std::size_t i = 0; i < triangles.size(); i += 3) {
		auto p1 = Point(coords[2 * triangles[i]], coords[2 * triangles[i] + 1]);
		auto p2 = Point(coords[2 * triangles[i + 1]], coords[2 * triangles[i + 1] + 1]);
		auto p3 = Point(coords[2 * triangles[i + 2]], coords[2 * triangles[i + 2] + 1]);
		auto p_middle = Point((p1.x + p2.x + p3.x) / 3.0, (p1.y + p2.y + p3.y) / 3.0);

		if (isOutsideConcaveHull(p_middle)) { continue; }//the (noise)  triangles generated by the convex hull can be situated beyond the largest radius.
		if (isPointInsideCurve(p_middle)) { continue; }//this (noise)  triangle is skipped because it lies beyond the black hole ISCO.
		driehoeken.push_back(driehoek(nt, p1, p2, p3, redshift_field[triangles[i]], redshift_field[triangles[i + 1]], redshift_field[triangles[i + 2]]));
		nt++;
	}
	//createAdjancenyMatrix();
}

/**
===================================================================================================================================
* @brief Checks if a  triangle is inside the ISCO region (those are filtered out).
*
*@param[in] Three points q11,q2,q3 of a triangle and a Point qm the centroid of a tiangle
*
*
* @return a boolean (true if triangle ovelaps the ISCO region)
=====================================================================================================================================*/
bool  Mesh::isPointInsideCurve(const Point& qm) {

	int crossingsm = 0; // Number of times the ray crosses the curve

	// Iterate through each segment of the curve
	for (size_t i = 0; i < ISCO.size(); ++i) {
		const Point& p1 = ISCO[i];
		const Point& p2 = ISCO[(i + 1) % ISCO.size()]; // Wrap around for the last segment

		// Check if the ray from qm intersects the segment (p1, p2)
		if ((p1.y > qm.y) != (p2.y > qm.y) &&
			qm.x < (p2.x - p1.x) * (qm.y - p1.y) / (p2.y - p1.y) + p1.x) {
			crossingsm++;
		}
	}
	// If the number of crossings is odd, the point is inside the curve
	if ( (crossingsm % 2 == 1))return true;
	
	return false;
}




std::vector<std::size_t> Mesh::get_triangles() {
	return triangles;
}
std::vector<double> Mesh::get_tri_coordinates() {
	return coords;
}

bool Mesh::is_left(const Point& a, const Point& b, const Point& c) {
	return ((b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y)) > 0;
}

/**
===================================================================================================================================
* @brief hecks if a triangle is outside the concave hull (helps to remove "noise" triangles generated by the Delaunay triangulation)
*
*@param[in] Point p1
*
*
* @return a boolean (true if point is on or out the the concave hull )
=====================================================================================================================================*/
bool Mesh::isOutsideConcaveHull(const Point& point) {
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

