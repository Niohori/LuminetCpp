#pragma once
#ifndef MESH_H
#define MESH_H
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


using namespace dlib;
#define DIFFERENCE 0.0005
#define EQ(_x_,_y_) (((_x_-_y_<DIFFERENCE)&&(_y_-_x_<DIFFERENCE))?1:0)
// Define the maximum value for infinity
const double INF = std::numeric_limits<double>::max();

namespace meshes {
	// Point structure representing a 2D Point
	struct Point {
		double x, y;
		double epsilon = 0.1; // Threshold distance for equality check
		Point(double x_, double y_) : x(x_), y(y_) {}
		// Custom comparison operators for set
		bool operator<(const Point& other) const {
			return x < other.x || (x == other.x && y < other.y);
		}

		bool operator==(const Point& other) const {
			return std::abs(x - other.x) < epsilon && std::abs(y - other.y) < epsilon;
		}
		bool operator!=(const Point& other) const {
			return x != other.x || y != other.y;
		}
		double dist(const Point& other) const {
			return std::sqrt((x - other.x) * (x - other.x) + (y - other.y) * (y - other.y));
		}
	};

	struct driehoek {
		//driehoek() {};
		size_t ID;//identification -> not necessary for the moment being
		size_t ID1;//identification of neighbours Mesh triangle (driehoeken)
		size_t ID2;//identification of neighbours Mesh triangle (driehoeken)
		size_t ID3;//identification of neighbours Mesh triangle (driehoeken)

		size_t nNeighbors;//how many has he got (could be two)
		Point p1;//1st triange Point
		Point p2;//2nd triange Point
		Point p3;//3rd triange Point
		Point Centroid;//centroid triange Point

		double  Rs1;//redshift @ triange Point
		double  Rs2;//redshift @ triange Point
		double  Rs3;//redshift @ triange Point

		driehoek(size_t ID_, Point p1_, Point p2_, Point p3_, double rs1_, double rs2_, double rs3_) :
			ID(ID_),
			p1(p1_),
			p2(p2_),
			p3(p3_),
			Centroid(p1_),
			Rs1(rs1_),
			Rs2(rs2_),
			Rs3(rs3_)
		{
			Centroid = Point((p1.x + p2.x + p3.x / 3), (p1.y + p2.y + p3.y / 3));
			ID1 = std::numeric_limits<size_t>::quiet_NaN();
			ID2 = std::numeric_limits<size_t>::quiet_NaN();
			ID3 = std::numeric_limits<size_t>::quiet_NaN();
		};
		// Function to check if two triangles share an edge
		bool share_edge(const driehoek& other) const {
			std::set<Point> shared_vertices{ p1, p2, p3 };
			return static_cast<int>(shared_vertices.count(other.p1) + shared_vertices.count(other.p2) + shared_vertices.count(other.p3) == 2);
		}
		// Function to check which edge two triangles share
		std::set<Point> shared_edge(const driehoek& other) const {
			std::set<Point> shared_vertices{ p1, p2, p3 };
			std::set<Point> common_vertices;

			for (const Point& vertex : { other.p1, other.p2, other.p3 }) {
				if (shared_vertices.count(vertex) > 0) {
					common_vertices.insert(vertex);
				}
			}

			return common_vertices;
		}
	};

	struct Segment {
		Point p1;
		Point p2;
		int ID;;//identification of Mesh triangle (driehoeken)
		int ID1 = -1;// = std::numeric_limits<size_t>::quiet_NaN();;//identification of neighbours Mesh triangle (driehoeken)
		int ID2 = -1;// = std::numeric_limits<size_t>::quiet_NaN();;//identification of neighbours Mesh triangle (driehoeken)
		std::map<size_t, bool> IDs;
		size_t nNeighbors;// = 0;//how many has he got (could be two)
		double epsilon = 0.1001; // Threshold distance for equality check

		//bool check_equal(Point, Point);
		//Segment(); ;
		double disk = 0.1;
		Segment(size_t ID_, Point p1_, Point p2_) :
			ID(ID_), p1(p1_), p2(p2_)
		{
			ID1 = -1;
			ID2 = -1;
			nNeighbors = 0;
		}
		// Define a function to check if two points are approximately equal
		bool pointsApproximatelyEqual(const Point& a, const Point& b) const {
			double dx = a.x - b.x;
			double dy = a.y - b.y;
			return std::sqrt(dx * dx + dy * dy) < epsilon;
		}
		bool operator==(const Segment& other) const {
			return p1.x == other.p1.x && p1.y == other.p1.y && p2.x == other.p2.x && p2.y == other.p2.y;
		}
		bool operator!=(const Segment& other) const {
			return p1.x != other.p1.x || p1.y != other.p1.y || p2.x != other.p2.x || p2.y != other.p2.y;
		}
		// Function to check if two segments share any common points
		bool shareCommonPoints(const Segment& other) const {
			return (pointsApproximatelyEqual(p1, other.p1) || pointsApproximatelyEqual(p1, other.p2) || pointsApproximatelyEqual(p2, other.p1) || pointsApproximatelyEqual(p2, other.p2));
		}
		// Function to check if two segments share any common points
		std::vector<Point> sharedPoints(const Segment& other) const {
			std::vector<Point> shared;
			if (pointsApproximatelyEqual(p1, other.p1) || pointsApproximatelyEqual(p1, other.p2)) shared.push_back(p1);
			if (pointsApproximatelyEqual(p2, other.p1) || pointsApproximatelyEqual(p2, other.p2)) shared.push_back(p2);
			return shared;
		}
	};

	//@see https://stackoverflow.com/questions/33333363/built-in-mod-vs-custom-mod-function-improve-the-performance-of-modulus-op/33333636#33333636
	inline size_t fast_mod(const size_t i, const size_t c) {
		return i >= c ? i % c : i;
	}

	// Kahan and Babuska summation, Neumaier variant; accumulates less FP error
	inline double sum(const std::vector<double>& x) {
		double sum = x[0];
		double err = 0.0;

		for (size_t i = 1; i < x.size(); i++) {
			const double k = x[i];
			const double m = sum + k;
			err += std::fabs(sum) >= std::fabs(k) ? sum - m + k : k - m + sum;
			sum = m;
		}
		return sum + err;
	}

	inline double dist(
		const double ax,
		const double ay,
		const double bx,
		const double by) {
		const double dx = ax - bx;
		const double dy = ay - by;
		return dx * dx + dy * dy;
	}

	inline double circumradius(
		const double ax,
		const double ay,
		const double bx,
		const double by,
		const double cx,
		const double cy) {
		const double dx = bx - ax;
		const double dy = by - ay;
		const double ex = cx - ax;
		const double ey = cy - ay;

		const double bl = dx * dx + dy * dy;
		const double cl = ex * ex + ey * ey;
		const double d = dx * ey - dy * ex;

		const double x = (ey * bl - dy * cl) * 0.5 / d;
		const double y = (dx * cl - ex * bl) * 0.5 / d;

		if ((bl > 0.0 || bl < 0.0) && (cl > 0.0 || cl < 0.0) && (d > 0.0 || d < 0.0)) {
			return x * x + y * y;
		}
		else {
			return std::numeric_limits<double>::max();
		}
	}

	inline bool orient(
		const double px,
		const double py,
		const double qx,
		const double qy,
		const double rx,
		const double ry) {
		return (qy - py) * (rx - qx) - (qx - px) * (ry - qy) < 0.0;
	}

	inline std::pair<double, double> circumcenter(
		const double ax,
		const double ay,
		const double bx,
		const double by,
		const double cx,
		const double cy) {
		const double dx = bx - ax;
		const double dy = by - ay;
		const double ex = cx - ax;
		const double ey = cy - ay;

		const double bl = dx * dx + dy * dy;
		const double cl = ex * ex + ey * ey;
		const double d = dx * ey - dy * ex;

		const double x = ax + (ey * bl - dy * cl) * 0.5 / d;
		const double y = ay + (dx * cl - ex * bl) * 0.5 / d;

		return std::make_pair(x, y);
	}

	struct compare {
		std::vector<double> const& coords;
		double cx;
		double cy;

		bool operator()(std::size_t i, std::size_t j) {
			const double d1 = dist(coords[2 * i], coords[2 * i + 1], cx, cy);
			const double d2 = dist(coords[2 * j], coords[2 * j + 1], cx, cy);
			const double diff1 = d1 - d2;
			const double diff2 = coords[2 * i] - coords[2 * j];
			const double diff3 = coords[2 * i + 1] - coords[2 * j + 1];

			if (diff1 > 0.0 || diff1 < 0.0) {
				return diff1 < 0;
			}
			else if (diff2 > 0.0 || diff2 < 0.0) {
				return diff2 < 0;
			}
			else {
				return diff3 < 0;
			}
		}
	};

	inline bool in_circle(
		const double ax,
		const double ay,
		const double bx,
		const double by,
		const double cx,
		const double cy,
		const double px,
		const double py) {
		const double dx = ax - px;
		const double dy = ay - py;
		const double ex = bx - px;
		const double ey = by - py;
		const double fx = cx - px;
		const double fy = cy - py;

		const double ap = dx * dx + dy * dy;
		const double bp = ex * ex + ey * ey;
		const double cp = fx * fx + fy * fy;

		return (dx * (ey * cp - bp * fy) -
			dy * (ex * cp - bp * fx) +
			ap * (ex * fy - ey * fx)) < 0.0;
	}

	constexpr double EPSILON = std::numeric_limits<double>::epsilon();
	constexpr std::size_t INVALID_INDEX = std::numeric_limits<std::size_t>::max();

	inline bool check_pts_equal(double x1, double y1, double x2, double y2) {
		return std::fabs(x1 - x2) <= EPSILON &&
			std::fabs(y1 - y2) <= EPSILON;
	}

	// monotonically increases with real angle, but doesn't need expensive trigonometry
	inline double pseudo_angle(const double dx, const double dy) {
		const double p = dx / (std::abs(dx) + std::abs(dy));
		return (dy > 0.0 ? 3.0 - p : 1.0 + p) / 4.0; // [0..1)
	}

	struct MeshPoint {
		std::size_t i;
		double x;
		double y;
		std::size_t t;
		std::size_t prev;
		std::size_t next;
		bool removed;
	};

	class Mesh {
	public:
		//Mesh(void);
		//Mesh(std::vector<double> const&, std::vector<double> const&, std::vector<double> const&);
		Mesh(std::vector<double> const&, std::vector<double> const&, std::vector<double> const&, const std::vector<double>&, const std::vector<double>& yisco);
		//Mesh(const std::vector<double> , const std::vector<std::size_t> );
		Mesh(const Mesh& other);
		// Assignment operator
		Mesh& operator=(const Mesh& other) {
			if (this != &other) {
				coords = other.coords;
				triangles = other.triangles;
				driehoeken = other.driehoeken;
				x_coords = other.x_coords;
				y_coords = other.y_coords;
				redshift_field = other.redshift_field;
				xISCO = other.xISCO;
				yISCO = other.yISCO;
				xIscoMax = other.xIscoMax;
				xIscoMin = other.xIscoMin;
				yIscoMax = other.yIscoMax;
				yIscoMin = other.yIscoMin;
				ISCO = other.ISCO;
				IscoEpsilon = other.IscoEpsilon;
				concaveHull = other.concaveHull;
				nConcaveHull = other.nConcaveHull;
			}
			return *this;
		}

		std::vector<std::size_t> getMesh();
		std::vector<double> getCoords();
		void triangulate();
		double get_hull_area();
		std::vector<double> get_hull_coords();
		std::vector<size_t> get_hull_points();

		double edge_length(size_t e);
		size_t get_interior_point(size_t e);
		std::vector<double> concavehull(double chi_factor = 0.1);
		std::vector<std::size_t> get_triangles();
		std::vector<double> get_tri_coordinates();

	public://variables
		std::vector<double> coords;
		std::vector<std::size_t> triangles;
		std::vector<std::size_t> halfedges;
		std::vector<std::size_t> hull_prev;
		std::vector<std::size_t> hull_next;
		std::vector<std::size_t> hull_tri;
		std::size_t hull_start;
	public://methods
		std::size_t legalize(std::size_t a);
		std::size_t hash_key(double x, double y) const;
		std::size_t add_triangle(
			std::size_t i0,
			std::size_t i1,
			std::size_t i2,
			std::size_t a,
			std::size_t b,
			std::size_t c);
		void link(std::size_t a, std::size_t b);
		size_t  next_halfedge(size_t e);

		size_t prev_halfedge(size_t e);
		std::pair<Segment, Segment> findConvexHullSegments(const std::vector<double>&, const std::vector<double>&);

		bool isOnPolygon(const Point&, const Point&, const Point&);
		bool oneIsOnPolygon(const Point&);
		bool is_left(const Point& a, const Point& b, const Point& c);
		bool isOutsideConcaveHull(const Point&);
		//void makeISCO(const std::vector<double>& , const std::vector<double>& );
		void makeISCO();
		void make_driehoeken();

	public://variables
		std::vector<std::size_t> m_hash;
		double m_center_x;
		double m_center_y;
		std::size_t m_hash_size;
		std::vector<std::size_t> m_edge_stack;
		std::vector<double>  x_coords;
		std::vector<double>  y_coords;
		std::vector<double> redshift_field;
		std::vector<driehoek> driehoeken;
		std::vector<Point> concaveHull;
		int nConcaveHull;
		std::vector<double>  xISCO;
		std::vector<double>  yISCO;
		double xIscoMax;// = *std::max_element(xisco.begin(), xisco.end());
		double xIscoMin;// = *std::min_element(xisco.begin(), xisco.end());
		double yIscoMax;// = *std::max_element(xisco.begin(), xisco.end());
		double yIscoMin;// = *std::min_element(xisco.begin(), xisco.end());
		double IscoEpsilon=0.0;// = 0;
		std::vector<Point> ISCO;
	};

} //namespace meshes
#endif