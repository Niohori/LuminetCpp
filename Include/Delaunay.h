#pragma once
#ifndef DELAUNAY_H
#define DELAUNAY_H
#include <algorithm>
#include <cmath>
#include <exception>
#include <iostream>
#include <limits>
#include <memory>
#include <utility>
#include <vector>

#include <algorithm>
#include <iostream>
#include <limits>
#include <set>
#include <stdexcept>
#include <vector>
#define DIFFERENCE 0.0005
#define EQ(_x_,_y_) (((_x_-_y_<DIFFERENCE)&&(_y_-_x_<DIFFERENCE))?1:0)

namespace delaunay {
	// Point structure representing a 2D Point
	struct Point {
		double x, y;
		Point(double x_, double y_) : x(x_), y(y_) {}
	};
	bool operator ==(Point, Point);
	bool operator !=(Point, Point);

	struct driehoek {
		//driehoek() {};
		size_t ID;//identification -> not necessary for the moment being
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
		};
	};
	struct Segment {
		Point p1;
		Point p2;
		size_t ID;//identification of Delaunay triangle
		Segment();
		Segment(size_t ID_, Point p1_, Point p2_) :
			ID(ID_), p1(p1_), p2(p2_)
		{
			;
		}
	};
	bool operator ==(Segment, Segment);
	bool operator !=(Segment, Segment);

	//======================================================================================================================================
	/*size_t next_halfedge(size_t e) {
		return (e % 3 == 2) ? e - 2 : e + 1;
	}

	size_t prev_halfedge(size_t e) {
		return (e % 3 == 0) ? e + 2 : e - 1;
	}
	*/
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

	struct DelaunayPoint {
		std::size_t i;
		double x;
		double y;
		std::size_t t;
		std::size_t prev;
		std::size_t next;
		bool removed;
	};

	class Delaunay {
	public:

		Delaunay(std::vector<double> const&, std::vector<double> const&, std::vector<double> const&);
		void triangulate();
		std::vector<Segment>  get_iso_lines(const double&, const std::vector<double>&, const std::vector<double>&);

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
	private://methods
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
		std::vector<Segment> generateContourLines(const std::vector<driehoek>&, double);
		double interpolate(const Point& p1, const Point& p2, double z1, double z2, const Point& query);
		std::pair<Segment, Segment> findConvexHullSegments(const std::vector<double>&, const std::vector<double>&);

		bool isOnPolygon(const double&, const Point&, const Point&, const Point&, const std::vector<double>&, const std::vector<double>&, const double&, const double&);
		bool is_left(const Point& a, const Point& b, const Point& c);
		bool isOutsideConcaveHull(const Point&);
		//bool isOutsideConcaveHull(const Point&, const std::vector<Point>&);
		//bool is_point_inside_hull(const Point&, const std::vector<Point>&);

		bool isNoise(Segment);

	private://variables
		std::vector<std::size_t> m_hash;
		double m_center_x;
		double m_center_y;
		std::size_t m_hash_size;
		std::vector<std::size_t> m_edge_stack;
		std::vector<double> const& x_coords;
		std::vector<double>  y_coords;
		std::vector<double> redshift_field;
		std::vector<driehoek> driehoeken;
		std::vector<Point> concaveHull;
		int nConcaveHull;
	};
} //namespace delaunay
#endif