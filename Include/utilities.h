#pragma once
#ifndef UTILITIES_H
#define UTILITIES_H

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

const double M_PI = 3.14159265358979323846;
//INI settings are wrapped in struct
struct irs_solver_params {
	unsigned initial_guesses = 12;
	unsigned midpoint_iterations = 12;
	double	times_inbetween = 2;// amount of times to double the precision of an isoredshift line when improving
	unsigned retry_angular_precision = 15;//angular precision to calculate isoradials with when improving solutions
	double	min_periastron = 3.01;// minimum distance to black hole(must be strictly larger than 3M), in units of black hole mass(photon sphere is at 3M)
	bool use_ellipse = true;
	unsigned retry_tip = 50;
	unsigned initial_radial_precision = 15;
	bool plot_inbetween = false;// plot isoredshifts while improving them
	double angular_margin = 0.3;
};

struct angular_properties {
	double start_angle = 0.0;
	double end_angle = M_PI;
	unsigned angular_precision = 500;
	bool mirror = true;
};

struct ir_params {
	double start_angle = 0.0;
	double	end_angle = M_PI;
	unsigned angular_precision = 500;
	bool	mirror = true;// if True, calculates only half of the isoradial and mirrors it
	double angular_margin = 0.3;
};

struct plot_params {
	bool plot_isoredshifts_inbetween = false;
	bool save_plot = false;
	bool plot_ellipse = false;
	bool plot_core = true;
	bool redshift = true;
	std::string	linestyle = "-";
	double	linewidth = 1.;
	std::string key = "";
	std::string face_color = "black";
	std::string line_color = "white";
	std::string text_color = "white";
	double alpha = 1.;
	bool show_grid = false;
	bool legend = false;
	bool orig_background = false;
	bool plot_disk_edges = false;
	std::pair<double, double>	ax_lim = { -100, 100 };
	std::string title = "Isoradials for R =";
};

struct solver_params {
	unsigned initial_guesses = 12;
	unsigned midpoint_iterations = 20;
	bool plot_inbetween = false;// plot isoredshifts while improving them
	double min_periastron = 3.001;// minimum distance to black hole, in units of black hole mass(photon sphere is at 3M)
	bool use_ellipse = true;
};

struct  find_redshift_params {
	bool force_redshift_solution = false;
	unsigned max_force_iter = 5;
};

class CSVRow
{
public:
	std::string_view operator[](std::size_t index) const
	{
		return std::string_view(&m_line[m_data[index] + 1], m_data[index + 1] - (m_data[index] + 1));
	}
	std::size_t size() const
	{
		return m_data.size() - 1;
	}
	void readNextRow(const std::string& line)
		//void readNextRow(std::istream& str)
	{
		//std::getline(str, m_line);
		m_line = line;
		//std::cout << m_line << std::endl;
		m_data.clear();
		m_data.emplace_back(-1);
		std::string::size_type pos = 0;
		while ((pos = m_line.find(',', pos)) != std::string::npos)
		{
			m_data.emplace_back(pos);
			++pos;
		}
		// This checks for a trailing comma with no data after it.
		pos = m_line.size();
		m_data.emplace_back(pos);
	}
private:
	std::string         m_line;
	std::vector<int>    m_data;
};
/*
std::istream& operator>>(std::istream& str, CSVRow& data)
{
	data.readNextRow(str);
	return str;
}*/

#endif