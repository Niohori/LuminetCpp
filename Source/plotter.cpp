#include "Plotter.h"
Plotter::Plotter()
{
	screenWidth = GetSystemMetrics(SM_CXSCREEN);;
	screenHeight = GetSystemMetrics(SM_CYSCREEN);

	std::cout << "Screen Resolution: " << screenWidth << "x" << screenHeight << std::endl;

	g = new Dislin;
	x_max = -10000000000.0;
	x_min = 10000000000.0;
	y_max = -10000000000.0;
	y_min = 10000000000.0;
	Npoints = 0;
	g->metafl("cons");
	//g->scrmod("revers");
	g->disini();
	g->pagera();
	g->complx();
	g->nochek();
	g->axspos(500, 0 * screenHeight - 200);
	//g->axspos(1250, 1600);
	int scrheight = std::min(screenWidth, screenHeight) - 20;
	g->axslen(scrheight, scrheight);
	//g->axslen(1200, 2800);
	g->chacod("ISO1");
	g->winfnt("Times New Roman");
	ic = g->intrgb(0., 0., 0.);
	g->axsbgd(ic);
	//g->setclr(0);
}
Plotter::~Plotter() {
	// Destructor implementation
	g->disfin();
}

void Plotter::plot_iso_redshifts(const double& bhMass, const double& inclination, const std::multimap<double, std::vector<meshes::Point > >& isolines_, const double& x_max_, const double& x_min_, const double& max_rs, const double& min_rs, const std::pair<std::vector<double>, std::vector<double>>& isco_, const std::pair<std::vector<double>, std::vector<double>>& concavehull, const bool& loop) {
	//Dislin g;

	//==================== Creating the concave hull ===========================================================
	std::multimap<double, std::vector<meshes::Point > > isolines = isolines_;
	double factor = 2.0;
	for (auto& iso : isolines) {
		for (int p = 0; p < iso.second.size(); p++) {
			iso.second[p].x /= (x_max_ / factor);
			iso.second[p].y /= (x_max_ / factor);
		}
	}
	std::pair<std::vector<double>, std::vector<double>> isco = isco_;
	for (int i = 0; i < isco.first.size(); i++) {
		isco.first[i] /= (x_max_ / factor);
		isco.second[i] /= (x_max_ / factor);
	}
	double N_redshifts = isolines.size();
	// Counting the number of different keys
	std::vector<double> uniqueKeys;
	int nIsolines = 0;
	for (const auto& pair : isolines) {
		uniqueKeys.push_back(pair.first);
		nIsolines++;
	}
	std::cout << "Number of Isolines = " << nIsolines << std::endl;
	rgbVector = convertToRGB(normalize_vector(uniqueKeys));
	//rgbVector = convertToRGB(uniqueKeys);
	//x_max = std::min(x_max_, 60.0);

	x_max = 1.1;
	x_min = -x_max;
	y_max = x_max;
	y_min = -x_max;
	if (loop) {
		g->erase();
	};
	g->titlin("Redshift isolines", 1);
	std::string bhmass_as_string = std::to_string(int(bhMass));
	std::string tit1 = "Mass  = " + bhmass_as_string + "Msun";//static_cast<char>('\u00B0');
	g->titlin(&tit1[0], 2);
	// Convert double to string with 2-digit precision
	std::string inclination_as_string = std::to_string(inclination);
	size_t dotPos = inclination_as_string.find('.');
	if (dotPos != std::string::npos && dotPos + 0 < inclination_as_string.size()) {
		inclination_as_string = inclination_as_string.substr(0, dotPos + 0); // keep 2 digits after the dot
	}
	if (inclination_as_string.size() == 1) inclination_as_string = "00" + inclination_as_string;
	if (inclination_as_string.size() == 2) inclination_as_string = "0" + inclination_as_string;
	std::string tittel = "inclination = " + inclination_as_string + static_cast<char>(186);//static_cast<char>('\u00B0');
	g->titlin(&tittel[0], 3);
	g->setclr(0);

	g->graf(x_min, x_max, x_min, (x_max - x_min) / 10, y_min, y_max, y_min, (y_max - y_min) / 10);
	g->setrgb(0.0, 0.0, 0.0);
	g->color("fore");
	g->title();
	// Set the color map based on the color parameter

	g->linwid(3);
	double color = 0.5;

	auto clr = g->intrgb(0.0, 0.0, 0.0);
	double color_i = 1.0 / isolines.size();
	size_t cc = 0;
	for (auto iso : isolines) {
		auto segments = iso.second;
		auto cl = rgbVector[cc];
		g->setrgb(std::get<0>(cl), std::get<1>(cl), std::get<2>(cl));
		//g->setrgb(cc * color_i, 0, 1.0 - cc * color_i);
		/*	for (auto segment : segments) {
				g->rline(segment.p1.x, segment.p1.y, segment.p2.x, segment.p2.y);
			}*/
		if (segments.size() != 0) {
			for (auto i = 0; i < segments.size() - 1; i++) {
				if (std::pow(segments[i].x - segments[i + 1].x, 2) + std::pow(segments[i].y - segments[i + 1].y, 2) > 0.002) { continue; }
				//if (segments[i].y*segments[i + 1].y<0 ) { continue; }
				g->rline(segments[i].x, segments[i].y, segments[i + 1].x, segments[i + 1].y);
			}
		}
		cc++;
	}

	g->setrgb(0.1, 0.1, 0.1);
	g->lintyp(0);
	//===================================== Plotting the concave hull  ===========================================================================
	//g->hsymbl(1);
	//for (size_t i = 0; i < concavehull.first.size(); i++) {
	//	g->rlsymb(21, concavehull.first[i], concavehull.second[i]);
	//}

	//===================================== Plotting the Isco  ===========================================================================
	//g->setrgb(0.0, 1.0, 0.0);
	g->curve(&isco.first[0], &isco.second[0], isco.first.size());

	if (loop) {
		g->endgrf();
		g->sendbf();
	};
	//g.disfin();
}

void Plotter::plot_isoradials(std::vector<double>& xx, std::vector<double>& yy, std::vector<double>& xx_g, std::vector<double>& yy_g) {
	//Dislin g;
	Npoints = xx.size();

	for (unsigned i = 0; i < Npoints; i++)
	{
		if (xx[i] <= x_min)x_min = xx[i];
		if (xx[i] >= x_max)x_max = xx[i];
		if (yy[i] <= y_min)y_min = yy[i];
		if (yy[i] > y_max)y_max = yy[i];
	}
	x_max *= 1.1;
	x_min *= 1.1;
	y_max *= 1.1;
	y_min *= 1.1;
	if (x_max < y_max)x_max = y_max;
	if (x_min < y_min)y_min = x_min;
	if (x_max > y_max)y_max = x_max;
	if (x_min > y_min)x_min = y_min;

	g->metafl("cons");
	//g->scrmod("revers");
	g->disini();
	g->pagera();
	g->complx();
	g->axspos(450, 1800);
	g->axslen(1200, 1200);

	g->name("X-axis", "x");
	g->name("Y-axis", "y");

	g->labdig(-1, "x");
	g->ticks(9, "x");
	g->ticks(10, "y");

	g->titlin("Test of a bare Isoradial", 1);
	g->titlin("(no redshift)", 3);

	//int ic = g->intrgb(0.95, 0.95, 0.95);
	int ic = g->intrgb(0., 0., 0.);
	g->axsbgd(ic);

	g->graf(x_min, x_max, x_min, (x_max - x_min) / 10, y_min, y_max, y_min, (y_max - y_min) / 10);
	//g->setrgb(0.7, 0.7, 0.7);
	g->setrgb(0.0, 0.0, 0.0);
	g->grid(1, 1);
	g->winfnt("Times New Roman Bold");
	//g->psfont("Courier-Bold");
	g->color("fore");
	g->height(50);
	g->title();

	g->color("red");
	g->linwid(3);
	g->curve(&xx[0], &yy[0], Npoints);
	g->color("blue");
	g->linwid(3);
	g->curve(&xx_g[0], &yy_g[0], Npoints);

	//g->disfin();
}

void Plotter::plot_isoradials(double inclination, std::vector<double>& xx_, std::vector<double>& yy_, std::vector<double>& xx_g_, std::vector<double>& yy_g_, std::vector<double>& rs, std::vector<double>& rs_g, const double& Xmax, bool loop) {
	//Dislin g;

	std::vector<double> xx = xx_;
	std::vector<double> yy = yy_;
	std::vector<double> xx_g = xx_g_;
	std::vector<double> yy_g = yy_g_;
	Npoints = xx.size();
	rgbVector = convertToRGB(normalize_vector(rs));
	rgbVector_g = convertToRGB(normalize_vector(rs_g));
	//rgbVector = convertToRGBbis(rs, 5000.0);
	//rgbVector_g = convertToRGBbis(rs_g, 5000.0);
	double factor = 1;
	for (unsigned i = 0; i < Npoints; i++)
	{
		if (Xmax == 0.0)break;
		xx[i] /= (Xmax / factor);
		yy[i] /= (Xmax / factor);
		xx_g[i] /= (Xmax / factor);
		yy_g[i] /= (Xmax / factor);
	}
	//g->ticks(10, "y");
	if (loop) {
		g->erase();
	};
	x_max = 1.05;
	x_min = -x_max;
	y_max = x_max;
	y_min = x_min;
	// Convert double to string with 2-digit precision
	std::string inclination_as_string = std::to_string(inclination);
	size_t dotPos = inclination_as_string.find('.');
	if (dotPos != std::string::npos && dotPos + 0 < inclination_as_string.size()) {
		inclination_as_string = inclination_as_string.substr(0, dotPos + 0); // keep 2 digits after the dot
	}
	if (inclination_as_string.size() == 1) inclination_as_string = "00" + inclination_as_string;
	if (inclination_as_string.size() == 2) inclination_as_string = "0" + inclination_as_string;
	std::string tittel = "inclination = " + inclination_as_string + static_cast<char>(186);//static_cast<char>('\u00B0');
	g->titlin(&tittel[0], 3);
	g->setclr(0);

	g->graf(x_min, x_max, x_min, (x_max - x_min) / 10, y_min, y_max, y_min, (y_max - y_min) / 10);
	//g->setrgb(0.7, 0.7, 0.7);
	g->setrgb(0.0, 0.0, 0.0);
	g->color("fore");
	g->title();
	// Set the color map based on the color parameter
	size_t cc = 0;
	for (size_t i = 0; i < Npoints; i++) {
		g->setrgb(std::get<0>(rgbVector[i]), std::get<1>(rgbVector[i]), std::get<2>(rgbVector[i]));
		g->rlcirc(xx[i], yy[i], 0.005);
		g->setrgb(std::get<0>(rgbVector_g[i]), std::get<1>(rgbVector_g[i]), std::get<2>(rgbVector_g[i]));
		g->rlcirc(xx_g[i], yy_g[i], 0.005);
	}
	if (loop) {
		g->endgrf();
		g->sendbf();
	};
	//g->disfin();
}

// Function to normalize the redshiftfactor in a range 0 to 100
//stretches it by taking the span of min and max
std::vector<double> Plotter::normalize_vector(std::vector<double> avector) {
	std::vector<double> normalized_vector;
	double max_rs = *std::max_element(std::begin(avector), std::end(avector));
	double min_rs = *std::min_element(std::begin(avector), std::end(avector));
	for (auto value : avector) {
		value = 100 * (value - min_rs) / (max_rs - min_rs);
		//std::cout << value << std::endl;
		normalized_vector.push_back(value);
	}
	return normalized_vector;
}

// Function to convert a value to RGB format
std::vector<std::tuple<double, double, double> > Plotter::convertToRGB(const std::vector<double>& avector) {
	std::vector<std::tuple<double, double, double> > colors;
	double kelvin = 50000;
	for (auto value : avector) {
		if (value == 0.0) {
			colors.push_back(std::make_tuple(0, 0, 0));
			continue;
		}
		std::vector<double> rgb = BHphysics::kelvinToRGB(kelvin / value);
		// Map the value to the RGB range
		double red = rgb[0];
		double green = rgb[1];
		double blue = rgb[2];
		colors.push_back(std::make_tuple(red, green, blue));
	}
	return colors;
}

// Function to convert a value to RGB format
std::vector<std::tuple<double, double, double> > Plotter::convertToRGBbis(const std::vector<double>& avector, const double& temperature) {
	std::vector<std::tuple<double, double, double> > colors;
	for (auto value : avector) {
		std::vector<double> rgb = BHphysics::wavelengthToRGB(temperature / value, 100.0);
		// Map the value to the RGB range
		double red = rgb[0];
		double green = rgb[1];
		double blue = rgb[2];
		colors.push_back(std::make_tuple(red, green, blue));
	}
	return colors;
}

std::vector<double> Plotter::flatten_matrix(const std::vector<std::vector<double> >& amatrix) {
	std::vector<double> zarr;
	for (const auto& row : amatrix) {
		zarr.insert(zarr.end(), row.begin(), row.end());
	}

	return zarr;
}

void Plotter::plot_BlackHole(double inclination, std::vector<double>& xx_, std::vector<double>& yy_, std::vector<double>& fluxes, std::vector<double>& xx_S, std::vector<double>& yy_S, std::vector<double>& fluxesS, const bool& loop) {
	//Dislin g;
	std::vector<double> x_shad;
	std::vector<double> y_shad;
	for (double al = 0; al < 2.0 * 3.14159; al = al + 0.01) {
		double r = 5.19695;
		x_shad.push_back(r * std::cos(al));
		y_shad.push_back(r * std::sin(al));
	}
	std::vector<double> xx = xx_;
	std::vector<double> yy = yy_;

	std::vector<double> xxS = xx_S;
	std::vector<double> yyS = yy_S;

	Npoints = xx.size();
	double x_max = -1e8;
	double x_min = 1e8;
	double y_max = -1e8;
	double y_min = 1e8;
	double factor = 1;
	for (unsigned i = 0; i < Npoints; i++)
	{
		if (x_max < xx[i]) { x_max = xx[i]; }
	}
	double magnifier = 1.5;
	x_max /= magnifier;
	x_min = -x_max;
	y_max = x_max;
	y_min = x_min;
	//for (unsigned i = 0; i < Npoints; i++)
	//{
	//	xx[i] /= x_max;
	//	yy[i] /= x_max;
	//	xxS[i] /= x_max;
	//	yyS[i] /= x_max;
	//}
	//for (int i = 0; i < x_shad.size(); i++) {
	//	x_shad[i] /= x_max;;
	//	y_shad[i] /= x_max;;

	//}
	//x_max = 1.;
	//x_min = -x_max;
	//y_max = x_max;
	//y_min = x_min;
	if (loop) {
		g->erase();
	};
	// Convert double to string with 2-digit precision
	//std::string inclination_as_string = std::to_string(inclination);
	//size_t dotPos = inclination_as_string.find('.');
	//if (dotPos != std::string::npos && dotPos + 0 < inclination_as_string.size()) {
	//	inclination_as_string = inclination_as_string.substr(0, dotPos + 0); // keep 2 digits after the dot
	//}
	//if (inclination_as_string.size() == 1) inclination_as_string = "00" + inclination_as_string;
	//if (inclination_as_string.size() == 2) inclination_as_string = "0" + inclination_as_string;
	//std::string tittel = "inclination = " + inclination_as_string + static_cast<char>(186);//static_cast<char>('\u00B0');
	//g->titlin(&tittel[0], 3);
	g->setclr(0);

	g->graf(x_min, x_max, x_min, (x_max - x_min) / 10, y_min, y_max, y_min, (y_max - y_min) / 10);
	//g->setrgb(0.7, 0.7, 0.7);
	g->setrgb(0.0, 0.0, 0.0);
	g->color("fore");
	g->title();
	//g->hsymbl(5);
	//for (size_t i = 0; i < x_shad.size(); i++) {
	//	;
	//	g->setrgb(0, 1, 0);
	//	//std::cout << x_shad[i] << " , " << y_shad[i] << std::endl;
	//	g->rlsymb(21, x_shad[i], y_shad[i]);
	//}
	g->hsymbl(4);
	// Set the color map based on the color parameter
	for (size_t i = 0; i < Npoints; i++) {
		double c = fluxesS[i];
		g->setrgb(c, c, c);
		//g->setrgb(0, 1, 0);
		g->rlsymb(21, xxS[i], yyS[i]);
	}

	for (size_t i = 0; i < Npoints; i++) {
		double c = fluxes[i];
		g->setrgb(c, c, c);
		//g->rlcirc(xx[i], yy[i], 0.01);
		g->setrgb(1, 0, 0);
		g->rlsymb(21, xx[i], yy[i]);
	}
	if (loop) {
		g->endgrf();
		g->sendbf();
	};
}

void Plotter::plot_Eq13(const std::vector<double> & peri, const std::vector<double> & val) {
	Npoints = peri.size();
	double x_max = -1e8;
	double x_min = 1e8;
	double y_max = -1e8;
	double y_min = 1e8;
	double factor = 1;

		for (unsigned i = 0; i < Npoints; i++)
		{
			if (x_max < peri[i]) { x_max = peri[i]; }
			if (x_min > peri[i]) { x_min = peri[i]; }
			if (y_max < val[i]) { y_max = val[i]; }
			if (y_min > val[i]) { y_min = val[i]; }
		}
	
	x_max =140.0;
	x_min = 0.0;
	y_min = -20.0;
	y_max = 20.0;

	g->metafl("cons");
	g->scrmod("revers");
	g->disini();
	g->pagera();
	g->complx();
	g->axspos(450, 1800);
	g->axslen(1200, 1200);

	g->name("Periastrons", "x");
	g->name("Eq13 value", "y");

	g->labdig(-1, "x");
	g->ticks(9, "x");
	g->ticks(10, "y");
	ic = g->intrgb(0.95, 0.95, 0.95);
	g->axsbgd(ic);
	g->graf(x_min, x_max, x_min, (x_max - x_min) / 10, y_min, y_max, y_min, (y_max - y_min) / 10);
	g->setrgb(0, 0, 1.0);
	g->grid(1, 1);

	g->color("red");

		g->curve(&peri[0], &val[0], Npoints);

}