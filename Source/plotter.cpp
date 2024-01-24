#include "Plotter.h"
Plotter::Plotter()
{
	g = new Dislin;
	x_max = -10000000000.0;
	x_min = 10000000000.0;
	y_max = -10000000000.0;
	y_min = 10000000000.0;
	Npoints = 0;
}
Plotter::~Plotter() {
	// Destructor implementation
	g->disfin();
}

void Plotter::plot(std::vector<double>& xx, std::vector<double>& yy, std::vector<double>& xx_g, std::vector<double>& yy_g) {
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

void Plotter::plot(double inclination,std::vector<double>& xx, std::vector<double>& yy, std::vector<double>& xx_g, std::vector<double>& yy_g, std::vector<double>& rs, std::vector<double>& rs_g) {
	//Dislin g;
	Npoints = xx.size();

	rgbVector = convertToRGB(normalize_vector(rs));
	rgbVector_g = convertToRGB(normalize_vector(rs_g));

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

	//g->name("X-axis", "x");
	//g->name("Y-axis", "y");

	//g->labdig(-1, "x");
	//g->ticks(9, "x");
	//g->ticks(10, "y");
	// Convert double to string with 2-digit precision
	std::string inclination_as_string  = std::to_string(inclination);
	size_t dotPos = inclination_as_string.find('.');
	if (dotPos != std::string::npos && dotPos + 0 < inclination_as_string.size()) {
		inclination_as_string = inclination_as_string.substr(0, dotPos + 0); // keep 2 digits after the dot
	}
	if (inclination_as_string.size() == 1) inclination_as_string = "00"+ inclination_as_string;
	if (inclination_as_string.size() == 2) inclination_as_string = "0" + inclination_as_string;
	char degreeSymbol = 176;
	std::string tittel = "inclination = " + inclination_as_string +  static_cast<char>(186);//static_cast<char>('\u00B0');
	g->titlin(&tittel[0], 3);
	//g->titlin("(with redshift)", 3);

	//int ic = g->intrgb(0.95, 0.95, 0.95);
	int ic = g->intrgb(0., 0., 0.);
	g->axsbgd(ic);
	g->setclr(0);

	g->graf(x_min, x_max, x_min, (x_max - x_min) / 10, y_min, y_max, y_min, (y_max - y_min) / 10);
	//g->setrgb(0.7, 0.7, 0.7);
	g->setrgb(0.0, 0.0, 0.0);
	//g->grid(1, 1);
	g->chacod("ISO1");
	g->winfnt("Times New Roman");
	//g->psfont("Courier-Bold");
	g->color("fore");
	g->height(50);
	g->title();
	// Set the color map based on the color parameter
	for (size_t i = 0; i < Npoints; i++) {
		g->setrgb(std::get<0>(rgbVector[i]), std::get<1>(rgbVector[i]), std::get<2>(rgbVector[i]));
		g->rlcirc(xx[i], yy[i], 0.1);
		g->setrgb(std::get<0>(rgbVector_g[i]), std::get<1>(rgbVector_g[i]), std::get<2>(rgbVector_g[i]));
		g->rlcirc(xx_g[i], yy_g[i], 0.1);
	}

	//g->disfin();
}

// Function to normalize the redshift in a range 0 to 100
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
std::vector<std::tuple<double, double, double> > Plotter::convertToRGB(std::vector<double> avector) {
	std::vector<std::tuple<double, double, double> > colors;
	for (auto value:avector) {
		// Map the value to the RGB range
		double red = (value) / 100.;
		double green = (100 - value) / 100;
		double blue = 0.0;
		colors.push_back(std::make_tuple(red, green, blue));
	}
	return colors;
}

/*

// Function to convert a value to RGB format
std::vector<std::tuple<double, double, double> > Plotter::convertToRGB(std::vector<double> avector) {
	std::vector<std::tuple<double, double, double> > colors;
	double normalizedValue;
	double hue=0.0;
	// Set saturation and value to 1 for full intensity and brightness
	double saturation = 1.0;
	double brightness = 1.0;
	// Convert HSV to RGB
	double C = brightness * saturation;
	double X = C * (1 - std::abs(std::fmod(hue / 60.0, 2) - 1));
	double m = brightness - C;
	double r =0.;
	double g = 0.;
	double b = 0.;
	for (auto value:avector) {
		normalizedValue = value / 100.0;
		// Map the normalized value to the hue in the HSV color model
		hue = 0.0 + normalizedValue * 120.0;  // 0 is red, 120 is green

		// Map the value to the RGB range
		if (hue >= 0 && hue < 60) {
			r = C;
			g = X;
			b = 0;
		}
		else if (hue >= 60 && hue < 120) {
			r = X;
			g = C;
			b = 0;
		}
		else if (hue >= 120 && hue < 180) {
			r = 0;
			g = C;
			b = X;
		}
		else if (hue >= 180 && hue < 240) {
			r = 0;
			g = X;
			b = C;
		}
		else if (hue >= 240 && hue < 300) {
			r = X;
			g = 0;
			b = C;
		}
		else {
			r = C;
			g = 0;
			b = X;
		}
		colors.push_back(std::make_tuple(r, g, b));
	}
	return colors;
}*/