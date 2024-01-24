//Unoffial C++ version of Luminet project of bgmeulem
//https://github.com/bgmeulem/Luminet/

#include "BlackHole.h"
#include "BlackHolePhysics.h"
#include "plotter.h"
#include <conio.h>
#include <iostream>
//#include <discpp.h>
/*
************************************************************************************************************************************
The plotting has been moved to the Hary Plotter class.
************************************************************************************************************************************
*/
/*void plot_IR(std::vector<double>&, std::vector<double>&, std::vector<double>&, std::vector<double>&);//forward declaration Isoradials
void plot_RS(std::vector<double>&, std::vector<double>&, std::vector<double>&, std::vector<double>&, std::vector<double>&, std::vector<double>&);

// Function to convert a value to RGB format
std::vector<std::tuple<double, double, double> >convertToRGB(std::vector<double>);
std::vector<double> normalize_vector(std::vector<double>);*/
int main() {
	BHphysics BHp;

	// Test the translated functions
	double periastron = 30.;// provide periastron value ;
	double bh_mass = 1.;// provide black hole mass value
	double radius = 100.0 * bh_mass;// 20. provide radius value
	double angle = 20.;//provide angle value ;
	double incl = 76.;// provide inclination value ;
	double acc = 5.;// provide accretion value
	double redshift_factor = 2.;// provide redshift factor ;
	double b_ = 1.0;// provide b_ value ;

	/*double result_phi_inf = BHp.phi_inf(periastron, bh_mass);

	double result_mu = BHp.mu(periastron, bh_mass);

	double result_ellipse = BHp.ellipse(radius, angle, incl);

	double result_flux_intrinsic = BHp.flux_intrinsic(radius, acc, bh_mass);

	double result_flux_observed = BHp.flux_observed(radius, acc, bh_mass, redshift_factor);

	double result_redshift_factor = BHp.redshift_factor(radius, angle, incl, bh_mass, b_);

	std::cout << "Result phi_inf: " << result_phi_inf << std::endl;
	std::cout << "Result mu: " << result_mu << std::endl;
	std::cout << "Result ellipse: " << result_ellipse << std::endl;
	std::cout << "Result flux_intrinsic: " << result_flux_intrinsic << std::endl;
	std::cout << "Result flux_observed: " << result_flux_observed << std::endl;
	std::cout << "Result result_redshift_factor: " << result_redshift_factor << std::endl;

	//Test data reading
	BlackHole BH(bh_mass);
	BH.sample_Sources(20000);// , "D:/dev/Cpp/Luminet/Points/points_incl=85b.csv", "D:/dev/Cpp/Luminet/Points/points_secondary_incl=85b.csv");
	*/
	//test single isoradial
	double inclination = 80.0;//in degrees
	double distance = 30.0;
	Isoradial* IR;
	Isoradial* IRg;//ghost

	//setting up the plotting
	std::pair<std::vector<double>, std::vector<double>> bare_isoradials ;//temporary functions for debugging
	std::vector<double> xx;
	std::vector<double> yy;
	std::pair<std::vector<double>, std::vector<double>> bare_ghost_isoradials;//temporary functions for debugging
	std::vector<double> xx_ghost ;
	std::vector<double> yy_ghost;
	std::vector<double> redshift_factors;
	std::vector<double> redshift_factors_ghost ;
	Plotter plotter;
	char escape;
	std::cout << "press esc to exit! " << std::endl;
	bool loop = true;
	std::vector<double> inclinations;
	for (double p = 10; p < 81.0; p += 1.0) {
		inclinations.push_back(p);
	}
	size_t s = inclinations.size();
	/*for (size_t i = 0; i < s; i++) {
		inclinations.push_back(inclinations[s - i - 1]);
	}*/
	int sign = 1;
	int count = 0;
	unsigned index = 0;
	double theta = 0.0;
	while (true)
	{

	//inclination = { 80.0 };
	//for (auto theta : inclinations) {
		theta = inclinations[count];
		IR = new Isoradial (distance * bh_mass, theta * M_PI / 180, bh_mass, 0);
		IR->calculate();
		IRg= new Isoradial(distance * bh_mass, theta* M_PI / 180, bh_mass, 1);//ghost
		IRg->calculate();
		bare_isoradials = IR->get_bare_isoradials();//temporary functions for debugging
		xx = std::get<0>(bare_isoradials);
		yy = std::get<1>(bare_isoradials);
		bare_ghost_isoradials = IRg->get_bare_isoradials();//temporary functions for debugging
		xx_ghost = std::get<0>(bare_ghost_isoradials);
		yy_ghost = std::get<1>(bare_ghost_isoradials);
		redshift_factors = IR->get_redshift_factors();
		redshift_factors_ghost = IRg->get_redshift_factors();

		//plotter.plot(xx, yy, xx_ghost,yy_ghost);
		plotter.plot(theta, xx, yy, xx_ghost, yy_ghost, redshift_factors, redshift_factors_ghost, loop);
		count = count + sign;
		if (count == -1) {
			sign = 1;
			count = 0;
		}
		if (count == inclinations.size()) {
			sign = -1;
			count = inclinations.size()-1;
		}
		std::cout << count << std::endl;
	}
	/*escape = _getch();
	if (escape == 27)
		break;
			std::cout << "exited: " << std::endl;
	}*/


	return 0;
}
/*
void plot_IR(std::vector<double>& xx, std::vector<double>& yy, std::vector<double>& xx_g, std::vector<double>& yy_g) {
	Dislin g;
	int n = xx.size();

	for (unsigned i = 0; i < n; i++)
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

	g.metafl("cons");
	//g.scrmod("revers");
	g.disini();
	g.pagera();
	g.complx();
	g.axspos(450, 1800);
	g.axslen(1200, 1200);

	g.name("X-axis", "x");
	g.name("Y-axis", "y");

	g.labdig(-1, "x");
	g.ticks(9, "x");
	g.ticks(10, "y");

	g.titlin("Test of a bare Isoradial", 1);
	g.titlin("(no redshift)", 3);

	//int ic = g.intrgb(0.95, 0.95, 0.95);
	int ic = g.intrgb(0., 0., 0.);
	g.axsbgd(ic);

	g.graf(x_min, x_max, x_min, (x_max - x_min) / 10, y_min, y_max, y_min, (y_max - y_min) / 10);
	//g.setrgb(0.7, 0.7, 0.7);
	g.setrgb(0.0, 0.0, 0.0);
	g.grid(1, 1);
	g.winfnt("Times New Roman Bold");
	//g.psfont("Courier-Bold");
	g.color("fore");
	g.height(50);
	g.title();

	g.color("red");
	g.linwid(3);
	//g.curve(xray, yray, n);
	g.curve(&xx[0], &yy[0], n);
	g.color("blue");
	g.linwid(3);
	g.curve(&xx_g[0], &yy_g[0], n);

	g.disfin();
}

void plot_RS(std::vector<double>& xx, std::vector<double>& yy, std::vector<double>& xx_g, std::vector<double>& yy_g, std::vector<double>& rs, std::vector<double>& rs_g) {
	Dislin g;
	int n = xx.size();

	// Vector to store RGB color values
	std::vector<std::tuple<double, double, double>> rgbVector = convertToRGB(normalize_vector(rs));
	std::vector<std::tuple<double, double, double>> rgbVector_g = convertToRGB(normalize_vector(rs_g));

	double x_max = -10000000000.0;
	double x_min = 10000000000.0;
	double y_max = -10000000000.0;
	double y_min = 10000000000.0;

	for (unsigned i = 0; i < n; i++)
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

	g.metafl("cons");
	//g.scrmod("revers");
	g.disini();
	g.pagera();
	g.complx();
	g.axspos(450, 1800);
	g.axslen(1200, 1200);

	g.name("X-axis", "x");
	g.name("Y-axis", "y");

	g.labdig(-1, "x");
	g.ticks(9, "x");
	g.ticks(10, "y");

	g.titlin("Test of a bare Isoradial", 1);
	g.titlin("(with redshift)", 3);

	//int ic = g.intrgb(0.95, 0.95, 0.95);
	int ic = g.intrgb(0., 0., 0.);
	g.axsbgd(ic);

	g.graf(x_min, x_max, x_min, (x_max - x_min) / 10, y_min, y_max, y_min, (y_max - y_min) / 10);
	//g.setrgb(0.7, 0.7, 0.7);
	g.setrgb(0.0, 0.0, 0.0);
	g.grid(1, 1);
	g.winfnt("Times New Roman Bold");
	//g.psfont("Courier-Bold");
	g.color("fore");
	g.height(50);
	g.title();
	// Set the color map based on the color parameter
	for (size_t i = 0; i < n; i++) {
		g.setrgb(std::get<0>(rgbVector[i]), std::get<1>(rgbVector[i]), std::get<2>(rgbVector[i]));
		g.rlcirc(xx[i], yy[i], 0.1);
		g.setrgb(std::get<0>(rgbVector_g[i]), std::get<1>(rgbVector_g[i]), std::get<2>(rgbVector_g[i]));
		g.rlcirc(xx_g[i], yy_g[i], 0.1);
	}

	g.disfin();
}

// Function to normalize the redshift in a range 0 to 100
//stretches it by taking the span of min and max
std::vector<double> normalize_vector(std::vector<double> avector) {
	std::vector<double> normalized_vector;
	double max_rs = *std::max_element(std::begin(avector), std::end(avector));
	double min_rs = *std::min_element(std::begin(avector), std::end(avector));
	for (auto value : avector) {
		value = 100*(value-min_rs)/ (max_rs-min_rs);
		//std::cout << value << std::endl;
		normalized_vector.push_back(value);
	}
	return normalized_vector;
}

// Function to convert a value to RGB format
std::vector<std::tuple<double, double, double> >convertToRGB(std::vector<double> avector) {
	std::vector<std::tuple<double, double, double> > colors;
	for (size_t i = 0; i < avector.size(); i++) {
		// Map the value to the RGB range
		double red = (avector[i]) / 100.;
		double green = (100 - avector[i]) / 100;
		double blue = 0.0;
		colors.push_back(std::make_tuple(red, green, blue));
	}
	return colors;
}*/