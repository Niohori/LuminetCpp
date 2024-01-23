//Unoffial C++ version of Luminet project of bgmeulem
//https://github.com/bgmeulem/Luminet/


#include "BlackHole.h"
#include "BlackHolePhysics.h"
#include <discpp.h>

void plot_IR(std::vector<double>&, std::vector<double>&, std::vector<double>&, std::vector<double>&);
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
	double inclination =80.0;//in degrees
	double distance = 30.0;
	Isoradial IR(distance *bh_mass, inclination *M_PI/180, bh_mass,0);
	IR.calculate();
	Isoradial IRg(distance * bh_mass, inclination * M_PI / 180, bh_mass, 1);//ghost
	IRg.calculate();
	std::pair<std::vector<double>, std::vector<double>> bare_isoradials= IR.get_bare_isoradials();//temporary functions for debugging
	std::vector<double> xx = std::get<0>(bare_isoradials);
	std::vector<double> yy = std::get<1>(bare_isoradials);
	std::pair<std::vector<double>, std::vector<double>> bare_ghost_isoradials = IRg.get_bare_isoradials();//temporary functions for debugging
	std::vector<double> xx_ghost = std::get<0>(bare_ghost_isoradials);
	std::vector<double> yy_ghost = std::get<1>(bare_ghost_isoradials);
	/*std::cout << std::endl << std::endl << std::endl << "Passed to main with size " << xx.size() << std::endl;
	for (size_t i = 0; i < xx.size(); i++) {
		std::cout << i << "):  (x, y) = (" << xx[i] << ", " << yy[i] << ")" << std::endl;
	}*/
	plot_IR(xx, yy, xx_ghost,yy_ghost);

	return 0;
}

void plot_IR(std::vector<double>& xx, std::vector<double>& yy, std::vector<double>& xx_g, std::vector<double>& yy_g) {
	Dislin g;
	int n =xx.size();
	
	//double xray[400], yray[400];

	double x_max = -10000000000.0;
	double x_min = 10000000000.0;
	double y_max = -10000000000.0;
	double y_min = 10000000000.0;
	/*
	
	int n = 400;
	for (unsigned i = 0; i < n; i++)
	{
		yray[i] = xx[i];
		xray[i] = yy[i];
		if (xray[i] <= x_min)x_min = xray[i];
		if (xray[i] >= x_max)x_max = xray[i];
		if (yray[i] <= y_min)y_min = yray[i];
		if (yray[i] > y_max)y_max = yray[i];
	}*/
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

	g.graf(x_min, x_max, x_min, (x_max-x_min)/10, y_min, y_max, y_min,  (y_max - y_min) / 10);
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