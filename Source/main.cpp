//Unoffial C++ version of Luminet project of bgmeulem
//https://github.com/bgmeulem/Luminet/


#include "BlackHole.h"
#include "BlackHolePhysics.h"
#include <discpp.h>

void plot_IR(std::vector<double>&, std::vector<double>&);
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
	Isoradial IR(30.0*bh_mass,80*M_PI/180, bh_mass,0);
	IR.calculate();
	std::pair<std::vector<double>, std::vector<double>> bare_isoradials= IR.get_bare_isoradials();//temporary functions for debugging
	std::vector<double> xx = std::get<0>(bare_isoradials);
	std::vector<double> yy = std::get<1>(bare_isoradials);
	/*for (int i = 0; i < xx.size(); i++) {
		std::cout << i << "):  (x, y) = (" << xx[i] << ", " << yy[i] << ")" << std::endl;
	}*/
	plot_IR(yy, xx);

	return 0;
}

void plot_IR(std::vector<double>& xx, std::vector<double>& yy) {
	Dislin g;
	int n = 400;
	double xray[400], yray[400], x2[10], y2[10];

	double x_max = -10000000000.0;
	double x_min = 10000000000.0;
	double y_max = -10000000000.0;
	double y_min = 10000000000.0;
	for (unsigned i = 0; i < n; i++)
	{
		yray[i] = xx[i];
		xray[i] = yy[i];
		if (xray[i] <= x_min)x_min = xray[i];
		if (xray[i] >= x_max)x_max = xray[i];
		if (yray[i] <= y_min)y_min = yray[i];
		if (yray[i] > y_max)y_max = yray[i];
	}

	g.metafl("cons");
	g.scrmod("revers");
	g.disini();
	g.pagera();
	g.complx();
	g.axspos(450, 1800);
	g.axslen(2200, 1200);

	g.name("X-axis", "x");
	g.name("Y-axis", "y");

	g.labdig(-1, "x");
	g.ticks(9, "x");
	g.ticks(10, "y");

	g.titlin("Test of a bare Isoradial", 1);
	g.titlin("(no redshift)", 3);


	int ic = g.intrgb(0.95, 0.95, 0.95);
	g.axsbgd(ic);

	g.graf(x_min, x_max, x_min, (x_max-x_min)/10, y_min, y_max, y_min,  (y_max - y_min) / 10);
	g.setrgb(0.7, 0.7, 0.7);
	g.grid(1, 1);
	g.winfnt("Times New Roman Bold");
	//g.psfont("Courier-Bold");
	g.color("fore");
	g.height(50);
	g.title();

	g.color("red");
	g.linwid(3);
	g.curve(xray, yray, n);

	g.disfin();
}