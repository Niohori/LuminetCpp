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
TESTS.
************************************************************************************************************************************
*/
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
	/*
	* *********************************************************************************************************
	* The code for animation has been commented out for debugging purposes.
	* BUG to investigate: inclinations > 80° give broken isoradials 
	* * *********************************************************************************************************
	*/
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
	for (double p = 1; p < 180.0; p += 1.0) {
		inclinations.push_back(p);
	}
	//inclinations = { 5.0 }; loop = false;


	int sign = 1;
	int count = 0;
	unsigned index = 0;
	double theta = 0.0;
	while (true)
	{

	//inclination = { 80.0 };
	//for (auto theta : inclinations) {
		theta = inclinations[count];
		//std::cout << "theta = "<< theta * M_PI / 180 << std::endl;
		IR = new Isoradial (distance * bh_mass, theta * M_PI / 180, bh_mass, 0);
		IR->calculate();
		 //std::cout << "Calculating Ghost image" << std::endl;
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
		if (!loop) { break; };
	}
	escape = _getch();
	if (escape == 27) { std::cout << "exited: " << std::endl; };


	return 0;
}