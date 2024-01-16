//Unoffial C++ version of Luminet project of bgmeulem
//https://github.com/bgmeulem/Luminet/


#include "BlackHole.h"
#include "BlackHolePhysics.h"


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

	double result_phi_inf = BHp.phi_inf(periastron, bh_mass);

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

	
	//BH.calc_isoredshifts();
	return 0;
}