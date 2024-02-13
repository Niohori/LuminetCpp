//Unofficial C++ version of Luminet project of bgmeulem
//https://github.com/bgmeulem/Luminet/

#include "Tests.h"




int main() {
	Tests test;
	unsigned test_choice =3;
	const unsigned type_of_plot =2;
	const double inclination = 89.0;

	const double bh_mass = 1.;// provide black hole mass value
	const double radius = 30.0;// 20. provide radius value
	

	int order = 0;

	switch (test_choice) {
	case 1:
		test.test_functions(bh_mass, radius);
		break;
	case 2:
		test.test_BH_rendering(bh_mass, radius);
		break;

	case 3:
		test.test_iso_radials(bh_mass, radius, inclination,type_of_plot);
		break;
	case 4:
		test.test_iso_redshifts(bh_mass, radius, inclination,order,type_of_plot);
		break;
	default:
		// code to be executed if
		// expression doesn't match any constant
		break;
	};

	return 0;
};

