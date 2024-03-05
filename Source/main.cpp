//Unofficial C++ version of Luminet project of bgmeulem
//https://github.com/bgmeulem/Luminet/

#include "Tests.h"




int main() {
	Tests test;
	unsigned test_choice =1;
	const unsigned type_of_plot = 3;
	const double inclination = 85.0;

	const double bh_mass =1.;// provide black hole mass value
	const double radius =70.0*bh_mass;//  provide radius value > 6M
	
	int order = 0;
	switch (test_choice) {
	case 1:
		/*
		************************************************************************************************************************************
		MAKE A CHOICE with type_of_plot:
								1= picture of a bare isoradial (no redshift) (make sure you choose the desired inclination)
								2= picture of an isoradial (with redshift) (make sure you choose the desired inclination)
								3= picture of an animated isoradial (with redshift) (standard: step=1° and between 0 and 180°)
		************************************************************************************************************************************
*/
		test.test_iso_radials(bh_mass, radius, inclination,type_of_plot);
		break;
	case 2:
		/*
		************************************************************************************************************************************
		MAKE A CHOICE with type_of_plot:
								1= picture of an isoredshift (make sure you choose the desired inclination)
								2= picture of an animated isoredshifts  (standard: step=1° and between 0 and 180°)
		************************************************************************************************************************************
*/
		test.test_iso_redshifts(bh_mass, radius, inclination,order,type_of_plot);
		break;
	case 3:
		test.test_accretion_disk(bh_mass, radius, inclination, order, type_of_plot);
		break;
	default:
		break;
	};
	return 0;
};

