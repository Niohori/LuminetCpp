
#include "Tests.h"


	/*
************************************************************************************************************************************
TEMPORARY MAIN: USED FOR TESTINGS.
************************************************************************************************************************************



/*
************************************************************************************************************************************
SOME TESTS OF INDIVIDUAL FUNCTIONS
************************************************************************************************************************************
	*/
void Tests::test_functions(const double& bh_mass, const double& radius) {
	//INITIALIZATION OF SOME VARIABLES
	BHphysics BHp;
	double periastron = 30.;// provide periastron value ;
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
};

/************************************************************************************************************************************
TEST OF READINGG THE ACCRETION DISK DATA
************************************************************************************************************************************/
void Tests::test_BH_rendering(const double& bh_mass, const double& radius) {
	BlackHole BH(bh_mass);
	BH.sample_Sources(20000);// , "D:/dev/Cpp/Luminet/Points/points_incl=85b.csv", "D:/dev/Cpp/Luminet/Points/points_secondary_incl=85b.csv");
};




/************************************************************************************************************************************
TESTS OF PLOTTING SINGLE ISORADIALS
************************************************************************************************************************************/

void Tests::test_iso_radials(const double& bh_mass, const double& radius, const double& inclination, const unsigned& my_switch) {
	//setting up the plotting
	Plotter plotter;
	char escape;
	//std::cout << "press esc to exit! " << std::endl;
	bool loop = true;
	std::vector<double> inclinations;//will contain the inclinations for the animated vis (in degrees)
	double step = 1.0;//the resolution (in degrees )of animated views
	//double inclination = 80.0;//in degrees
	double distance = 30.0;
	std::pair<std::vector<double>, std::vector<double>> bare_isoradials;//temporary functions for debugging
	std::vector<double> xx;
	std::vector<double> yy;
	std::pair<std::vector<double>, std::vector<double>> bare_ghost_isoradials;//temporary functions for debugging
	std::vector<double> xx_ghost;
	std::vector<double> yy_ghost;
	std::vector<double> redshift_factors;
	std::vector<double> redshift_factors_ghost;
	int sign = 1;//used for the loop (animation)
	int count = 0;//used for the loop (animation)
	unsigned index = 0;//used for the loop (animation)
	double theta = 0.0;//used for the loop (animation)
	std::vector<Isoradial*> IR;;
	std::vector<Isoradial*> IRg;
	//end of setting
	/*
			************************************************************************************************************************************
			MAKE A CHOICE with my_switch:
									1= picture of a bare isoradial (no redshift) (make sure you choose the desired inclination)
									2= picture of an isoradial (with redshift) (make sure you choose the desired inclination)
									3= picture of an animated isoradial (with redshift) (standard: step=1� and between 0 and 180�)
			************************************************************************************************************************************
	*/
	//my_switch = 2;
	//inclination = 89.1;//in degrees @90� gives a strange result: to investigate the maths
	switch (my_switch) {
	case 1:
		loop = false;
		IR.push_back(new Isoradial(distance * bh_mass, inclination * M_PI / 180, bh_mass, 0));
		IRg.push_back(new Isoradial(distance * bh_mass, inclination * M_PI / 180, bh_mass, 1));
		theta = inclination;
		//IR = new Isoradial(distance * bh_mass, theta * M_PI / 180, bh_mass, 0);
		IR[0]->calculate();
		//std::cout << "Calculating Ghost image" << std::endl;
		//IRg = new Isoradial(distance * bh_mass, theta * M_PI / 180, bh_mass, 1);//ghost
		IRg[0]->calculate();
		bare_isoradials = IR[0]->get_bare_isoradials();//temporary functions for debugging
		xx = std::get<0>(bare_isoradials);
		yy = std::get<1>(bare_isoradials);
		bare_ghost_isoradials = IRg[0]->get_bare_isoradials();//temporary functions for debugging
		xx_ghost = std::get<0>(bare_ghost_isoradials);
		yy_ghost = std::get<1>(bare_ghost_isoradials);
		plotter.plot_isoradials(xx, yy, xx_ghost, yy_ghost);
		break;
	case 2:
		loop = false;
		IR.push_back(new Isoradial(distance * bh_mass, inclination * M_PI / 180, bh_mass, 0));
		IRg.push_back(new Isoradial(distance * bh_mass, inclination * M_PI / 180, bh_mass, 1));
		theta = inclination;
		IR[0]->calculate();
		IRg[0]->calculate();
		bare_isoradials = IR[0]->get_bare_isoradials();//temporary functions for debugging
		xx = std::get<0>(bare_isoradials);
		yy = std::get<1>(bare_isoradials);
		bare_ghost_isoradials = IRg[0]->get_bare_isoradials();//temporary functions for debugging
		xx_ghost = std::get<0>(bare_ghost_isoradials);
		yy_ghost = std::get<1>(bare_ghost_isoradials);
		redshift_factors = IR[0]->get_redshift_factors();
		redshift_factors_ghost = IRg[0]->get_redshift_factors();
		plotter.plot_isoradials(theta, xx, yy, xx_ghost, yy_ghost, redshift_factors, redshift_factors_ghost, loop);
		break;

	case 3:
		loop = true;
		for (double p = 1; p < 180.0; p += step) {
			inclinations.push_back(p);
			IR.push_back(new Isoradial(distance * bh_mass, p * M_PI / 180, bh_mass, 0));
			IRg.push_back(new Isoradial(distance * bh_mass, p * M_PI / 180, bh_mass, 1));
		}
		while (true)
		{
			theta = inclinations[count];
			IR[count]->calculate();
			IRg[count]->calculate();
			bare_isoradials = IR[count]->get_bare_isoradials();//temporary functions for debugging
			xx = std::get<0>(bare_isoradials);
			yy = std::get<1>(bare_isoradials);
			bare_ghost_isoradials = IRg[count]->get_bare_isoradials();//temporary functions for debugging
			xx_ghost = std::get<0>(bare_ghost_isoradials);
			yy_ghost = std::get<1>(bare_ghost_isoradials);
			redshift_factors = IR[count]->get_redshift_factors();
			redshift_factors_ghost = IRg[count]->get_redshift_factors();
			plotter.plot_isoradials(theta, xx, yy, xx_ghost, yy_ghost, redshift_factors, redshift_factors_ghost, loop);
			count = count + sign;
			if (count == -1) {
				sign = 1;
				count = 0;
			}
			if (count == inclinations.size()) {
				sign = -1;
				count = inclinations.size() - 1;
			}
			if (!loop) { break; };
			//escape = _getch();
			//if (escape == 27) { std::cout << "exited: " << std::endl; };
		}

		break;
	default:
		// code to be executed if
		// expression doesn't match any constant
		break;
	};
};



/************************************************************************************************************************************
TESTS OF PLOTTING SINGLE ISORADIALS
************************************************************************************************************************************/

void Tests::test_iso_redshifts(const double& bh_mass, const double& radius, const double& inclination, const int& order,const unsigned& my_switch) {
	std::map<double, std::pair<int, Isoradial*> > IR;
	double step = 10.0;//the resolution (in degrees )of animated views
	for (double radius =20;radius < 50.0; radius += step) {
		IR[radius]=std::make_pair(order, new Isoradial(radius* bh_mass, inclination* M_PI / 180, bh_mass, order));
	}
	for (auto ir:IR){
		ir.second.second->calculate();
	}
	double redshift = 0.25;
	IsoRedShift Irs(inclination, redshift, bh_mass,  IR);
	Irs.improve();
	std::pair<std::map<double, std::pair<std::vector<double>, std::vector<double>>>, std::map<double, std::pair<std::vector<double>, std::vector<double>>> > solutions = Irs.split_co_on_solutions();
	std::map<double, std::pair<std::vector<double>, std::vector<double>>> good = solutions.first;
	std::map<double, std::pair<std::vector<double>, std::vector<double>>> bad = solutions.second;
	for (auto asolution : good) {
		double rs = asolution.first;
		std::vector<double> x = asolution.second.first;
		std::vector<double> y = asolution.second.second;
		std::cout << "Redshift ? = " << rs << std::endl;
		for (size_t i = 0; i < x.size(); i++) {
			std::cout << "\t\t\t (x,y) = (" << x[i] << ", " << y[i] << ")" << std::endl;

		}
	}
	Plotter plot;
	plot.plot_redshifts(inclination,good);
};