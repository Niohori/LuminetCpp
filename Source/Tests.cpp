#include "Tests.h"

/*
************************************************************************************************************************************
TEMPORARY MAIN: USED FOR TESTINGS.
************************************************************************************************************************************/

void print_ISCO(const std::pair< std::vector<double>, const std::vector<double> > Isco_curve, const double& inclination) {
	// Convert double to string with 2-digit precision
	std::string inclination_as_string = std::to_string(inclination);
	size_t dotPos = inclination_as_string.find('.');
	if (dotPos != std::string::npos && dotPos + 0 < inclination_as_string.size()) {
		inclination_as_string = inclination_as_string.substr(0, dotPos + 0); // keep 2 digits after the dot
	}

	std::string filename = "data_dump_of_ISCO_curve_" + inclination_as_string + ".txt";

	// Open the file for writing
	std::ofstream outputFile(filename);

	// Check if the file is opened successfully
	if (outputFile.is_open()) {
		outputFile << "x;y;";
		for (int i = 0; i < Isco_curve.first.size(); i++) {
			outputFile << std::endl << Isco_curve.first[i] << ";" << Isco_curve.second[i] << ";";
		}
		// Close the file
		outputFile.close();

		std::cout << "ISCO-curve has been written to " << filename << std::endl;
	}
	else {
		std::cerr << "Error opening file: " << filename << std::endl;
	}

	std::cout << "end now " << std::endl;
}

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
	double distance = radius;
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
	std::vector<std::unique_ptr<Isoradial> > IR;;
	std::vector<std::unique_ptr<Isoradial> > IRg;
	//end of setting
	/*
			************************************************************************************************************************************
			MAKE A CHOICE with my_switch:
									1= picture of a bare isoradial (no redshift) (make sure you choose the desired inclination)
									2= picture of an isoradial (with redshift) (make sure you choose the desired inclination)
									3= picture of an animated isoradial (with redshift) (standard: step=1° and between 0 and 180°)
			************************************************************************************************************************************
	*/
	//my_switch = 2;
	//inclination = 89.1;//in degrees @90° gives a strange result: to investigate the maths
	switch (my_switch) {
	case 1:

		loop = false;
		IR.push_back(std::make_unique<Isoradial>(distance * bh_mass, inclination * M_PI / 180, bh_mass, 0));
		IRg.push_back(std::make_unique<Isoradial>(distance * bh_mass, inclination * M_PI / 180, bh_mass, 1));
		theta = inclination;
		//IR = new Isoradial(distance * bh_mass, theta * M_PI / 180, bh_mass, 0);
		IR[0]->calculate(false);
		//std::cout << "Calculating Ghost image" << std::endl;
		//IRg = new Isoradial(distance * bh_mass, theta * M_PI / 180, bh_mass, 1);//ghost
		IRg[0]->calculate(false);
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
		IR.push_back(std::make_unique<Isoradial>(distance * bh_mass, inclination * M_PI / 180, bh_mass, 0));
		IRg.push_back(std::make_unique<Isoradial>(distance * bh_mass, inclination * M_PI / 180, bh_mass, 1));
		theta = inclination;
		IR[0]->calculate(false);
		IRg[0]->calculate(false);
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
			IR.push_back(std::make_unique<Isoradial>(distance * bh_mass, p * M_PI / 180, bh_mass, 0));
			IRg.push_back(std::make_unique<Isoradial>(distance * bh_mass, p * M_PI / 180, bh_mass, 1));
		}
		while (true)
		{
			theta = inclinations[count];
			IR[count]->calculate(false);
			IRg[count]->calculate(false);
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
TESTS OF PLOTTING ISOREDSHIFTS
************************************************************************************************************************************/

/************************************************************************************************************************************
Version 1: all isoradials with rs
************************************************************************************************************************************/

void Tests::test_iso_redshifts(const double& bh_mass, const double& radius, const double& inclination, const int& order, const unsigned& my_switch) {
	const double lower_radius = 6.0;
	const double upper_radius = 50.0;
	const size_t _n_radii_ = 20;
	const size_t _n_angles_ = size_t(360. / 2.);
	double redshift_treschold_ = 0.1;
	bool loop = false;
	std::vector<double> X;
	std::vector<double> Y;
	char escape;
	//std::cout << "press esc to exit! " << std::endl;
	//bool loop = true;
	std::vector<double> inclinations;//will contain the inclinations for the animated vis (in degrees)

	//double inclination = 80.0;//in degrees
	double distance = radius;

	int sign = 1;//used for the loop (animation)
	int count = 0;//used for the loop (animation)
	unsigned index = 0;//used for the loop (animation)
	double theta = 0.0;//used for the loop (animation)
	std::vector<std::unique_ptr<IsoRedShift> > Irs;;

	//std::cout << "Inclination : " << inclination << std::endl;
	Isoradial Ir(0, inclination * M_PI / 180, bh_mass, 0);
	Ir.calculate_ISCO_region();
	auto Isco_curve = Ir.get_ISCO_curve();

	//print_ISCO(Isco_curve, inclination );

	//IsoRedShift Irs(inclination * M_PI / 180, bh_mass, lower_radius, upper_radius, _n_radii_, _n_angles_, redshift_treschold_);
	////std::multimap<double, std::vector<delaunay::Segment> > isolines = Irs.get_isolines_2(15, Isco_curve.first, Isco_curve.second);
	double xmax = 0.0;
	double ymax = 0.0;
	double xmin = 0.0;
	double ymin = 0.0;
	double rsmax = 0.0;
	double rsmin = 0.0;
	//Plotter plot;
	//plot.plot_iso_redshifts(inclination, isolines, xmax, xmin, rsmax, rsmin, false);
	std::multimap<double, std::vector<delaunay::Segment> > isolines;
	std::vector< std::multimap<double, std::vector<delaunay::Segment> > > allIsolines;
	Plotter plot;
	//==================================================================
	//==================================================================
	//==================================================================

	double STEP = 1.0;//the resolution (in degrees )of animated views
	double nIsolines =10;//numebr of isolines wanted

	//==================================================================
	//==================================================================
	//==================================================================
	std::vector<double> xCoordinates;
	std::vector<double> yCoordinates;
	std::vector<double> redshifts;
	double  t0 = 0.0;// std::chrono::high_resolution_clock::now();
	auto t1 = std::chrono::high_resolution_clock::now();
	auto t2 = std::chrono::high_resolution_clock::now();

	switch (my_switch) {
	case 1:

		loop = false;
		t1 = std::chrono::high_resolution_clock::now();
		Irs.push_back(std::make_unique<IsoRedShift>(inclination * M_PI / 180, bh_mass, lower_radius, upper_radius, _n_radii_, _n_angles_, redshift_treschold_));
		t2 = std::chrono::high_resolution_clock::now();
		std::cout << "Average time (s) for the redshift objectts creation: " << (t2 - t1).count() / 1e9 << std::endl;
		std::cout << "==================================================================" << std::endl << std::endl;
		theta = inclination;
		//IR = new Isoradial(distance * bh_mass, theta * M_PI / 180, bh_mass, 0);
		isolines = Irs[0]->get_isolines(nIsolines, Isco_curve.first, Isco_curve.second);
		xmax = Irs[0]->x_max;
		ymax = Irs[0]->y_max;
		xmin = Irs[0]->x_min;
		ymin = Irs[0]->y_min;
		rsmax = Irs[0]->redshift_max;
		rsmin = Irs[0]->redshift_min;
		t1 = std::chrono::high_resolution_clock::now();
		plot.plot_iso_redshifts(inclination, isolines, xmax, xmin, rsmax, rsmin, false);
		t2 = std::chrono::high_resolution_clock::now();
		std::cout << "Average time (s) for the isolines calculation: " << (t2 - t1).count() / 1e9 << std::endl;
		count = 0;
		break;
	case 2:

		loop = true;

		for (double p = 1; p < 180.0; p += STEP) {
			t1 = std::chrono::high_resolution_clock::now();
			if (p == 90.0)continue;
			inclinations.push_back(p);
			Irs.push_back(std::make_unique<IsoRedShift>(p * M_PI / 180, bh_mass, lower_radius, upper_radius, _n_radii_, _n_angles_, redshift_treschold_));
			t2 = std::chrono::high_resolution_clock::now();
			t0 += (t2 - t1).count() / 1e9;
			std::cout << "Generated redshift object for "<< p << "°" <<std::endl;
		}
		std::cout << "Average time (s) for the redshift objectts creation: " << t0 / inclinations.size() << std::endl;
		std::cout << "==================================================================" << std::endl << std::endl;
		t0 = 0;
		for (int i = 0; i < inclinations.size(); i++) {
			t1 = std::chrono::high_resolution_clock::now();
			Isoradial Ir(0, inclination * M_PI / 180, bh_mass, 0);
			Ir.calculate_ISCO_region();
			auto Isco_curve = Ir.get_ISCO_curve();
			isolines = Irs[i]->get_isolines(nIsolines, Isco_curve.first, Isco_curve.second);
			allIsolines.push_back(isolines);
			t2 = std::chrono::high_resolution_clock::now();
			t0 += (t2 - t1).count() / 1e9;
			std::cout << "Calculated redshift isolines  for " << inclinations[i] << "°" << std::endl;
		}
		std::cout << "Average time (s) for the isolines calculation: " << t0 / inclinations.size() << std::endl;
		count = 0;
		while (true)
		{
			//Isoradial Ir(0, inclination * M_PI / 180, bh_mass, 0);
			//Ir.calculate_ISCO_region();
			//auto Isco_curve = Ir.get_ISCO_curve();

			//print_ISCO(Isco_curve, inclinations[count]);

			theta = inclinations[count];
			auto anIsolines = allIsolines[count];
			xmax = Irs[count]->x_max;
			ymax = Irs[count]->y_max;
			xmin = Irs[count]->x_min;
			ymin = Irs[count]->y_min;
			rsmax = Irs[count]->redshift_max;
			rsmin = Irs[count]->redshift_min;

			plot.plot_iso_redshifts(inclinations[count], anIsolines, xmax, xmin, rsmax, rsmin, true);

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
	case 3:

		loop = false;
		Irs.push_back(std::make_unique<IsoRedShift>(inclination * M_PI / 180, bh_mass, lower_radius, upper_radius, _n_radii_, _n_angles_, redshift_treschold_));

		theta = inclination;
		//IR = new Isoradial(distance * bh_mass, theta * M_PI / 180, bh_mass, 0);
		isolines = Irs[0]->get_isolines(nIsolines, Isco_curve.first, Isco_curve.second);
		xmax = Irs[0]->x_max;
		ymax = Irs[0]->y_max;
		xmin = Irs[0]->x_min;
		ymin = Irs[0]->y_min;
		rsmax = Irs[0]->redshift_max;
		rsmin = Irs[0]->redshift_min;
		xCoordinates = Irs[0]->xCoordinates;
		yCoordinates = Irs[0]->yCoordinates;
		redshifts = Irs[0]->redshifts;
		plot.plot_iso_redshifts(inclination, xCoordinates, yCoordinates, redshifts, xmax, xmin, rsmax, rsmin);
		break;

	default:
		// code to be executed if
		// expression doesn't match any constant
		break;
	};
};