#include "Tests.h"

/*
************************************************************************************************************************************
TEMPORARY MAIN: USED FOR TESTINGS.
************************************************************************************************************************************/

void print_eq13(const std::vector<double> perias, const std::vector<double> vals, const double& inclination, const double& phi, const double& radius, const double& alfa) {
	// Convert double to string with 2-digit precision
	std::string inclination_as_string = std::to_string(inclination);
	size_t dotPos = inclination_as_string.find('.');
	if (dotPos != std::string::npos && dotPos + 0 < inclination_as_string.size()) {
		inclination_as_string = inclination_as_string.substr(0, dotPos + 0); // keep 2 digits after the dot
	}

	std::string filename = "data_dump_of_Eq13_" + inclination_as_string + ".txt";

	// Open the file for writing
	std::ofstream outputFile(filename);

	// Check if the file is opened successfully
	if (outputFile.is_open()) {
		outputFile << "Radius=" << radius << ";Phi=" << phi << ";Alpha=" << alfa << ";Inclination=" << inclination << "; ";
		for (int i = 0; i < vals.size(); i++) {
			outputFile << std::endl << perias[i] << ";" << vals[i] << ";";
		}
		// Close the file
		outputFile.close();

		std::cout << "Equation has been written to " << filename << std::endl;
	}
	else {
		std::cerr << "Error opening file: " << filename << std::endl;
	}

	std::cout << "end now " << std::endl;
}

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
	std::pair<std::vector<double>, std::vector<double>> bare_ghost_isoradials;//temporary functions for debugging
	std::vector< std::pair<std::vector<double>, std::vector<double> > > bareIsoradialsContainer;
	std::vector< std::pair<std::vector<double>, std::vector<double> > > bareGhostIsoradialsContainer;
	std::vector<double> xx;
	std::vector<double> yy;

	std::vector<double> xx_ghost;
	std::vector<double> yy_ghost;
	std::vector<double> redshift_factors;
	std::vector<double> redshift_factors_ghost;
	std::vector< std::vector<double> >  redshiftFactorsContainer;
	std::vector< std::vector<double> >  redshiftFactorsGhostContainer;
	int sign = 1;//used for the loop (animation)
	int count = 0;//used for the loop (animation)
	unsigned index = 0;//used for the loop (animation)
	double theta = 0.0;//used for the loop (animation)
	std::vector<std::unique_ptr<Isoradial> > IR;;
	std::vector<std::unique_ptr<Isoradial> > IRg;
	double Xmax;
	std::vector<double> XmaxX;
	auto t1 = std::chrono::high_resolution_clock::now();
	auto t2 = std::chrono::high_resolution_clock::now();
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
		Xmax = IR[0]->getMaxX();
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
		Xmax = IR[0]->getMaxX();
		plotter.plot_isoradials(theta, xx, yy, xx_ghost, yy_ghost, redshift_factors, redshift_factors_ghost, Xmax, loop);
		break;

	case 3:
		loop = true;
		std::cout << std::endl << "================================== Calculating isoradials ...  ===========================================" << std::endl << std::endl;
		t1 = std::chrono::high_resolution_clock::now();
		for (double p = 1; p < 180.0; p += step) {
			inclinations.push_back(p);
		}
		inclinations.clear();
		for (double p = 0.0; p < 180.0; p += 1.0) {
			double an = 90.0 * std::sin(p/180.0*M_PI );
			if (p > 90) {
				an = 180.0 - 90.0 * std::sin(p / 180.0 * M_PI);
			}
			inclinations.push_back(an);
		}

		IR = std::vector<std::unique_ptr<Isoradial> >(inclinations.size());
		IRg = std::vector<std::unique_ptr<Isoradial> >(inclinations.size());
		XmaxX = std::vector<double>(inclinations.size());
		bareIsoradialsContainer = std::vector< std::pair<std::vector<double>, std::vector<double> > >(inclinations.size());
		bareGhostIsoradialsContainer = std::vector< std::pair<std::vector<double>, std::vector<double> > >(inclinations.size());
		redshiftFactorsContainer = std::vector< std::vector<double> >(inclinations.size());
		redshiftFactorsGhostContainer = std::vector< std::vector<double> >(inclinations.size());
		parallel_for(0, inclinations.size(), [&](int i) {
			//for (int i = 0; i < inclinations.size(); i++) {
			IR[i] = std::make_unique<Isoradial>(distance * bh_mass, inclinations[i] * M_PI / 180, bh_mass, 0);
			IRg[i] = std::make_unique<Isoradial>(distance * bh_mass, inclinations[i] * M_PI / 180, bh_mass, 1);
			theta = inclinations[i];
			IR[i]->calculate(false);
			IRg[i]->calculate(false);
			bare_isoradials = IR[i]->get_bare_isoradials();//temporary functions for debugging
			bare_ghost_isoradials = IRg[i]->get_bare_isoradials();
			bareIsoradialsContainer[i] = bare_isoradials;
			bareGhostIsoradialsContainer[i] = bare_ghost_isoradials;
			redshift_factors = IR[i]->get_redshift_factors();
			redshift_factors_ghost = IRg[i]->get_redshift_factors();
			redshiftFactorsContainer[i] = redshift_factors;
			redshiftFactorsGhostContainer[i] = redshift_factors_ghost;
			XmaxX[i] = IR[i]->getMaxX();
			}
		);
		Xmax = *std::max_element(XmaxX.begin(), XmaxX.end());
		t2 = std::chrono::high_resolution_clock::now();
		std::cout << "Time (s) for " << inclinations.size() << " isoradials: " << (t2 - t1).count() / 1e9 << " s." << " (" << (t2 - t1).count() / 1e6 / inclinations.size() << " ms/isoradial)." << std::endl;
		std::cout << std::endl << "=====================================================================================================================================================" << std::endl << std::endl;

		while (true)
		{
			theta = inclinations[count];
			//IR[count]->calculate(false);
			//IRg[count]->calculate(false);
			bare_isoradials = bareIsoradialsContainer[count];
			xx = std::get<0>(bare_isoradials);
			yy = std::get<1>(bare_isoradials);
			bare_ghost_isoradials = bareGhostIsoradialsContainer[count];
			xx_ghost = std::get<0>(bare_ghost_isoradials);
			yy_ghost = std::get<1>(bare_ghost_isoradials);
			redshift_factors = redshiftFactorsContainer[count];
			redshift_factors_ghost = redshiftFactorsGhostContainer[count];
			plotter.plot_isoradials(theta, xx, yy, xx_ghost, yy_ghost, redshift_factors, redshift_factors_ghost, Xmax, loop);
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
			//std::this_thread::sleep_for(std::chrono::milliseconds(250));
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

void Tests::test_iso_redshifts(const double& bh_mass, const double& radius, const double& inclination, const int& order, const unsigned& my_switch) {
	//==================================================================
//==================================================================
//==================================================================
	double STEP = 1.0;//the resolution (in degrees )of animated views
	double nIsolines = 10;//number of isolines wanted
	const size_t _n_radii_ = 30;
	//==================================================================
	//==================================================================
	//==================================================================
	const double lower_radius = 6 * bh_mass;
	const double upper_radius = 50.0 * bh_mass;

	const size_t _n_angles_ = size_t(360. / 1.);
	double redshift_treschold_ = 0.1;
	bool loop = false;
	std::vector<double> X;
	std::vector<double> Y;
	std::vector<double> inclinations;//will contain the inclinations for the animated vis (in degrees)

	double distance = radius;

	int sign = 1;//used for the loop (animation)
	int count = 0;//used for the loop (animation)
	unsigned index = 0;//used for the loop (animation)
	double theta = 0.0;//used for the loop (animation)
	std::vector<std::unique_ptr<IsoRedShift> > Irs;;
	double xmax = 0.0;
	double ymax = 0.0;
	double xmin = 0.0;
	double ymin = 0.0;
	double rsmax = 0.0;
	double rsmin = 0.0;
	std::multimap<double, std::vector<meshes::Point> > isolines;
	std::vector< std::multimap<double, std::vector<meshes::Point> > > allIsolines;
	Plotter plot;

	std::vector<double> xCoordinates;
	std::vector<double> yCoordinates;
	std::vector<double> redshifts;
	auto t1 = std::chrono::high_resolution_clock::now();
	auto t2 = std::chrono::high_resolution_clock::now();

	std::vector<std::pair<std::vector<double>, std::vector<double>> > Iscos;
	std::vector< std::pair<std::vector<double>, std::vector<double>> > concaveHulls;
	double Xmax;
	std::vector<double> XmaxX;
	switch (my_switch) {
	case 1:

		loop = false;
		t1 = std::chrono::high_resolution_clock::now();

		Irs.push_back(std::make_unique<IsoRedShift>(inclination * M_PI / 180, bh_mass, lower_radius, upper_radius, _n_radii_, _n_angles_, redshift_treschold_));

		t2 = std::chrono::high_resolution_clock::now();
		std::cout << "Average time (s) for the redshift objects creation: " << (t2 - t1).count() / 1e6 << " ms." << std::endl;
		std::cout << "==================================================================" << std::endl << std::endl;
		theta = inclination;
		//IR = new Isoradial(distance * bh_mass, theta * M_PI / 180, bh_mass, 0);
		t1 = std::chrono::high_resolution_clock::now();
		isolines = Irs[0]->get_isolines(nIsolines);// , Isco_curve.first, Isco_curve.second);
		concaveHulls.push_back(Irs[0]->get_ConcaveHull());
		Iscos.push_back(Irs[0]->get_ISCO_curve());
		Xmax = Irs[0]->getMaxX();
		t2 = std::chrono::high_resolution_clock::now();
		std::cout << "Time (s) for the isolines calculation: " << (t2 - t1).count() / 1e6 << " ms." << std::endl;
		rsmax = Irs[0]->redshift_max;
		rsmin = Irs[0]->redshift_min;
		t1 = std::chrono::high_resolution_clock::now();
		plot.plot_iso_redshifts(bh_mass, inclination, isolines, Xmax, -Xmax, rsmax, rsmin, Iscos[0], concaveHulls[0], false);
		t2 = std::chrono::high_resolution_clock::now();
		std::cout << "Time (s) for the plotting: " << (t2 - t1).count() / 1e6 << " ms." << std::endl;
		count = 0;
		break;
	case 2:

		loop = true;
		t1 = std::chrono::high_resolution_clock::now();
		count = 0;
		for (double p = 1; p < 91.0; p += STEP) {
			if (p == 90.0)continue;
			inclinations.push_back(p);
			//Irs.push_back(std::make_unique<IsoRedShift>(p * M_PI / 180, bh_mass, lower_radius, upper_radius, _n_radii_, _n_angles_, redshift_treschold_));;
		}
		Irs = std::vector<std::unique_ptr<IsoRedShift> >(inclinations.size());
		parallel_for(0, inclinations.size(), [&](int i) {
			Irs[i] = std::make_unique<IsoRedShift>(inclinations[i] * M_PI / 180, bh_mass, lower_radius, upper_radius, _n_radii_, _n_angles_, redshift_treschold_);
			}
		);
		t2 = std::chrono::high_resolution_clock::now();
		std::cout << "Time (s) for redshift objects creation for " << inclinations.size() << " inclinations: " << (t2 - t1).count() / 1e9 << " s." << " (" << (t2 - t1).count() / 1e6 / inclinations.size() << " ms/inclination)." << std::endl;
		std::cout << std::endl << "================================== Calculating isolines ...  ===========================================" << std::endl << std::endl;
		t1 = std::chrono::high_resolution_clock::now();
		allIsolines = std::vector< std::multimap<double, std::vector<meshes::Point> > >(inclinations.size(), isolines);
		concaveHulls = std::vector<std::pair<std::vector<double>, std::vector<double>> >(inclinations.size());
		Iscos = std::vector<std::pair<std::vector<double>, std::vector<double>> >(inclinations.size());
		XmaxX = std::vector<double>(inclinations.size());
		parallel_for(0, inclinations.size(), [&](long i) {
			//for (int i = 0; i < inclinations.size(); i++) {
			std::multimap<double, std::vector<meshes::Point> > iso = Irs[i]->get_isolines(nIsolines);;
			concaveHulls[i] = Irs[i]->get_ConcaveHull();
			Iscos[i] = Irs[i]->get_ISCO_curve();
			Xmax = Irs[i]->getMaxX();
			XmaxX[i] = Xmax;;
			allIsolines[i] = iso;
			}
		);
		t2 = std::chrono::high_resolution_clock::now();
		std::cout << "Time (s) for " << nIsolines << " isolines calculation for " << inclinations.size() << " inclinations: " << (t2 - t1).count() / 1e9 << " s." << " (" << (t2 - t1).count() / 1e6 / inclinations.size() / nIsolines << " ms/isoline/inclination)." << std::endl;
		std::cout << "=========================================================================================================" << std::endl << std::endl;
		Xmax = *std::max_element(XmaxX.begin(), XmaxX.end());
		count = 0;
		while (true)
		{
			theta = inclinations[count];
			std::multimap<double, std::vector<meshes::Point> >  anIsolines = allIsolines[count];
			rsmax = Irs[count]->redshift_max;
			rsmin = Irs[count]->redshift_min;

			plot.plot_iso_redshifts(bh_mass, inclinations[count], anIsolines, Xmax, -Xmax, rsmax, rsmin, Iscos[count], concaveHulls[count], true);

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
			std::this_thread::sleep_for(std::chrono::milliseconds(100));
		}
		break;

	default:
		// code to be executed if
		// expression doesn't match any constant
		break;
	};
};

void Tests::test_accretion_disk(const double& bh_mass, const double& radius, const double& inclination, const int& order, const unsigned& my_switch) {
	accretiondisk::AccretionDisk aDisk(bh_mass, inclination * M_PI / 180, 6 * bh_mass, radius, 1000);;
	double maxFluxP = aDisk.getMaxFluxPrimary();
	double maxFluxS = aDisk.getMaxFluxSecundary();
	double minFluxP = aDisk.getMinFluxPrimary();
	double minFluxS = aDisk.getMinFluxSecundary();
	double power_scale = .9;
	std::cout << "Accretion Disk: calculations done " << std::endl;
	Plotter plot;
	while (true)
	{
		aDisk.updateDisk( 1.0);
	std::vector<accretiondisk::Particle> primaryParticles = aDisk.getPrimaryImage();
	std::vector<accretiondisk::Particle> secundaryParticles = aDisk.getSecundaryImage();

	std::vector<double> fluxesPrimary;
	std::vector<double> xPrimary;
	std::vector<double> yPrimary;
	std::vector<double> fluxesSecundary;
	std::vector<double> xSecundary;
	std::vector<double> ySecundary;
	for (const auto& point : primaryParticles) {
		double fl = std::pow((std::abs(point.Fo) + minFluxP) / (maxFluxP + minFluxP), power_scale);
		fluxesPrimary.push_back(std::pow((std::abs(point.Fo) + minFluxP) / (maxFluxP + minFluxP), power_scale));//Enhances "contrast"
		xPrimary.push_back(point.x);
		yPrimary.push_back(point.y);
		//std::cout << "(x,y): (" << point.x << ", " << point.y << ") with flux :" << fl << std::endl;
	}
	for (const auto& point : secundaryParticles) {
		fluxesSecundary.push_back(std::pow((std::abs(point.Fo) + minFluxS) / (maxFluxS + minFluxS), power_scale));//Enhances "contrast"
		xSecundary.push_back(point.x);
		ySecundary.push_back(point.y);
	}
	plot.plot_BlackHole(inclination, xPrimary, yPrimary, fluxesPrimary, xSecundary, ySecundary, fluxesSecundary, true);
}
};

void Tests::test_display_eq13(const double& bh_mass, const double& radius, const double& incl, const int& order, const unsigned& my_switch) {
	std::vector<std::vector<double>> p;
	std::vector<std::vector<double>> v;
	double inclination = incl * M_PI / 180;
	double phi = 5 * M_PI / 180;
	double alpha = BHphysics::alpha(phi, inclination);
	double lower_bound = 2.0 * bh_mass;
	double upper_bound = 2.0 * radius;
	
	std::vector<double> periastrons = OperatorsOrder2::linspace(lower_bound, upper_bound, 1000);
	std::vector<double> valEq13;;
	for (auto p : periastrons) {
		double val = BHphysics::eq13(p, radius, alpha, bh_mass, inclination, order);
		valEq13.push_back(val);
		//if (val * val < 1e-4) { std::cout << "zero encountered" << std::endl; }
	}

	std::vector<double> valsEq13;
	for (auto p : periastrons) {
		double val = BHphysics::eq13(p, radius, alpha, bh_mass, inclination, order);
		valsEq13.push_back(val);
		//if (val * val < 1e-4) { std::cout << "zero encountered" << std::endl; }
	}
	print_eq13(periastrons, valsEq13, incl, phi * 180 / M_PI, radius, alpha * 180 / M_PI);
	Plotter plot;
	plot.plot_Eq13(periastrons, valEq13);
};