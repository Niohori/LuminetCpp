#include "IsoRedShift.h"

IsoRedShift::IsoRedShift() {
	//constructor;
};

IsoRedShift::~IsoRedShift() {
	//destructor;
};
IsoRedShift::IsoRedShift(const double& angle, const double& bh_mass, const double& _lower_radius_, const double& _upper_radius_, const size_t& _n_radii_, const size_t& _n_angles_, const double& redshift_treschold_)
	:
	theta_0(angle),
	rS(2.0 * bh_mass),
	rIsco(_lower_radius_),
	M(bh_mass),
	lower_radius(_lower_radius_),
	upper_radius(_upper_radius_),
	n_radii(_n_radii_), n_angles(_n_angles_),
	redshift_treshold(redshift_treschold_) {
	make_grid();
}
IsoRedShift::IsoRedShift(const double& angle, const double& redshift_, const double& bh_mass, const std::map<double, std::pair<int, Isoradial*> >& from_isoradials)
	:
	theta_0(angle),
	M(bh_mass),
	rS(2.0 * bh_mass),
	rIsco(6.0 * bh_mass),
	redshift(redshift_) {
}


/**
	===================================================================================================================================
	* @brief Creates the grid where redshiftfactors will be calculated
	*
	* @param[in] None
	*
	*
	* @return None Updates the class members.
	*
	=====================================================================================================================================*/
void IsoRedShift::make_grid() {
	if (std::acos(std::cos(theta_0)) < M_PI / 4 || std::acos(std::cos(theta_0)) > 3 * M_PI / 4) {
		radii = OperatorsOrder2::linspace(lower_radius, upper_radius, n_radii);
		std::cout << "Linear raddi" << std::endl;
	}
	else {
		radii = OperatorsOrder2::logspace(lower_radius, upper_radius, n_radii);
		std::cout << "Log raddi" << std::endl;
	}

	//angles = OperatorsOrder2::linspace(0, 2 * M_PI, n_angles);
	angles.clear();
}

std::pair<std::vector<double>, std::vector<double>> IsoRedShift::get_ConcaveHull() {
	std::vector<double> x_hull;
	std::vector<double> y_hull;
	for (std::size_t i = 0; i < ConcaveHull.size(); i += 2)
	{
		const double x = ConcaveHull[i].x;
		const double y = ConcaveHull[i].y;
		x_hull.push_back(x);
		y_hull.push_back(y);
	}
	return std::make_pair(x_hull, y_hull);
}

/**
	===================================================================================================================================
	* @brief Core pubilc function to get the redshift isolines (manages the Mesh and Isolines classes) 
	*
	* @param[in] size_t: number of isolines wanted (equaly space between min and max redshift factor
	*
	*
	* @return std::multimap<double, std::vector<meshes::Point> > a map (entry is an isoline value) with  vectors of points representing the isolines
	*
	=====================================================================================================================================*/
std::multimap<double, std::vector<meshes::Point> > IsoRedShift::get_isolines(const size_t number_of_isolines) {//}, const std::vector<double>& xIsco, const std::vector<double>& yIsco) {
	std::multimap<double, std::vector<meshes::Point> > IsoLines;
	double impact_parameter = 0.0;
	double redshift;
	x_max = -100000000000000.0;
	x_min = -x_max;
	y_max = x_max;
	y_min = x_min;
	redshift_max = -1000000000000.0;
	redshift_min = -redshift_max;
	//we first calculate the redshift for de polar grid, convert the polar to cartesian and store the triple(x,y,redshift) in a vector of struct "Point".
	//Isoradial  IR;
	{//scoping bracket (Isoradial object destroyed after usage...
		Isoradial IR(radii[0] * M, theta_0, M, 0);//calculate the ISCO for one of the radii (which one is of no importance
		IR.calculate();
		IR.calculate_ISCO_region();
		auto ISCO = IR.get_ISCO_curve();
		xIsco = ISCO.first;
		yIsco = ISCO.second;
	}//end of scoping bracket
	for (auto rad : radii) {
		Isoradial IR(rad * M, theta_0, M, 0);
		IR.calculate();
		double sizingFactor = 1.0;
		auto bare_isoradials = IR.get_bare_isoradials();
		std::vector<double> xx = std::get<0>(bare_isoradials);
		std::vector<double> yy = std::get<1>(bare_isoradials);
		auto redshift_factors = IR.get_redshift_factors();
		for (int i = 0; i < xx.size(); i++) {
			xx[i] /= sizingFactor;
			yy[i] /= sizingFactor;
			if (xx[i] >= x_max) x_max = xx[i];
			if (xx[i] < x_min) x_min = xx[i];
			if (yy[i] >= y_max) y_max = yy[i];
			if (yy[i] < y_min) y_min = yy[i];
			if (redshift_factors[i] >= redshift_max)redshift_max = redshift_factors[i];
			if (redshift_factors[i] < redshift_min)redshift_min = redshift_factors[i];
			xCoordinates.push_back(xx[i]);//TEMPORARY
			yCoordinates.push_back(yy[i]);;//TEMPORARY
			redshifts.push_back((redshift_factors[i]) * 1.0);//TEMPORARY
		}
	}
	// take equal boundaries for the x an y axis
	x_max = std::max(std::max(x_max, std::abs(x_min)) * 1.1, std::max(y_max, std::abs(y_min)) * 1.1);
	x_min = -x_max;
	maxCoordinate = x_max;
	std::vector<double> rsLevels = OperatorsOrder2::linspace(redshift_min, redshift_max, number_of_isolines);
	std::shared_ptr<meshes::Mesh> mesh = std::make_shared<meshes::Mesh>(xCoordinates, yCoordinates, redshifts, xIsco, yIsco);
	ConcaveHull = mesh->concaveHull;
	ISCO = mesh->ISCO;
	parallel_for(0, rsLevels.size(), [&](int i) {
		//for (int i = 0; i < rsLevels.size(); i++) {
		double iso = rsLevels[i];
		Isolines isoLines(mesh, iso);
		std::vector<meshes::Point> anIsoLine = isoLines.get_iso_lines();
		IsoLines.insert(std::make_pair(iso, anIsoLine));
		});
	return  IsoLines;
}

double IsoRedShift::calculateDistance(const meshes::Point& p1, const meshes::Point& p2) {
	// Function to calculate the Euclidean distance between two points
	return std::sqrt(std::pow(p2.x - p1.x, 2) + std::pow(p2.y - p1.y, 2));
}

double IsoRedShift::findSmallestDistance(const std::vector<meshes::Point>& dataset) {
	// Function to find the smallest distance between all points in the dataset
	size_t size = dataset.size();
	double minDistance = std::numeric_limits<double>::infinity();
	// Iterate through all pairs of points
	for (size_t i = 0; i < size - 1; ++i) {
		for (size_t j = i + 1; j < size; ++j) {
			double distance = calculateDistance(dataset[i], dataset[j]);
			// Update minimum distance if a smaller distance is found
			if (distance < minDistance) {
				minDistance = distance;
			}
		}
	}
	return minDistance;
}

std::pair<std::vector<double>, std::vector<double>> IsoRedShift::get_ISCO_curve() {
	std::vector<double> xIsco;
	std::vector<double> yIsco;
	for (std::size_t i = 0; i < ISCO.size(); i++)
	{
		xIsco.push_back(ISCO[i].x);
		yIsco.push_back(ISCO[i].y);
	};
	ISCO_boundary = std::make_pair(xIsco, yIsco);
	return ISCO_boundary;
}

double  IsoRedShift::getMaxX() {
	return std::abs(maxCoordinate);
}