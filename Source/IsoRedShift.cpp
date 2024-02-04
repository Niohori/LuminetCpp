#include "IsoRedShift.h"



void print_IsoLines(const std::vector<double>& X, const std::vector<double>& Y, const std::vector<double>& RS, const double& inclination) {
	// Convert double to string with 2-digit precision
	std::string inclination_as_string = std::to_string(inclination);
	size_t dotPos = inclination_as_string.find('.');
	if (dotPos != std::string::npos && dotPos + 0 < inclination_as_string.size()) {
		inclination_as_string = inclination_as_string.substr(0, dotPos + 0); // keep 2 digits after the dot
	}

	std::string filename = "data_dump_of_isoredshifts_" + inclination_as_string + ".txt";

	// Open the file for writing
	std::ofstream outputFile(filename);

	// Check if the file is opened successfully
	if (outputFile.is_open()) {
		outputFile << "x;y;rs;";
		for (int i = 0; i < X.size(); i++) {
			outputFile << std::endl << X[i] << ";" << Y[i] << ";" << RS[i] << ";" ;
		}
		// Close the file
		outputFile.close();

		std::cout << "IsoLines has been written to " << filename << std::endl;
	}
	else {
		std::cerr << "Error opening file: " << filename << std::endl;
	}

	std::cout << "end now " << std::endl;


}

IsoRedShift::IsoRedShift() {
	//constructor;
};

IsoRedShift::~IsoRedShift() {
	//destructor;
};
IsoRedShift::IsoRedShift(const double& angle, const double& bh_mass, const double& _lower_radius_, const double& _upper_radius_, const size_t& _n_radii_, const size_t& _n_angles_, const double& redshift_treschold_)
	:
	theta_0(angle),
	rS(2.0 * M),
	rIsco(3.0 * rS),
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
	rS(2.0 * M),
	rIsco(3.0 * rS),
	redshift(redshift_) {
}

void IsoRedShift::make_grid() {
	radii = OperatorsOrder2::linspace(lower_radius, upper_radius, n_radii);
	//angles = OperatorsOrder2::linspace(0, 2 * M_PI, n_angles);
	angles.clear();
}

std::multimap<double, std::vector<delaunay::Segment> > IsoRedShift::get_isolines(const size_t number_of_isolines,const std::vector<double>& xIsco, const std::vector<double>& yIsco) {
	std::multimap<double, std::vector<delaunay::Segment> > IsoLines;
	
	double impact_parameter = 0.0;

	//std::vector<delaunay::Point> all_points;
	double redshift;
	x_max = -100000000000000.0;
	x_min = -x_max;
	y_max = x_max;
	y_min = x_min;
	redshift_max = -1000000000000.0;
	redshift_min = -redshift_max;
	//we first calculate the redshift for de polar grid, convert the polar to cartesian and store the triple(x,y,redshift) in a vector of struct "cloud_points".
	Isoradial  IR;
	for (auto rad : radii) {
		Isoradial IR(rad * M, theta_0, M, 0);

		IR.calculate();
		double sizingFactor = 1.0;
		auto bare_isoradials = IR.get_bare_isoradials();//temporary functions for debugging
		std::vector<double> xx = std::get<0>(bare_isoradials);
		std::vector<double> yy = std::get<1>(bare_isoradials);
		auto redshift_factors = IR.get_redshift_factors();
		for (int i = 0; i < xx.size(); i++) {
			//if (std::abs(redshift_factors[i]) < 0.1) { continue; };
			//delaunay::Point apoint(xx[i], yy[i], redshift_factors[i]);
			//all_points.push_back(apoint);
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
			redshifts.push_back((redshift_factors[i])*1.0);//TEMPORARY
		}
	}
	//print_IsoLines(xCoordinates, yCoordinates, redshifts, theta_0* 180/M_PI );
	//std::cout << "Size x,y,rs : " << all_points.size() << std::endl;;

	// take equal boundaries for the x an y axis
	x_max = std::max(x_max, y_max);
	x_min = std::min(x_min, y_min);
	x_max = std::max(x_max, std::abs(x_min)) * 1.1;
	x_min = -x_max;// std::max(abs(x_max), abs(x_min)) * 1.1;
	//std::cout << "xmax,xmin= " << x_max << ", " << x_min << std::endl;
	std::vector<double> rsLevels = OperatorsOrder2::linspace(redshift_min, redshift_max, number_of_isolines);
	/*std::cout << "-------------------- " << rsLevels.size() << "  levels requested. " << std::endl;
	for (auto a : rsLevels) {

		std::cout << a << std::endl;
	}*/
	
	delaunay::Delaunay delaun(xCoordinates, yCoordinates, redshifts);
	/*std::vector<double> hull = delaun.concavehull(0.1);
	std::vector<double> isoline_values = OperatorsOrder2::linspace(redshift_min, redshift_max, 20);

	std::vector<double> x_hull;
	std::vector<double> y_hull;

	for (std::size_t i = 0; i < hull.size(); i += 2)
	{
		const double x = hull[i];
		const double y = hull[i + 1];

		x_hull.push_back(x);
		y_hull.push_back(y);
	}*/
	for (auto rs : rsLevels) {
		std::vector<delaunay::Segment> anIsoLine = delaun.get_iso_lines(rs, xIsco, yIsco);
		IsoLines.insert(std::make_pair(rs, anIsoLine));
		//std::cout << "For redshift level " << rs << " inserted an isoline of size " << anIsoLine.size() << std::endl;

	}
		

		return  IsoLines;
}

double IsoRedShift::calculateDistance(const delaunay::Point& p1, const delaunay::Point& p2) {
	// Function to calculate the Euclidean distance between two points
	return std::sqrt(std::pow(p2.x - p1.x, 2) + std::pow(p2.y - p1.y, 2));
}

double IsoRedShift::findSmallestDistance(const std::vector<delaunay::Point>& dataset) {
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

void IsoRedShift::improve() {
	;
}