#include "Plotter.h"
Plotter::Plotter()
{
	
	screenWidth = GetSystemMetrics(SM_CXSCREEN);;
	screenHeight = GetSystemMetrics(SM_CYSCREEN);

	std::cout << "Screen Resolution: " << screenWidth << "x" << screenHeight << std::endl;

	g = new Dislin;
	x_max = -10000000000.0;
	x_min = 10000000000.0;
	y_max = -10000000000.0;
	y_min = 10000000000.0;
	Npoints = 0;
	g->metafl("cons");
	//g->scrmod("revers");
	g->disini();
	g->pagera();
	g->complx();
	g->nochek();
	g->axspos(500, 0*screenHeight-200);
	//g->axspos(1250, 1600);
	int scrheight = std::min(screenWidth, screenHeight)-20;
	g->axslen(scrheight, scrheight);
	//g->axslen(1200, 2800);
	g->chacod("ISO1");
	g->winfnt("Times New Roman");
	ic = g->intrgb(0., 0., 0.);
	g->axsbgd(ic);
	//g->setclr(0);
}
Plotter::~Plotter() {
	// Destructor implementation
	g->disfin();
}


 
void Plotter::plot_iso_redshifts(const double& inclination, const std::multimap<double, std::vector<delaunay::Segment > >& isolines, const double& x_max_, const double& x_min_, const double& max_rs, const double& min_rs, const bool& loop) {
	//Dislin g;
	double N_redshifts = isolines.size();
	// Counting the number of different keys
	std::vector<double> uniqueKeys;

	for (const auto& pair : isolines) {
		uniqueKeys.push_back(pair.first);
	}

	//std::cout << "xmax" << x_max << "  xmin "<< x_min << std::endl;

	rgbVector = convertToRGB(normalize_vector(uniqueKeys));
	x_max = std::min(x_max_,60.0);
	
	//x_max = 1.0;
	x_min = -x_max;
	y_max = x_max;
	y_min = x_min;
	//x_max = std::max(std::abs(x_max), 40.0);
	//x_max = 40.0;
	x_min = -x_max;
	y_max = x_max;
	y_min = -x_max;
	if (loop) {
		g->erase();
		x_max = 50.0;
		x_min = -x_max;
		y_max = x_max;
		y_min = x_min;
		//std::cout << x_min << "   "<< x_max << std::endl;
	};
	g->titlin("Redshift", 1);
	// Convert double to string with 2-digit precision
	std::string inclination_as_string = std::to_string(inclination);
	size_t dotPos = inclination_as_string.find('.');
	if (dotPos != std::string::npos && dotPos + 0 < inclination_as_string.size()) {
		inclination_as_string = inclination_as_string.substr(0, dotPos + 0); // keep 2 digits after the dot
	}
	if (inclination_as_string.size() == 1) inclination_as_string = "00" + inclination_as_string;
	if (inclination_as_string.size() == 2) inclination_as_string = "0" + inclination_as_string;
	std::string tittel = "inclination = " + inclination_as_string + static_cast<char>(186);//static_cast<char>('\u00B0');
	g->titlin(&tittel[0], 3);
	g->setclr(0);

	g->graf(x_min, x_max, x_min, (x_max - x_min) / 10, y_min, y_max, y_min, (y_max - y_min) / 10);
	//g->setrgb(0.7, 0.7, 0.7);
	g->setrgb(0.0, 0.0, 0.0);
	g->color("fore");
	g->title();
	// Set the color map based on the color parameter

	g->linwid(3);
	double color = 0.5;
	
	auto clr = g->intrgb(0.0, 0.0, 0.0);
	double color_i = 1.0 / isolines.size();
	size_t cc = 0;
	//plotting the isoradials
	/*for (size_t i = 0; i < x_grid.size(); i++) {
		g.setrgb(0.5, 0.5, 0.5);
		g.rlcirc(x_grid[i], y_grid[i], 0.05);
	}*/

	for (auto iso : isolines) {

		auto segments = iso.second;;
		g->setrgb(cc * color_i, 0, 1.0 - cc * color_i);
		for (auto segment : segments) {

			g->rline(segment.p1.x, segment.p1.y, segment.p2.x, segment.p2.y);
		}
		cc++;
	}
	if (loop) {
		g->endgrf();
		g->sendbf();
	};
	//g.disfin();
}



void Plotter::plot_iso_redshifts(const double& inclination, const std::vector<double>& xVal, const std::vector<double>& yVal, const std::vector<double>& rsVal, const double& x_max_, const double& x_min_, const double& max_rs, const double& min_rs) {
	//Dislin g;
	//std::cout << "xmax, xmin : " << x_max_ << " ," << x_min_ << std::endl;
	//g->disfin();
	std::vector<double> xray;
	std::vector<double> yray;
	std::vector<double> x_grid;
	std::vector<double> y_grid;
	std::vector<std::vector< double> > redshift_grid;
	double xMax = *std::max_element(xVal.begin(), xVal.end());
	double yMax = *std::max_element(yVal.begin(), yVal.end());
	double xMin = *std::min_element(xVal.begin(), xVal.end());
	double yMin = *std::min_element(yVal.begin(), yVal.end());
	double rsMax = *std::max_element(rsVal.begin(), rsVal.end());
	double rsMin = *std::min_element(rsVal.begin(), rsVal.end());
	double grid_resolution =0.9;
	int N_grid = std::ceil((xMax - xMin) / grid_resolution);
	//create the grid
	x_grid = OperatorsOrder2::linspace(xMin, xMin + N_grid * grid_resolution, N_grid);
	y_grid = x_grid;
	redshift_grid = std::vector<std::vector<double> >(N_grid, std::vector<double>(N_grid, 0.0));// instantiated with zero's
	std::cout << "Size of the x_grid: " << y_grid.size() << " and y_grid  " << y_grid.size() << "and redshift_grid : " << redshift_grid.size() << " x " << redshift_grid[0].size() << std::endl;
	//std::cin.get();
	//std::cout << "x_grid.size() = " << x_grid.size() << "redshift_field.size() = " << redshift_grid.size() << std::endl;
	//we now start to fill the vectors x_grid an y_grid with the  (x,y) values of the cloud_points.
	//the redshifts are also store in the right row/col of the matrix redshift_grid
	for (size_t i = 0; i < xVal.size(); i++) {
		int I_x = std::floor((xVal[i] - xMin) / grid_resolution);
		if (I_x < 0)I_x = 0;
		if (I_x > N_grid - 1)I_x = N_grid - 1;
		int I_y = std::floor((yVal[i] - yMin) / grid_resolution);
		if (I_y < 0)I_y = 0;
		if (I_y > N_grid - 1)I_y = N_grid - 1;
		//I_y = N_grid - I_y - 1;//necessary as the isoline routine has a reversed axis
		//std::cout << "Size of the grid: " << x_grid.size() << " x " << y_grid.size() << " and trying to acces element " << I_x << " and " <<I_y << std::endl;
		x_grid[I_x] = xVal[i];
		y_grid[I_y] = yVal[i];
		double rs = rsVal[i]-1;
		//std::cout << rs << std::endl;
		redshift_grid[I_x][I_y] = rs;;
		//std::cout << "@(i,j)= @(" << I_x << ", " << I_y << ")  : " << "(x, y) = (" << x_grid[I_x] << ", " << y_grid[I_y] << ") with redshift = " << redshift_grid[I_x][I_y] << std::endl;
		//if (i > 20)break;;
	}
	int n = xVal.size();
	
	
	// Flatten the matrix into a one-dimensional array
	std::vector<double> zarr;
	for (const auto& row : redshift_grid) {
		zarr.insert(zarr.end(), row.begin(), row.end());
	}
	
	g->scrmod("revers");
	g->metafl("cons");
	g->disini();
	g->pagera();
	g->complx();
	g->nochek();
	g->axspos(500, 0 * screenHeight - 200);
	//g->axspos(1250, 1600);
	int scrheight = std::min(screenWidth, screenHeight) - 20;
	g->axslen(scrheight, scrheight);
	//g->axslen(1200, 2800);
	g->chacod("ISO1");
	g->winfnt("Times New Roman");
	int ic = g->intrgb(1., 1.,1.);
	g->axsbgd(ic);


	x_max = std::min(xMax, 60.0);


	//x_max = 1.0;
	x_min = -x_max;
	y_max = x_max;
	y_min = x_min;
	//x_max = std::max(std::abs(x_max), 40.0);
	//x_max = 40.0;
	x_min = -x_max;
	y_max = x_max;
	y_min = -x_max;
	
	g->titlin("Redshift", 1);
	// Convert double to string with 2-digit precision
	std::string inclination_as_string = std::to_string(inclination);
	size_t dotPos = inclination_as_string.find('.');
	if (dotPos != std::string::npos && dotPos + 0 < inclination_as_string.size()) {
		inclination_as_string = inclination_as_string.substr(0, dotPos + 0); // keep 2 digits after the dot
	}
	if (inclination_as_string.size() == 1) inclination_as_string = "00" + inclination_as_string;
	if (inclination_as_string.size() == 2) inclination_as_string = "0" + inclination_as_string;
	std::string tittel = "inclination = " + inclination_as_string + static_cast<char>(186);//static_cast<char>('\u00B0');
	g->titlin(&tittel[0], 3);
	//g->setclr(0);
	g->intax();
	g->graf(x_min, x_max, x_min, (x_max - x_min) / 10, y_min, y_max, y_min, (y_max - y_min) / 10);
	//g->setrgb(0.7, 0.7, 0.7);
	//g->setrgb(0.0, 0.0, 0.0);
	
	// Set the color map based on the color parameter

	g->linwid(3);
	double color = 0.5;

	
	size_t cc = 0;
	//plotting the isoradials
	/*for (size_t i = 0; i < x_grid.size(); i++) {
		g.setrgb(0.5, 0.5, 0.5);
		g.rlcirc(x_grid[i], y_grid[i], 0.05);
	}*/
	auto zlev= OperatorsOrder2::linspace(rsMin-1, rsMax-1, 5);
	

	for(int i = 0; i < zlev.size(); i++)
	{
		g->setclr(0);
		if (i == 2)
			g->labels("none", "contur");
		else
			g->labels("float", "contur");

		g->contur(&x_grid[0], N_grid, &y_grid[0], N_grid, zarr.data(), zlev[i]);
	}
	g->color("fore");
	g->title();
	
	//g.disfin();
}

void Plotter::plot_isoradials(std::vector<double>& xx, std::vector<double>& yy, std::vector<double>& xx_g, std::vector<double>& yy_g) {
	//Dislin g;
	Npoints = xx.size();

	for (unsigned i = 0; i < Npoints; i++)
	{
		if (xx[i] <= x_min)x_min = xx[i];
		if (xx[i] >= x_max)x_max = xx[i];
		if (yy[i] <= y_min)y_min = yy[i];
		if (yy[i] > y_max)y_max = yy[i];
	}
	x_max *= 1.1;
	x_min *= 1.1;
	y_max *= 1.1;
	y_min *= 1.1;
	if (x_max < y_max)x_max = y_max;
	if (x_min < y_min)y_min = x_min;
	if (x_max > y_max)y_max = x_max;
	if (x_min > y_min)x_min = y_min;

	g->metafl("cons");
	//g->scrmod("revers");
	g->disini();
	g->pagera();
	g->complx();
	g->axspos(450, 1800);
	g->axslen(1200, 1200);

	g->name("X-axis", "x");
	g->name("Y-axis", "y");

	g->labdig(-1, "x");
	g->ticks(9, "x");
	g->ticks(10, "y");

	g->titlin("Test of a bare Isoradial", 1);
	g->titlin("(no redshift)", 3);

	//int ic = g->intrgb(0.95, 0.95, 0.95);
	int ic = g->intrgb(0., 0., 0.);
	g->axsbgd(ic);

	g->graf(x_min, x_max, x_min, (x_max - x_min) / 10, y_min, y_max, y_min, (y_max - y_min) / 10);
	//g->setrgb(0.7, 0.7, 0.7);
	g->setrgb(0.0, 0.0, 0.0);
	g->grid(1, 1);
	g->winfnt("Times New Roman Bold");
	//g->psfont("Courier-Bold");
	g->color("fore");
	g->height(50);
	g->title();

	g->color("red");
	g->linwid(3);
	g->curve(&xx[0], &yy[0], Npoints);
	g->color("blue");
	g->linwid(3);
	g->curve(&xx_g[0], &yy_g[0], Npoints);

	//g->disfin();
}

void Plotter::plot_isoradials(double inclination, std::vector<double>& xx, std::vector<double>& yy, std::vector<double>& xx_g, std::vector<double>& yy_g, std::vector<double>& rs, std::vector<double>& rs_g, bool loop) {
	//Dislin g;
	Npoints = xx.size();

	rgbVector = convertToRGB(normalize_vector(rs));
	rgbVector_g = convertToRGB(normalize_vector(rs_g));

	for (unsigned i = 0; i < Npoints; i++)
	{
		if (xx[i] <= x_min)x_min = xx[i];
		if (xx[i] >= x_max)x_max = xx[i];
		if (yy[i] <= y_min)y_min = yy[i];
		if (yy[i] > y_max)y_max = yy[i];
	}
	x_max *= 1.1;
	x_min *= 1.1;
	y_max *= 1.1;
	y_min *= 1.1;
	/*if (x_max < y_max)x_max = y_max;
	if (x_min < y_min)y_min = x_min;
	if (x_max > y_max)y_max = x_max;
	if (x_min > y_min)x_min = y_min;*/

	x_max = *std::max_element(std::begin({ x_max, y_max }), std::end({ x_max, y_max }));
	x_min = -x_max;
	y_max = x_max;
	y_min = x_min;

	//g->name("X-axis", "x");
	//g->name("Y-axis", "y");

	//g->labdig(-1, "x");
	//g->ticks(9, "x");
	//g->ticks(10, "y");
	if (loop) {
		g->erase();
		x_max = 36.0;
		x_min = -x_max;
		y_max = x_max;
		y_min = x_min;
		//std::cout << x_min << "   "<< x_max << std::endl;
	};

	// Convert double to string with 2-digit precision
	std::string inclination_as_string = std::to_string(inclination);
	size_t dotPos = inclination_as_string.find('.');
	if (dotPos != std::string::npos && dotPos + 0 < inclination_as_string.size()) {
		inclination_as_string = inclination_as_string.substr(0, dotPos + 0); // keep 2 digits after the dot
	}
	if (inclination_as_string.size() == 1) inclination_as_string = "00" + inclination_as_string;
	if (inclination_as_string.size() == 2) inclination_as_string = "0" + inclination_as_string;
	std::string tittel = "inclination = " + inclination_as_string + static_cast<char>(186);//static_cast<char>('\u00B0');
	g->titlin(&tittel[0], 3);
	g->setclr(0);

	g->graf(x_min, x_max, x_min, (x_max - x_min) / 10, y_min, y_max, y_min, (y_max - y_min) / 10);
	//g->setrgb(0.7, 0.7, 0.7);
	g->setrgb(0.0, 0.0, 0.0);
	g->color("fore");
	g->title();
	// Set the color map based on the color parameter
	for (size_t i = 0; i < Npoints; i++) {
		g->setrgb(std::get<0>(rgbVector[i]), std::get<1>(rgbVector[i]), std::get<2>(rgbVector[i]));
		g->rlcirc(xx[i], yy[i], 0.2);
		g->setrgb(std::get<0>(rgbVector_g[i]), std::get<1>(rgbVector_g[i]), std::get<2>(rgbVector_g[i]));
		g->rlcirc(xx_g[i], yy_g[i], 0.2);
	}
	if (loop) {
		g->endgrf();
		g->sendbf();
	};
	//g->disfin();
}

// Function to normalize the redshiftfactor in a range 0 to 100
//stretches it by taking the span of min and max
std::vector<double> Plotter::normalize_vector(std::vector<double> avector) {
	std::vector<double> normalized_vector;
	double max_rs = *std::max_element(std::begin(avector), std::end(avector));
	double min_rs = *std::min_element(std::begin(avector), std::end(avector));
	for (auto value : avector) {
		value = 100 * (value - min_rs) / (max_rs - min_rs);
		//std::cout << value << std::endl;
		normalized_vector.push_back(value);
	}
	return normalized_vector;
}

// Function to convert a value to RGB format
std::vector<std::tuple<double, double, double> > Plotter::convertToRGB(std::vector<double> avector) {
	std::vector<std::tuple<double, double, double> > colors;
	for (auto value : avector) {
		// Map the value to the RGB range
		double red = (value) / 100.;
		double green = (100 - value) / 100;
		double blue = 0.0;
		red = 1;
		green = 1.0 - std::sqrt(2) / 200 * value;
		blue = green;
		colors.push_back(std::make_tuple(red, green, blue));
	}
	return colors;
}

std::vector<double> Plotter::flatten_matrix(const std::vector<std::vector<double> >& amatrix) {
	std::vector<double> zarr;
	for (const auto& row : amatrix) {
		zarr.insert(zarr.end(), row.begin(), row.end());
	}

	return zarr;
}
