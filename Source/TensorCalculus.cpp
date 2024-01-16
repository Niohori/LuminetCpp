#include "TensorCalculus.h""

OperatorsOrder2::OperatorsOrder2(int nGrid, double delta) :
	nGrid(nGrid),
	delta(delta)
{
	//std::cout << "Setting up 2nd-order derivative operators" << std::endl;
}
OperatorsOrder2::~OperatorsOrder2() {
	// Destructor implementation
}

std::vector<std::vector<std::vector<double>>> OperatorsOrder2::Add(const std::vector<std::vector<std::vector<double>>>& vec1, const double& fac1, const std::vector<std::vector<std::vector<double>>>& vec2, const double& fac2)
{
	std::vector<std::vector<std::vector<double>>> result = vec1;
	for (unsigned i = 0; i < nGrid; i++) {
		for (unsigned j = 0; j < nGrid; j++) {
			for (unsigned k = 0; k < nGrid; k++) {
				result[i][j][k] = fac1 * vec1[i][j][k] + fac2 * vec2[i][j][k];
			}
		}
	}

	return result;
}
std::vector<std::vector<std::vector<double>>> OperatorsOrder2::laplace(const std::vector<std::vector<std::vector<double>>>& fct) {
	int nGrid = this->nGrid;
	double ood2 = 1.0 / (this->delta * this->delta);

	std::vector<std::vector<std::vector<double>>> lap(nGrid, std::vector<std::vector<double>>(nGrid, std::vector<double>(nGrid, 0.0)));

	for (int i = 1; i < nGrid - 1; ++i) {
		for (int j = 1; j < nGrid - 1; ++j) {
			for (int k = 1; k < nGrid - 1; ++k) {
				double dfddx = (fct[i + 1][j][k] - 2.0 * fct[i][j][k] + fct[i - 1][j][k]);
				double dfddy = (fct[i][j + 1][k] - 2.0 * fct[i][j][k] + fct[i][j - 1][k]);
				double dfddz = (fct[i][j][k + 1] - 2.0 * fct[i][j][k] + fct[i][j][k - 1]);
				lap[i][j][k] = ood2 * (dfddx + dfddy + dfddz);
			}
		}
	}

	return lap;
}

std::tuple<std::vector<std::vector<std::vector<double>>>, std::vector<std::vector<std::vector<double>>>, std::vector<std::vector<std::vector<double>>>> OperatorsOrder2::gradient(const std::vector<std::vector<std::vector<double>>>& fct) {
	int nGrid = this->nGrid;
	double oo2d = 0.5 / this->delta;

	std::vector<std::vector<std::vector<double>>> grad_x(nGrid, std::vector<std::vector<double>>(nGrid, std::vector<double>(nGrid, 0.0)));
	std::vector<std::vector<std::vector<double>>> grad_y(nGrid, std::vector<std::vector<double>>(nGrid, std::vector<double>(nGrid, 0.0)));
	std::vector<std::vector<std::vector<double>>> grad_z(nGrid, std::vector<std::vector<double>>(nGrid, std::vector<double>(nGrid, 0.0)));

	for (int i = 1; i < nGrid - 1; ++i) {
		for (int j = 1; j < nGrid - 1; ++j) {
			for (int k = 1; k < nGrid - 1; ++k) {
				grad_x[i][j][k] = oo2d * (fct[i + 1][j][k] - fct[i - 1][j][k]);
				grad_y[i][j][k] = oo2d * (fct[i][j + 1][k] - fct[i][j - 1][k]);
				grad_z[i][j][k] = oo2d * (fct[i][j][k + 1] - fct[i][j][k - 1]);
			}
		}
	}

	return std::make_tuple(grad_x, grad_y, grad_z);
}

std::vector<std::vector<std::vector<double>>> OperatorsOrder2::divergence(const std::vector<std::vector<std::vector<double>>>& v_x,
	const std::vector<std::vector<std::vector<double>>>& v_y,
	const std::vector<std::vector<std::vector<double>>>& v_z) {
	int nGrid = this->nGrid;
	double oo2d = 0.5 / this->delta;

	std::vector<std::vector<std::vector<double>>> div(nGrid, std::vector<std::vector<double>>(nGrid, std::vector<double>(nGrid, 0.0)));

	for (int i = 1; i < nGrid - 1; ++i) {
		for (int j = 1; j < nGrid - 1; ++j) {
			for (int k = 1; k < nGrid - 1; ++k) {
				double v_xx = v_x[i + 1][j][k] - v_x[i - 1][j][k];
				double v_yy = v_y[i][j + 1][k] - v_y[i][j - 1][k];
				double v_zz = v_z[i][j][k + 1] - v_z[i][j][k - 1];
				div[i][j][k] = oo2d * (v_xx + v_yy + v_zz);
			}
		}
	}

	return div;
}

std::tuple<std::vector<std::vector<std::vector<double>>>, std::vector<std::vector<std::vector<double>>>, std::vector<std::vector<std::vector<double>>>> OperatorsOrder2::gradDiv(const std::vector<std::vector<std::vector<double>>>& A_x,
	const std::vector<std::vector<std::vector<double>>>& A_y,
	const std::vector<std::vector<std::vector<double>>>& A_z) {
	int nGrid = this->nGrid;
	double ood2 = 1.0 / (this->delta * this->delta);

	std::vector<std::vector<std::vector<double>>> grad_div_x(nGrid, std::vector<std::vector<double>>(nGrid, std::vector<double>(nGrid, 0.0)));
	std::vector<std::vector<std::vector<double>>> grad_div_y(nGrid, std::vector<std::vector<double>>(nGrid, std::vector<double>(nGrid, 0.0)));
	std::vector<std::vector<std::vector<double>>> grad_div_z(nGrid, std::vector<std::vector<double>>(nGrid, std::vector<double>(nGrid, 0.0)));

	for (int i = 1; i < nGrid - 1; ++i) {
		for (int j = 1; j < nGrid - 1; ++j) {
			for (int k = 1; k < nGrid - 1; ++k) {
				double xx_A_x = A_x[i + 1][j][k] - 2.0 * A_x[i][j][k] + A_x[i - 1][j][k];
				double xy_A_y = A_y[i + 1][j + 1][k] - A_y[i + 1][j - 1][k] - A_y[i - 1][j + 1][k] + A_y[i - 1][j - 1][k];
				double xz_A_z = A_z[i + 1][j][k + 1] - A_z[i + 1][j][k - 1] - A_z[i - 1][j][k + 1] + A_z[i - 1][j][k - 1];
				grad_div_x[i][j][k] = ood2 * (xx_A_x + 0.25 * xy_A_y + 0.25 * xz_A_z);

				double yx_A_x = A_x[i + 1][j + 1][k] - A_x[i + 1][j - 1][k] - A_x[i - 1][j + 1][k] + A_x[i - 1][j - 1][k];
				double yy_A_y = A_y[i][j + 1][k] - 2.0 * A_y[i][j][k] + A_y[i][j - 1][k];
				double yz_A_z = A_z[i][j + 1][k + 1] - A_z[i][j + 1][k - 1] - A_z[i][j - 1][k + 1] + A_z[i][j - 1][k - 1];
				grad_div_y[i][j][k] = ood2 * (0.25 * yx_A_x + yy_A_y + 0.25 * yz_A_z);

				double zx_A_x = A_x[i + 1][j][k + 1] - A_x[i + 1][j][k - 1] - A_x[i - 1][j][k + 1] + A_x[i - 1][j][k - 1];
				double zy_A_y = A_y[i][j + 1][k + 1] - A_y[i][j + 1][k - 1] - A_y[i][j - 1][k + 1] + A_y[i][j - 1][k - 1];
				double zz_A_z = A_z[i][j][k + 1] - 2.0 * A_z[i][j][k] + A_z[i][j][k - 1];
				grad_div_z[i][j][k] = ood2 * (0.25 * zx_A_x + 0.25 * zy_A_y + zz_A_z);
			}
		}
	}

	return std::make_tuple(grad_div_x, grad_div_y, grad_div_z);
}

std::tuple<double, double, double> OperatorsOrder2::partialDerivs(const std::vector<std::vector<std::vector<double>>>& fct, int i, int j, int k) {
	double oo2d = 0.5 / this->delta;
	double fx = 0.0;
	double fy = 0.0;
	double fz = 0.0;

	//compute_f_x
	if (i == 0) {
		fx = oo2d * (-3.0 * fct[i][j][k] + 4.0 * fct[i + 1][j][k] - fct[i + 2][j][k]);
	}
	else if (i == nGrid - 1) {
		fx = oo2d * (3.0 * fct[i][j][k] - 4.0 * fct[i - 1][j][k] + fct[i - 2][j][k]);
	}
	else {
		fx = oo2d * (fct[i + 1][j][k] - fct[i - 1][j][k]);
	}

	//compute_f_y
	if (j == 0) {
		fy = oo2d * (-3.0 * fct[i][j][k] + 4.0 * fct[i][j + 1][k] - fct[i][j + 2][k]);
	}
	else if (j == nGrid - 1) {
		fy = oo2d * (3.0 * fct[i][j][k] - 4.0 * fct[i][j - 1][k] + fct[i][j - 2][k]);
	}
	else {
		fy = oo2d * (fct[i][j + 1][k] - fct[i][j - 1][k]);
	}

	//compute_f_z
	if (k == 0) {
		fz = oo2d * (-3.0 * fct[i][j][k] + 4.0 * fct[i][j][k + 1] - fct[i][j][k + 2]);
	}
	else if (k == nGrid - 1) {
		fz = oo2d * (3.0 * fct[i][j][k] - 4.0 * fct[i][j][k - 1] + fct[i][j][k - 2]);
	}
	else {
		fz = oo2d * (fct[i][j][k + 1] - fct[i][j][k - 1]);
	}

	return std::make_tuple(fx, fy, fz);
}

std::vector<double> OperatorsOrder2::linspace(const double& start_in, const double& end_in, const int& num_in)
{
	std::vector<double> linspaced;

	double start = static_cast<double>(start_in);
	double end = static_cast<double>(end_in);
	double num = static_cast<double>(num_in);

	if (num == 0) { return linspaced; }
	if (num == 1)
	{
		linspaced.push_back(start);
		return linspaced;
	}

	double delta = (end - start) / (num - 1);

	for (int i = 0; i < num - 1; ++i)
	{
		linspaced.push_back(start + delta * i);
	}
	linspaced.push_back(end); // I want to ensure that start and end
	// are exactly the same as the input
	return linspaced;
}

std::pair<std::vector<double>, std::vector<double>> OperatorsOrder2::polar_to_cartesian_lists(std::vector<double> radii, std::vector<double> angles, double rotation) {
	std::vector<double> x;
	std::vector<double> y;

	for (size_t i = 0; i < radii.size(); ++i) {
		double x_ = radii[i] * std::cos(angles[i] + rotation);
		double y_ = radii[i] * std::sin(angles[i] + rotation);
		x.push_back(x_);
		y.push_back(y_);
	}

	return std::make_pair(x, y);
}

std::pair<double, double> OperatorsOrder2::polar_to_cartesian_single(double th, double radius, double rotation) {
	double x = radius * cos(th + rotation);
	double y = radius * sin(th + rotation);
	return std::make_pair(x, y);
}

std::pair<double, double> OperatorsOrder2::cartesian_to_polar(double x, double y) {
	double r = std::sqrt(x * x + y * y);
	double theta = std::atan2(y, x);
	return std::make_pair(theta, r);
}

double OperatorsOrder2::get_angle_around(const std::vector<double>& p1, const std::vector<double>& p2) {
		double cx = p1[0];
		double cy = p1[1];

		std::vector<double> p2_ = { p2[0] - cx, p2[1] - cy };

		double angleCenter = OperatorsOrder2::cartesian_to_polar(cx, cy).first;
		double theta = (angleCenter > M_PI) ? (M_PI - angleCenter) : angleCenter;

		double cosTheta = std::cos(theta);
		double sinTheta = std::sin(theta);

		std::vector<std::vector<double>> rotationMatrix = {
			{cosTheta, -sinTheta},
			{sinTheta, cosTheta}
		};

		std::vector<double> rotatedP2 = OperatorsOrder2::matrixMultiply(rotationMatrix, p2_);
		double angleTargetAroundCenter = OperatorsOrder2::cartesian_to_polar(rotatedP2[0], rotatedP2[1]).first;

		return angleTargetAroundCenter;
	}
std::vector<double> OperatorsOrder2::matrixMultiply(const std::vector<std::vector<double>>& matrix,
	const std::vector<double>& vector) {
	std::vector<double> result(matrix.size(), 0.0);

	for (size_t i = 0; i < matrix.size(); ++i) {
		for (size_t j = 0; j < matrix[i].size(); ++j) {
			result[i] += matrix[i][j] * vector[j];
		}
	}
	return result;
}


/*int main() {
	// Example usage
	VectorXd radii(3);
	radii << 1.0, 2.0, 3.0;

	VectorXd angles(3);
	angles << 0.0, M_PI_2, M_PI;

	auto [x_result, y_result] = polar_to_cartesian_lists(radii, angles, 0.0);

	for (int i = 0; i < x_result.size(); ++i) {
		std::cout << "Point " << i + 1 << ": (" << x_result(i) << ", " << y_result(i) << ")\n";
	}

	return 0;
}*/