#include"AccretionDisk.h"

using namespace accretiondisk;

AccretionDisk::AccretionDisk(void) {
	;
}
AccretionDisk::AccretionDisk(const double& mass, const double& incl, const double& innerradius, const double& outerradius, const int& nparticles) :
	M(mass),
	inclination(incl),
	innerRadius(innerradius),
	outerRadius(outerradius),
	nParticles(nparticles) {
	maxFluxPrimary = -1e15;
	minFluxPrimary = 1e15;
	maxFluxSecundary = -1e15;
	minFluxSecundary = 1e15;
	createDisk();
};



/**
===================================================================================================================================
* @brief The AccretionDisk class creates a disk of with nParticles scattered randomely over the inner and the outeredge
*
* @param[in] none
*
*
* @return None: feeds the class members std::vector<Particle> primaryParticles,	std::vector<Particle> secundaryParticles,
		double maxFluxPrimary,  double minFluxPrimary, double maxFluxSecundary,	double minFluxSecundary;
=====================================================================================================================================*/
void AccretionDisk::createDisk() {
	//std::cout << "BHMass = "<< M  <<  "   radii : " << innerRadius << " and " << outerRadius << std::endl;
	primaryParticles = std::vector<Particle>(nParticles);
	secundaryParticles = std::vector<Particle>(nParticles);
	std::random_device rd;  // Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<> disRadius(innerRadius, outerRadius);
	std::uniform_real_distribution<> disPhi(0.0, 2 * M_PI);
	double FsMax = 1.0;// 1.146 / 10000;
	for (int p = 0; p < nParticles; p++) {
		double radius = disRadius(gen);
		double phi = disPhi(gen);
		double alpha = BHphysics::alpha(phi, inclination);
		if (phi > M_PI) { alpha *= -1; }
		double bPrimary = BHphysics::calc_impact_parameter(radius, inclination, alpha, M, solver_params_.midpoint_iterations,
			solver_params_.plot_inbetween, 0, M * solver_params_.min_periastron, solver_params_.initial_guesses, solver_params_.use_ellipse);
		double bSecundary = BHphysics::calc_impact_parameter(radius, inclination, alpha, M, solver_params_.midpoint_iterations,
			solver_params_.plot_inbetween, 1, M * solver_params_.min_periastron, solver_params_.initial_guesses, solver_params_.use_ellipse);
		double velocity = std::sqrt(M / std::pow(radius, 3));
		double primaryRedshiftFactor = BHphysics::redshift_factor(radius, alpha, inclination, M, bPrimary);
		double secundaryRedshiftFactor = BHphysics::redshift_factor(radius, alpha, inclination, M, bSecundary);
		double FluxIntrisic = BHphysics::flux_intrinsic(radius, accretionRate, M) / FsMax;
		double primaryFluxObserved = BHphysics::flux_observed(radius, accretionRate, M, primaryRedshiftFactor) / FsMax;
		double secundaryFluxObserved = BHphysics::flux_observed(radius, accretionRate, M, secundaryRedshiftFactor) / 2;
		std::pair<double, double> xy = OperatorsOrder2::polar_to_cartesian_single(alpha, bPrimary, 0.0);
		double x_P = xy.first;
		double y_P = -xy.second;
		xy = OperatorsOrder2::polar_to_cartesian_single(alpha, bSecundary, 0.0);
		double x_S = xy.first;
		double y_S = xy.second;
		accretiondisk::Particle primary(radius, phi, x_P, y_P, FluxIntrisic, primaryFluxObserved, primaryRedshiftFactor, velocity);
		accretiondisk::Particle secundary(radius, phi, x_S, y_S, FluxIntrisic, secundaryFluxObserved, secundaryRedshiftFactor, velocity);
		primaryParticles[p] = primary;
		secundaryParticles[p] = secundary;
		if (maxFluxPrimary < primaryFluxObserved)maxFluxPrimary = primaryFluxObserved;
		if (minFluxPrimary > primaryFluxObserved)minFluxPrimary = primaryFluxObserved;
		if (maxFluxSecundary < secundaryFluxObserved)maxFluxSecundary = primaryFluxObserved;
		if (minFluxSecundary > secundaryFluxObserved)minFluxSecundary = primaryFluxObserved;
	}
}
double AccretionDisk::reduceAngle(const double& theta_) {
	// Reduce angle to [0, 2*pi)
	double theta = fmod(theta_, 2 * M_PI);

	// Ensure angle is positive
	if (theta < 0) {
		theta += 2 * M_PI;
	}

	// Reduce angle to [0, pi)
	if (theta >= M_PI) {
		theta -= 2 * M_PI;
	}

	return theta;
}
/**
===================================================================================================================================
* @brief Recalculates the accretion disk (mimics rotation of the particles)
*
* @param[in] doubledeltaT: the timelapse, set it to 1
*
*
* @return None: feeds the class members std::vector<Particle> primaryParticles,	std::vector<Particle> secundaryParticles,
		double maxFluxPrimary,  double minFluxPrimary, double maxFluxSecundary,	double minFluxSecundary;
=====================================================================================================================================*/
void AccretionDisk::updateDisk( const double& deltaT) {
	parallel_for(0, nParticles, [&](int p) {
		//for (int p = 0; p < nParticles; p++) {
		double radiusP = primaryParticles[p].r;
		double radiusS = secundaryParticles[p].r;
		double velocityP = std::sqrt(M / std::pow(radiusP, 3));
		double velocityS = std::sqrt(M / std::pow(radiusS, 3));
		double omegaP = std::sqrt(M / std::pow(radiusP, 3)) / radiusP;
		double omegaS = std::sqrt(M / std::pow(radiusS, 3)) / radiusS;
		/*double phiP = reduceAngle(primaryParticles[p].phi + omegaP * deltaT);
		double phiS = reduceAngle(secundaryParticles[p].phi + omegaS * deltaT);*/
		double phiP = primaryParticles[p].phi + omegaP * deltaT;
		double phiS = secundaryParticles[p].phi + omegaS * deltaT;
		double alphaP = BHphysics::alpha(phiP, inclination);
		double alphaS = BHphysics::alpha(phiP, inclination);
		if (phiP > M_PI) { alphaP *= -1; }
		if (phiS > M_PI) { alphaS *= -1; }
		double bPrimary = BHphysics::calc_impact_parameter(radiusP, inclination, alphaP, M, solver_params_.midpoint_iterations,
			solver_params_.plot_inbetween, 0, M * solver_params_.min_periastron, solver_params_.initial_guesses, solver_params_.use_ellipse);
		double bSecundary = BHphysics::calc_impact_parameter(radiusS, inclination, alphaS, M, solver_params_.midpoint_iterations,
			solver_params_.plot_inbetween, 1, M * solver_params_.min_periastron, solver_params_.initial_guesses, solver_params_.use_ellipse);

		double primaryRedshiftFactor = BHphysics::redshift_factor(radiusP, alphaP, inclination, M, bPrimary);
		double secundaryRedshiftFactor = BHphysics::redshift_factor(radiusS, alphaS, inclination, M, bSecundary);
		double FluxIntrisic = BHphysics::flux_intrinsic(radiusP, accretionRate, M);
		double primaryFluxObserved = BHphysics::flux_observed(radiusP, accretionRate, M, primaryRedshiftFactor);
		double secundaryFluxObserved = BHphysics::flux_observed(radiusS, accretionRate, M, secundaryRedshiftFactor) / 2;
		std::pair<double, double> xy = OperatorsOrder2::polar_to_cartesian_single(alphaP, bPrimary, 0.0);
		double x_P = xy.first;
		double y_P = -xy.second;
		xy = OperatorsOrder2::polar_to_cartesian_single(alphaS, bSecundary, 0.0);
		double x_S = xy.first;
		double y_S = xy.second;
		accretiondisk::Particle primary(radiusP, phiP, x_P, y_P, FluxIntrisic, primaryFluxObserved, primaryRedshiftFactor, velocityP);
		accretiondisk::Particle secundary(radiusS, phiS, x_S, y_S, FluxIntrisic, secundaryFluxObserved, secundaryRedshiftFactor, velocityS);
		primaryParticles[p] = primary;
		secundaryParticles[p] = secundary;
		}
	);
}

std::vector<Particle> AccretionDisk::getPrimaryImage() {
	return primaryParticles;
}
std::vector<Particle> AccretionDisk::getSecundaryImage() {
	return secundaryParticles;
}
double AccretionDisk::getMaxFluxPrimary() {
	return maxFluxPrimary;
}
double AccretionDisk::getMinFluxPrimary() {
	return minFluxPrimary;
}
double AccretionDisk::getMaxFluxSecundary() {
	return maxFluxSecundary;
}
double AccretionDisk::getMinFluxSecundary() {
	return minFluxSecundary;
}