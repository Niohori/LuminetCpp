#include "BlackHolePhysics.h"



BHphysics::BHphysics() {
    ;
}


    int BHphysics::find_index_sign_change_indices(const std::vector<double>& y){
        int index = -1;
        std::vector<int> sign;
        for (unsigned i = 0; i < y.size(); i++) {
            int s = 1;
            if (y[i] < 0.0) { s = -1; }
            sign.push_back(s);
        }
        std::vector<int> diff;
        for (unsigned i = 0; i < sign.size()-1; i++) {
            diff.push_back(sign[i+1]-sign[i]);
        }
        for (unsigned i = 0; i < diff.size() - 1; i++) {
            if (diff[i] != 0) { index = i; break; }
        }
        return index;
}




double BHphysics::calc_q(double periastron, double bh_mass, double tol = 1e-3) {
    // Convert Periastron distance P to the variable Q
    // (easier to work with)
    
    // If periastron - 2. * bh_mass < tol:
    //     // limit for small values
    //     return .5 * (periastron - 2. * bh_mass) * (periastron + 6. * bh_mass)
    // If 1/periastron < tol:
    //     // limit for large values
    //     return periastron
    // If periastron <= 2 * bh_mass:
    //     throw std::invalid_argument("Non-physical periastron found (P <= 2M, aka the photon sphere)."
    //                                 "If you want to calculate non-physical values, you should implement the mpmath library");

    // Q is complex if P < 2M = r_s
    return std::sqrt((periastron - 2. * bh_mass) * (periastron + 6. * bh_mass));
}


double BHphysics::calc_b_from_periastron(double periastron, double bh_mass, double tol = 1e-5) {
    // Get impact parameter b from Periastron distance P
    // limits give no substantial speed improvement
    // if (std::abs(periastron) < tol) {  // could physically never happen
    //     std::cout << "tolerance exceeded for calc_b_from_P(P_=" << periastron << ", M=" << bh_mass << ", tol=" << tol << std::endl;
    //     return std::sqrt(3 * periastron * periastron);
    // }
    // WARNING: the paper most definitely has a typo here. The fraction on the right-hand side equals b², not b.
    // Just fill in u_2 in equation 3, and you'll see. Only this way do the limits P -> 3M and P >> M hold true,
    // as well as the value for b_c
    return std::sqrt(periastron * periastron * periastron / (periastron - 2. * bh_mass));  // the impact parameter
}

double BHphysics::k(double periastron, double bh_mass) {
    // Calculate modulus of elliptic integral
    double q = calc_q(periastron, bh_mass);
    // adding limits does not substantially improve speed, nor stability
    // if (q < 10e-3) {  // numerical stability
    //     return std::sqrt(0.5);
    // }
    // WARNING: Paper has an error here. There should be brackets around the numerator.
    return std::sqrt((q - periastron + 6 * bh_mass) / (2 * q));  // the modulus of the elliptic integral
}

double BHphysics::k2(double periastron, double bh_mass, double tol = 1e-6) {
    // Calculate the squared modulus of elliptic integral
    double q = calc_q(periastron, bh_mass);
    // adding limits does not substantially improve speed
    // if (1 / periastron <= tol) {
    //     // limit of P -> inf, Q -> P
    //     return 0.;
    // }
    // WARNING: Paper has an error here. There should be brackets around the numerator.
    return (q - periastron + 6 * bh_mass) / (2 * q);  // the modulus of the ellipitic integral
}

double BHphysics::zeta_inf(double periastron, double bh_mass, double tol = 1e-6) {
    // Calculate Zeta_inf for elliptic integral F(Zeta_inf, k)
    double q = calc_q(periastron, bh_mass);  // Assuming calc_q is implemented elsewhere
    double arg = (q - periastron + 2 * bh_mass) / (q - periastron + 6 * bh_mass);
    double z_inf = std::asin(std::sqrt(arg));
    return z_inf;
}

double BHphysics::zeta_r(double periastron, double r, double bh_mass) {
    // Calculate the elliptic integral argument Zeta_r for a given value of P and r
    double q = calc_q(periastron, bh_mass);  // Assuming calc_q is implemented elsewhere
    double a = (q - periastron + 2 * bh_mass + (4 * bh_mass * periastron) / r) / (q - periastron + (6 * bh_mass));
    double s = std::asin(std::sqrt(a));
    return s;
}

double BHphysics::cos_gamma(double _a, double incl, double tol= 1e-5) {
    // Calculate the cos of the angle gamma
    if (std::abs(incl) < tol) {
        return 0;
    }
    return std::cos(_a) / std::sqrt(std::cos(_a) * std::cos(_a) + 1 / (std::tan(incl) * std::tan(incl)));  // real
}

double BHphysics::cos_alpha(double phi, double incl) {
    // Returns cos(angle) alpha in observer frame given angles phi (black hole frame) and inclination (black hole frame)
    return std::cos(phi) * std::cos(incl) / std::sqrt((1 - std::sin(incl) * std::sin(incl) * std::cos(phi) * std::cos(phi)));
}

double BHphysics::alpha(double phi, double incl) {
    // Returns observer coordinate of photon given phi (BHF) and inclination (BHF)
    return std::acos(cos_alpha(phi, incl));
}

std::vector<double> BHphysics::filter_periastrons(const std::vector<double>& periastron, double bh_mass, double tol) {
    // Removes instances where P == 2*M
    // returns indices where this was the case
    std::vector<double> result;
    for (auto e : periastron) {
        if (std::abs(e - 2. * bh_mass) > tol) {
            result.push_back(e);
        }
    }
    return result;
}


double BHphysics::eq13(double periastron, double ir_radius, double ir_angle, double bh_mass, double incl, int n, double tol) {
    // Relation between radius (where photon was emitted in accretion disk), a and P.
    // P can be converted to b, yielding the polar coordinates (b, a) on the photographic plate

    double z_inf = zeta_inf(periastron, bh_mass, tol);
    double q = calc_q(periastron, bh_mass);
    double m_ = k2(periastron, bh_mass, tol);
    double ell_inf = std::ellint_1(std::sqrt(m_),z_inf);  //incomplete elliptic integral
    double g = std::acos(cos_gamma(ir_angle, incl, tol));

    double ellips_arg;
    if (n) {
        double ell_k = std::ellint_1(std::sqrt(m_),M_PI/2);//complete elliptic integral
        ellips_arg = (g - 2. * n * M_PI) / (2. * std::sqrt(periastron / q)) - ell_inf + 2. * ell_k;
    } else {
        ellips_arg = g / (2. * std::sqrt(periastron / q)) + ell_inf;
    }
    double sn= boost::math::jacobi_sn(std::sqrt(m_), ellips_arg); //Jacobi elliptic function sn TO CHECK:arguments
    double term1 = -(q - periastron + 2. * bh_mass) / (4. * bh_mass * periastron);
    double term2 = ((q - periastron + 6. * bh_mass) / (4. * bh_mass * periastron)) * sn * sn;

    return 1. - ir_radius * (term1 + term2);
}

std::tuple<std::vector<double>, std::vector<double>, int>BHphysics::midpoint_method(
    const std::function<double(double, double, double, double, double, int, double)> func,
    const std::unordered_map<std::string, double>& args,
    const std::vector<double>& x,
    const std::vector<double>& y,
    int index_of_sign_change)
 {
    std::vector<double> new_x = x;
    std::vector<double> new_y = y;

    double x0 = new_x[index_of_sign_change];
    double x1 = new_x[index_of_sign_change + 1];
    double inbetween_x = (x0 + x1) / 2.0;
    new_x.insert(new_x.begin() + index_of_sign_change + 1, inbetween_x);

    double y0 = new_y[index_of_sign_change];
    double y1 = new_y[index_of_sign_change + 1];
    double inbetween_solution =  func(inbetween_x, args.at("ir_radius"), args.at("ir_angle"), args.at("bh_mass"), args.at("incl"), args.at("n"), args.at("tol"));
    new_y.insert(new_y.begin() + index_of_sign_change + 1, inbetween_solution);

    int ind_of_sign_change_ = index_of_sign_change + (y0 * inbetween_solution < 0 ? 0 : 1);
    return std::make_tuple(new_x, new_y, ind_of_sign_change_);
}

double BHphysics::improve_solutions_midpoint(
    const std::function<double(double, double, double, double, double, int, double)>& func,
    const std::unordered_map<std::string, double>& args,
    const std::vector<double>& x,
    const std::vector<double>& y,
    int index_of_sign_change,
    int iterations
) {
    int index_of_sign_change_ = index_of_sign_change;
    std::vector<double> new_x = x;
    std::vector<double> new_y = y;

    for (int iteration = 0; iteration < iterations; ++iteration) {
        auto [updated_x, updated_y, updated_ind] = BHphysics::midpoint_method(func, args, new_x, new_y, index_of_sign_change_);
        new_x = updated_x;
        new_y = updated_y;
        index_of_sign_change_ = updated_ind;
    }

    return new_x[index_of_sign_change_];
}

void BHphysics::get_plot(const std::vector<double>& X, const std::vector<double>& Y, const std::vector<double>& solution, const double& radius) {
    /*fig = plt.figure()
        plt.title("Eq13(P)\nr={}, a={}".format(radius, round(_alpha, 5)))
        plt.xlabel('P')
        plt.ylabel('Eq13(P)')
        plt.axhline(0, color='black')
        plt.plot(X, Y)
        plt.scatter(solution, 0, color='red')
        return plt*/

}
double BHphysics::calc_periastron(double _r, double incl, double _alpha,     double bh_mass, 
    int midpoint_iterations=100, bool plot_inbetween=false, int n=0,    double min_periastron=1, int initial_guesses=20) 
{
    /*
    Given a value for r (BH frame) and alpha (BH/observer frame), calculate the corresponding periastron value
        This periastron can be converted to an impact parameter b, yielding the observer frame coordinates (b, alpha).
        Does this by generating range of periastron values, evaluating eq13 on this range and using a midpoint method
        to iteratively improve which periastron value solves equation 13.
        The considered initial periastron range must not be lower than min_periastron (i.e. the photon sphere),
        otherwise non-physical solutions will be found. These are interesting in their own right (the equation yields
        complex solutions within radii smaller than the photon sphere!), but are for now outside the scope of this project.
        Must be large enough to include solution, hence the dependency on the radius (the bigger the radius of the
        accretion disk where you want to find a solution, the bigger the periastron solution is, generally)

        Args:
            _r (float): radius on the accretion disk (BH frame)
            incl (float): inclination of the black hole
            _alpha: angle along the accretion disk (BH frame and observer frame)
            bh_mass (float): mass of the black hole
            midpoint_iterations (int): amount of midpoint iterations to do when searching a periastron value solving eq13
            plot_inbetween (bool): plot
            */
   
    std::vector<double> periastron_range = OperatorsOrder2::linspace(min_periastron, 2.0 * _r, initial_guesses);
    std::vector<double> y_values;

    for (double periastron : periastron_range) {
        y_values.push_back(eq13(periastron, _r, _alpha, bh_mass, incl, n, 1e-6));
    }
    auto ind = find_index_sign_change_indices(y_values);
    double periastron_solution = (ind >= 0) ? periastron_range[ind] : NAN;

    if (!std::isnan(periastron_solution)) {
        std::unordered_map<std::string, double> args_eq13 = {{"ir_radius", _r},
                                                  {"ir_angle", _alpha},
                                                  {"bh_mass", bh_mass},
                                                  {"incl", incl},
                                                  {"n", static_cast<double>(n)},
                                                  {"tol", 1e-8}
        };

        periastron_solution = BHphysics::improve_solutions_midpoint(&BHphysics::eq13, args_eq13, periastron_range, y_values, ind, midpoint_iterations);
    }

    if (plot_inbetween) {
        // Implement plotting if necessary
    }

    return periastron_solution;
}

double BHphysics::calc_impact_parameter(double _r, double incl, double _alpha, double bh_mass, int midpoint_iterations, bool plot_inbetween, int n, double min_periastron, int initial_guesses, bool use_ellipse) {
    double periastron_solution = BHphysics::calc_periastron(_r, incl, _alpha, bh_mass, midpoint_iterations, plot_inbetween, n, min_periastron, initial_guesses);

    if (periastron_solution == NAN) {
        // No periastron was found
        throw std::invalid_argument("No solution was found for the periastron.");
    } else if (periastron_solution <= 2. * bh_mass) {
        // Periastron found is non-physical
        if (use_ellipse) {
            return BHphysics::ellipse(_r, _alpha, incl);
        } else {
            throw std::invalid_argument("No physical periastron solution found.");
        }
    } else {
        // Physical periastron found
        return BHphysics::calc_b_from_periastron(periastron_solution, bh_mass);
    }
}




double BHphysics::phi_inf(double periastron, double M) {
    double q = calc_q(periastron, M);
    double ksq = (q - periastron + 6. * M) / (2. * q);
    double z_inf = zeta_inf(periastron, M);
    double phi = 2. * (std::sqrt(periastron / q)) * (std::ellint_1(std::sqrt(ksq), M_PI/2.0) - std::ellint_1(std::sqrt(ksq),z_inf));
    return phi;
}

double BHphysics::mu(double periastron, double bh_mass) {
    return 2 * phi_inf(periastron, bh_mass) - M_PI;
}

double BHphysics::ellipse(double r, double a, double incl) {
    double g = std::acos(cos_gamma(a, incl));
    double b_ = r * sin(g);
    return b_;
}

double BHphysics::flux_intrinsic(double r, double acc, double bh_mass) {
    double r_ = r / bh_mass;
    double log_arg = ((std::sqrt(r_) + std::sqrt(3)) * (std::sqrt(6) - std::sqrt(3))) / ((std::sqrt(r_) - std::sqrt(3)) * (std::sqrt(6) + std::sqrt(3)));
    double f = (3. * bh_mass * acc / (8 * M_PI)) * (1 / ((r_ - 3) * std::pow(r, 2.5))) *
               (std::sqrt(r_) - std::sqrt(6) + std::pow(3, -0.5) * std::log10(log_arg));
    return f;
}

double BHphysics::flux_observed(double r, double acc, double bh_mass, double redshift_factor) {
    double flux_intr = BHphysics::flux_intrinsic(r, acc, bh_mass);
    return flux_intr / pow(redshift_factor, 4);
}
double BHphysics::redshift_factor(double radius, double angle, double incl, double bh_mass, double b_) {
    // WARNING: the paper is absolutely incomprehensible here. Equation 18 for the redshift completely
    // leaves out important factors. It should be:
    // 1 + z = (1 - Ω*b*cos(η)) * (-g_tt -2Ω*g_tϕ - Ω²*g_ϕϕ)^(-1/2)
    // The expressions for the metric components, Ω and the final result of Equation 19 are correct though
    // TODO perhaps implement other metrics? e.g. Kerr, where g_tϕ != 0
    // gff = (radius * np.sin(incl) * np.sin(angle)) ** 2
    // gtt = - (1 - (2. * M) / radius)
    double z_factor = (1. + sqrt(bh_mass / pow(radius, 3)) * b_ * sin(incl) * sin(angle)) *
                      pow((1 - 3. * bh_mass / radius), -0.5);
    return z_factor;
}
