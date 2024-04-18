#ifndef BRENT_H
#define BRENT_H

#include <cmath>
#include <cfloat>

// The return value of Minimize is the minimum of the function f.
// The location where f takes its minimum is returned in the variable minLoc.
// Notation and implementation based on Chapter 5 of Richard Brent's book
// "Algorithms for Minimization Without Derivatives".

template <class TFunction>
double Minimize
(
    TFunction& f,		// [in] objective function to minimize
    double leftEnd,     // [in] smaller value of bracketing interval
    double rightEnd,    // [in] larger value of bracketing interval
    double epsilon    // [in] stopping tolerance
)
{
   
    double d, e, m, p, q, r, tol, t2, u, v, w, fu, fv, fw, fx;
    static const double c = 0.5*(3.0 - sqrt(5.0));
    static const double SQRT_DBL_EPSILON = sqrt(DBL_EPSILON);
    
    double& a = leftEnd; double& b = rightEnd; 
    double x;

    v = w = x = a + c*(b - a); d = e = 0.0;
    fv = fw = fx = f(x);
	int counter = 0;
loop:
	counter++;
    m = 0.5*(a + b);
    tol = SQRT_DBL_EPSILON*fabs(x) + epsilon; t2 = 2.0*tol;
    // Check stopping criteria
    if (fabs(x - m) > t2 - 0.5*(b - a))
    {
        p = q = r = 0.0;
        if (fabs(e) > tol)
        {
            // fit parabola
            r = (x - w)*(fx - fv);
            q = (x - v)*(fx - fw);
            p = (x - v)*q - (x - w)*r;
            q = 2.0*(q - r);
            (q > 0.0) ? p = -p : q = -q;
            r = e; e = d;
        }
        if (fabs(p) < fabs(0.5*q*r) && p < q*(a - x) && p < q*(b - x))
        {
            // A parabolic interpolation step
            d = p/q;
            u = x + d;
            // f must not be evaluated too close to a or b
            if (u - a < t2 || b - u < t2)
                d = (x < m) ? tol : -tol;
        }
        else
        {
            // A golden section step
            e = (x < m) ? b : a;
            e -= x;
            d = c*e;
        }
        // f must not be evaluated too close to x
        if (fabs(d) >= tol)
            u = x + d;
        else if (d > 0.0)
            u = x + tol;
        else
            u = x - tol;
        fu = f(u);
        // Update a, b, v, w, and x
        if (fu <= fx)
        {
            (u < x) ? b = x : a = x;
            v = w; fv = fw; 
            w = x; fw = fx; 
            x = u; fx = fu;
        }
        else
        {
            (u < x) ? a = u : b = u;
            if (fu <= fw || w == x)
            {
                v = w; fv = fw; 
                w = u; fw = fu;
            }
            else if (fu <= fv || v == x || v == w)
            {
                v = u; fv = fu;
            }
        }
        goto loop;  // Yes, the dreaded goto statement. But the code here is faithful to Brent's orginal pseudocode.
    }
    return  x;
}

#endif