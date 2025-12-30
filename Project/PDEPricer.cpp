#include <cmath>

#include "PDEPricer.h"
#include "Matrix.h"


inline double dmax(double a, double b)
{
    return (a > b) ? a : b;
}

// Time and space grids

std::vector<double> PDEPricer::buildTimeGrid(double T, int n) const
{
    std::vector<double> t(n + 1, 0.0);
    double dt = T / n;
    for (int i = 0; i <= n; ++i) t[i] = i * dt;
    return t;
}

std::vector<double> PDEPricer::buildSpaceGrid(double xmin, double xmax, int m) const
{
    std::vector<double> x(m + 1, 0.0);
    double dx = (xmax - xmin) / m;
    for (int j = 0; j <= m; ++j) x[j] = xmin + j * dx;
    return x;
}

// Compute Call final payoff for each point in space 
std::vector<double> PDEPricer::terminalCall(const PDEParams& p,
                                            const std::vector<double>& xgrid) const
{
    std::vector<double> u(xgrid.size(), 0.0);

    const double F  = p.S0 * std::exp((p.r - p.q) * p.T);  
    const double df = std::exp(-p.r * p.T);              

    for (int j = 0; j < (int)xgrid.size(); ++j)
    {
        double S = F * std::exp(xgrid[j]);
        u[j] = df * dmax(S - p.K, 0.0);
    }
    return u;
}

// Set boundaries 

void PDEPricer::boundaryCall(const PDEParams& p,
                            const std::vector<double>& tgrid,
                            const std::vector<double>& xgrid,
                            std::vector<double>& lower,
                            std::vector<double>& upper) const
{
    const int n = (int)tgrid.size() - 1;

    lower.assign(n + 1, 0.0);
    upper.assign(n + 1, 0.0);

    // Left boundary is 0
    for (int i = 0; i <= n; ++i)
        lower[i] = 0.0;

    // Right boundary
    const double xmax = xgrid.back();
    const double uR = p.S0 * std::exp(-p.q * p.T) * std::exp(xmax)
                    - p.K  * std::exp(-p.r * p.T);

    for (int i = 0; i <= n; ++i)
        upper[i] = uR;
}

// Solve the PDE 

double PDEPricer::priceCall(const PDEParams& p) const
{
    const int n = p.n;
    const int m = p.m;

    const double dt = p.T / n;
    const double dx = (p.xmax - p.xmin) / m;

    std::vector<double> tgrid = buildTimeGrid(p.T, n);
    std::vector<double> xgrid = buildSpaceGrid(p.xmin, p.xmax, m);

    std::vector<double> uT = terminalCall(p, xgrid);

    std::vector<double> bcLower, bcUpper;
    boundaryCall(p, tgrid, xgrid, bcLower, bcUpper);

    // Solve for points inside the boundaries 
    const int M = m - 1;
    std::vector<double> U_next(M, 0.0);
    for (int j = 1; j <= m - 1; ++j)
        U_next[j - 1] = uT[j];

    // Coefficients
    const double a = 0.0;
    const double b = -0.5 * p.sigma * p.sigma;
    const double c =  0.5 * p.sigma * p.sigma;
    const double d = 0.0;

    const double Acoeff = c * (1.0/(dx*dx)) - b * (1.0/(2.0*dx));
    const double Bcoeff = a + c * (-2.0/(dx*dx));
    const double Ccoeff = c * (1.0/(dx*dx)) + b * (1.0/(2.0*dx));

    // Matrixes 
    Matrix P(M, M, 0.0);
    Matrix Q(M, M, 0.0);

    for (int k = 0; k < M; ++k)
    {
        if (k - 1 >= 0) P(k, k - 1) = -0.5 * dt * Acoeff;
        P(k, k) = 1.0 - 0.5 * dt * Bcoeff;
        if (k + 1 < M) P(k, k + 1) = -0.5 * dt * Ccoeff;

        if (k - 1 >= 0) Q(k, k - 1) = 0.5 * dt * Acoeff;
        Q(k, k) = 1.0 + 0.5 * dt * Bcoeff;
        if (k + 1 < M) Q(k, k + 1) = 0.5 * dt * Ccoeff;
    }

    // Compute backward
    for (int i = n - 1; i >= 0; --i)
    {
        std::vector<double> rhs = Q.multiply(U_next);

        const double left_i    = bcLower[i];
        const double right_i   = bcUpper[i];
        const double left_ip1  = bcLower[i + 1];
        const double right_ip1 = bcUpper[i + 1];

        rhs[0]     += 0.5 * dt * Acoeff * left_ip1;
        rhs[M - 1] += 0.5 * dt * Ccoeff * right_ip1;

        rhs[0]     -= (-0.5 * dt * Acoeff) * left_i;
        rhs[M - 1] -= (-0.5 * dt * Ccoeff) * right_i;

        U_next = P.solve(rhs);
    }

    // Interpolate to find x=0
 
    int j0 = 0;
    while (j0 < m && xgrid[j0 + 1] < 0.0)
        ++j0;

    double uj, ujp1;

    if (j0 == 0)
        uj = bcLower[0];
    else
        uj = U_next[j0 - 1];

    if (j0 + 1 == m)
        ujp1 = bcUpper[0];
    else
        ujp1 = U_next[j0];

    double xj = xgrid[j0];
    double xjp1 = xgrid[j0 + 1];
    double w = (0.0 - xj) / (xjp1 - xj);

    return (1.0 - w) * uj + w * ujp1;   
}
