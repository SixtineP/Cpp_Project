#include "BlackScholes.h"
#include <cmath>

inline double dmax(double a, double b)
{
    return (a > b) ? a : b;
}

double normal_cdf(double x)
{
    return 0.5 * (1.0 + std::erf(x / std::sqrt(2.0)));
}

double black_scholes_call(double S0, double K, double r, double q, double sigma, double T)
{
    if (T <= 0.0)
        return dmax(S0 - K, 0.0);

    if (sigma <= 0.0)
        return dmax(S0 * std::exp(-q * T) - K * std::exp(-r * T), 0.0);

    double sqrtT = std::sqrt(T);
    double d1 = (std::log(S0 / K)
        + (r - q + 0.5 * sigma * sigma) * T)
        / (sigma * sqrtT);

    double d2 = d1 - sigma * sqrtT;

    return S0 * std::exp(-q * T) * normal_cdf(d1)
         - K * std::exp(-r * T) * normal_cdf(d2);
}

double black_scholes_put(double S0, double K, double r, double q, double sigma, double T)
{
    double call = black_scholes_call(S0, K, r, q, sigma, T);
    return call - S0 * std::exp(-q * T) + K * std::exp(-r * T);
}
