#ifndef BLACK_SCHOLES_H
#define BLACK_SCHOLES_H

double normal_cdf(double x);

double black_scholes_call(double S0, double K, double r, double q, double sigma, double T);
double black_scholes_put (double S0, double K, double r, double q, double sigma, double T);

#endif