#ifndef PDEPRICER_H
#define PDEPRICER_H

#include <vector>

struct PDEParams
{
    double S0;
    double K;
    double r;
    double q;
    double sigma;
    double T;

    double xmin;
    double xmax;

    int n; 
    int m; 
};

class PDEPricer
{
public:
    double priceCall(const PDEParams& p) const;

private:
    // Grids
    std::vector<double> buildTimeGrid(double T, int n) const;
    std::vector<double> buildSpaceGrid(double xmin, double xmax, int m) const;

    // Conditions
    std::vector<double> terminalCall(const PDEParams& p,
                                     const std::vector<double>& xgrid) const;

    void boundaryCall(const PDEParams& p,
                      const std::vector<double>& tgrid,
                      const std::vector<double>& xgrid,
                      std::vector<double>& lower,
                      std::vector<double>& upper) const;
};

#endif
