#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>

#include "Matrix.h"
#include "BlackScholes.h"
#include "PDEPricer.h"

int main()
{
    std::cout << std::fixed << std::setprecision(8);

    // Test Matrix
    bool ok = Matrix::Test(true);
    if (!ok) {
        std::cerr << "Matrix test failed\n";
        return 1;
    }

    std::cout << "\n";

    // Parameters
    PDEParams p;
    p.S0 = 100.0;
    p.K = 100.0;
    p.r = 0.02;
    p.q = 0.00;
    p.sigma = 0.20;
    p.T = 1.0;

    // Boundaries in x
    double width = 10.0 * p.sigma * std::sqrt(p.T);
    p.xmin = -width;
    p.xmax = +width;

    // Black-Scholes price
    double bs = black_scholes_call(p.S0, p.K, p.r, p.q, p.sigma, p.T);
    std::cout << "Black-Scholes Call = " << bs << "\n\n";

    PDEPricer pricer;

    // Values to test and output
    std::vector<int> Ns = {25, 50, 75, 100};
    std::vector<int> Ms = {50, 100, 150, 200};

    std::ofstream fileList("Convergence_results.csv");
    fileList << "n,m,pde,err\n";

    // PDE price and difference with BS
    for (int ni : Ns)
    {
        for (int mj : Ms)
        {
            p.n = ni;
            p.m = mj;

            double pde = pricer.priceCall(p);
            double err = std::fabs(pde - bs);

            std::cout << "n=" << std::setw(4) << p.n
                      << " m=" << std::setw(4) << p.m
                      << " PDE=" << pde
                      << " |err|=" << err << "\n";

            fileList << p.n << "," << p.m << "," << pde << "," << err << "\n";
        }

        std::cout << "----\n";
    }

    fileList.close();
    return 0;
}