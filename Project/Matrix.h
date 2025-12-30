#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <string>

class Matrix
{
private:
    int m_rows;
    int m_cols;
    std::vector<std::vector<double>> m;

public:
    Matrix();
    Matrix(int rows, int cols, double value = 0.0);

    int rows() const;
    int cols() const;

    double& operator()(int i, int j);
    double operator()(int i, int j) const;


    // Operations for PDE 

    std::vector<double> multiply(const std::vector<double>& v) const;
    std::vector<double> solve(const std::vector<double>& b, double eps = 1e-12) const;


    // Test Invertibility Matrix (question 1)

    static Matrix identity(int n);

    bool isInvertible(double eps = 1e-12) const;
    Matrix inverse(double eps = 1e-12) const;
    
    static bool Test(bool verbose = true);
};

#endif