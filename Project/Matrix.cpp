#include <iostream>
#include <cmath>

#include "Matrix.h"

void errorMessage(const std::string& msg)
{
    std::cout << "[ERROR] " << msg << std::endl;
}

// Define Matrix 

Matrix::Matrix() : m_rows(0), m_cols(0) {}

Matrix::Matrix(int rows, int cols, double value)
    : m_rows(rows), m_cols(cols), m(rows, std::vector<double>(cols, value)) {}

int Matrix::rows() const { return m_rows; }
int Matrix::cols() const { return m_cols; }

double& Matrix::operator()(int i, int j) { return m[i][j]; }
double  Matrix::operator()(int i, int j) const { return m[i][j]; }


// Multiply Matrix 

std::vector<double> Matrix::multiply(const std::vector<double>& v) const
{   
    if ((int)v.size() != m_cols)
    {
        errorMessage("Can't multiply - size error");
        return std::vector<double>();
    }

    std::vector<double> res(m_rows, 0.0);
    for (int i = 0; i < m_rows; i++)
    {
        for (int j = 0; j < m_cols; j++)
            res[i] += m[i][j] * v[j];
    }
    return res;
}


// Solve linear system 

std::vector<double> Matrix::solve(const std::vector<double>& b, double eps) const
{
    if (m_rows != m_cols)
    {
        errorMessage("solve: matrix not square");
        return std::vector<double>();
    }

    if ((int)b.size() != m_rows)
    {
        errorMessage("Can't multiply - size error");
        return std::vector<double>();
    }

    int n = m_rows;
    Matrix A = *this;
    std::vector<double> x = b;

    for (int i = 0; i < n; i++)
    {
        if (std::fabs(A(i, i)) < eps)
        {
            errorMessage("zero pivot");
            return std::vector<double>();
        }

        for (int k = i + 1; k < n; k++)
        {
            double factor = A(k, i) / A(i, i);
            for (int j = i; j < n; j++)
                A(k, j) -= factor * A(i, j);
            x[k] -= factor * x[i];
        }
    }

    for (int i = n - 1; i >= 0; i--)
    {
        double sum = x[i];
        for (int j = i + 1; j < n; j++)
            sum -= A(i, j) * x[j];

        if (std::fabs(A(i, i)) < eps)
        {
            errorMessage("zero pivot");
            return std::vector<double>();
        }

        x[i] = sum / A(i, i);
    }

    return x;
}

// Invertibility Test (question 1)

// Check if matrix is invertible

bool Matrix::isInvertible(double eps) const
{
    if (m_rows != m_cols) return false;

    int n = m_rows;
    Matrix A = *this;

    for (int i = 0; i < n; i++)
    {
        if (std::fabs(A(i, i)) < eps) return false; 

        for (int k = i + 1; k < n; k++)
        {
            double factor = A(k, i) / A(i, i);
            for (int j = i; j < n; j++)
                A(k, j) -= factor * A(i, j);
        }
    }

    return true;
}

// Compute the inverse

Matrix Matrix::identity(int n)
{
    Matrix I(n, n, 0.0);
    for (int i = 0; i < n; i++)
        I(i, i) = 1.0;
    return I;
}

Matrix Matrix::inverse(double eps) const
{
    if (m_rows != m_cols)
    {
        errorMessage("Matrix is not square");
        return Matrix();
    }

    int n = m_rows;
    Matrix A = *this;
    Matrix Inv = Matrix::identity(n);

    for (int i = 0; i < n; i++)
    {
        double piv = A(i, i);
        if (std::fabs(piv) < eps)
        {
            errorMessage("Zero pivot (matrix not invertible)");
            return Matrix();
        }

        for (int j = 0; j < n; j++)
        {
            A(i, j) /= piv;
            Inv(i, j) /= piv;
        }

        for (int k = 0; k < n; k++)
        {
            if (k == i) continue;
            double factor = A(k, i);
            for (int j = 0; j < n; j++)
            {
                A(k, j) -= factor * A(i, j);
                Inv(k, j) -= factor * Inv(i, j);
            }
        }
    }

    return Inv;
}

// Test

bool Matrix::Test(bool verbose)
{
    Matrix A(2, 2, 0.0);
    A(0,0) = 2; A(0,1) = 1;
    A(1,0) = 5; A(1,1) = 3;

    if (!A.isInvertible())
    {
        errorMessage("Matrix not invertible");
        return false;
    }

    if (verbose)
    {
        Matrix inv = A.inverse();
        std::cout << "Matrice is invertible" << std::endl;
         std::cout << "Inverse matrix:\n";
        for (int i = 0; i < inv.rows(); ++i)
        {
            for (int j = 0; j < inv.cols(); ++j)
                std::cout << inv(i, j) << " ";
            std::cout << "\n";
        }
    }

    return true;
}
