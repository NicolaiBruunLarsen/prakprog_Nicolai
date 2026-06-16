#pragma once
#include "vector.h"
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <stdexcept>

namespace pp {

struct matrix {
    int nrows = 0;
    int ncols = 0;
    std::vector<double> data;

    matrix() = default;
    matrix(int n, int m) : nrows(n), ncols(m), data(n*m) {}

    int size1() const { return nrows; }
    int size2() const { return ncols; }

    // Professor-style indexing: A[i,j]. Requires C++23.
    double& operator[](int i, int j) { return data.at(i*ncols+j); }
    const double& operator[](int i, int j) const { return data.at(i*ncols+j); }

    // Also keep old-style indexing for convenience.
    double& operator()(int i, int j) { return (*this)[i,j]; }
    const double& operator()(int i, int j) const { return (*this)[i,j]; }

    void setid();
    matrix T() const;
    void print(const std::string& s="") const;
};

matrix identity(int n);
matrix diagonal(const vector& d);
matrix operator+(const matrix& A, const matrix& B);
matrix operator-(const matrix& A, const matrix& B);
matrix operator*(const matrix& A, double c);
matrix operator*(double c, const matrix& A);
matrix operator/(const matrix& A, double c);
vector operator*(const matrix& A, const vector& x);
matrix operator*(const matrix& A, const matrix& B);
bool approx(const matrix& A, const matrix& B, double acc=1e-9, double eps=1e-9);

} // namespace pp
