#pragma once
#include "vector.h"
#include <string>

namespace pp {

struct matrix {
private:
    int nrows_ = 0;
    int ncols_ = 0;
    std::vector<double> data_; 

public:
    matrix() = default;
    matrix(int n, int m);

    int size1() const; // rows
    int size2() const; // columns

    double& operator[](int i, int j);
    const double& operator[](int i, int j) const;

    double& operator()(int i, int j);
    const double& operator()(int i, int j) const;

    vector col(int j) const;
    void set_col(int j, const vector& v);

    void setid();
    matrix T() const;
    matrix transpose() const;

    void print(const std::string& s = "") const;

    matrix& operator+=(const matrix& B);
    matrix& operator-=(const matrix& B);
    matrix& operator*=(double c);
    matrix& operator/=(double c);
};

matrix operator+(matrix A, const matrix& B);
matrix operator-(matrix A, const matrix& B);
matrix operator*(matrix A, double c);
matrix operator*(double c, matrix A);
matrix operator/(matrix A, double c);

vector operator*(const matrix& A, const vector& x);
matrix operator*(const matrix& A, const matrix& B);

matrix identity(int n);
matrix diagonal(const vector& w);
matrix random_matrix(int n, int m, double scale = 1.0, unsigned seed = 1);
matrix random_symmetric_matrix(int n, double scale = 1.0, unsigned seed = 1);

bool approx(const matrix& A, const matrix& B, double acc = 1e-9, double eps = 1e-9);

} // namespace pp
