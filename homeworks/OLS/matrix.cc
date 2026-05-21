#include "matrix.h"
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <random>

namespace pp {

matrix::matrix(int n, int m) : nrows_(n), ncols_(m), data_(n*m) {}

int matrix::size1() const { return nrows_; }
int matrix::size2() const { return ncols_; }

double& matrix::operator[](int i, int j) {
    return data_.at(i*ncols_ + j);
}

const double& matrix::operator[](int i, int j) const {
    return data_.at(i*ncols_ + j);
}

double& matrix::operator()(int i, int j) {
    return (*this)[i,j];
}

const double& matrix::operator()(int i, int j) const {
    return (*this)[i,j];
}

vector matrix::col(int j) const {
    vector v(nrows_);
    for (int i = 0; i < nrows_; ++i) v[i] = (*this)[i,j];
    return v;
}

void matrix::set_col(int j, const vector& v) {
    if (v.size() != nrows_) throw std::runtime_error("set_col: size mismatch");
    for (int i = 0; i < nrows_; ++i) (*this)[i,j] = v[i];
}

void matrix::setid() {
    if (nrows_ != ncols_) throw std::runtime_error("setid: matrix must be square");
    for (int i = 0; i < nrows_; ++i)
        for (int j = 0; j < ncols_; ++j)
            (*this)[i,j] = (i == j ? 1.0 : 0.0);
}

matrix matrix::T() const {
    matrix B(ncols_, nrows_);
    for (int i = 0; i < nrows_; ++i)
        for (int j = 0; j < ncols_; ++j)
            B[j,i] = (*this)[i,j];
    return B;
}

matrix matrix::transpose() const {
    return T();
}

void matrix::print(const std::string& s) const {
    if (!s.empty()) std::cout << s << '\n';
    for (int i = 0; i < nrows_; ++i) {
        for (int j = 0; j < ncols_; ++j) {
            std::cout << std::setw(16) << (*this)[i,j];
        }
        std::cout << '\n';
    }
}

matrix& matrix::operator+=(const matrix& B) {
    if (size1() != B.size1() || size2() != B.size2()) throw std::runtime_error("matrix +=: size mismatch");
    for (int i = 0; i < nrows_; ++i)
        for (int j = 0; j < ncols_; ++j)
            (*this)[i,j] += B[i,j];
    return *this;
}

matrix& matrix::operator-=(const matrix& B) {
    if (size1() != B.size1() || size2() != B.size2()) throw std::runtime_error("matrix -=: size mismatch");
    for (int i = 0; i < nrows_; ++i)
        for (int j = 0; j < ncols_; ++j)
            (*this)[i,j] -= B[i,j];
    return *this;
}

matrix& matrix::operator*=(double c) {
    for (double& x : data_) x *= c;
    return *this;
}

matrix& matrix::operator/=(double c) {
    for (double& x : data_) x /= c;
    return *this;
}

matrix operator+(matrix A, const matrix& B) { return A += B; }
matrix operator-(matrix A, const matrix& B) { return A -= B; }
matrix operator*(matrix A, double c) { return A *= c; }
matrix operator*(double c, matrix A) { return A *= c; }
matrix operator/(matrix A, double c) { return A /= c; }

vector operator*(const matrix& A, const vector& x) {
    if (A.size2() != x.size()) throw std::runtime_error("matrix*vector: size mismatch");
    vector y(A.size1());
    for (int i = 0; i < A.size1(); ++i) {
        double sum = 0;
        for (int j = 0; j < A.size2(); ++j) sum += A[i,j]*x[j];
        y[i] = sum;
    }
    return y;
}

matrix operator*(const matrix& A, const matrix& B) {
    if (A.size2() != B.size1()) throw std::runtime_error("matrix*matrix: size mismatch");
    matrix C(A.size1(), B.size2());
    for (int i = 0; i < C.size1(); ++i) {
        for (int j = 0; j < C.size2(); ++j) {
            double sum = 0;
            for (int k = 0; k < A.size2(); ++k) sum += A[i,k]*B[k,j];
            C[i,j] = sum;
        }
    }
    return C;
}

matrix identity(int n) {
    matrix I(n,n);
    I.setid();
    return I;
}

matrix diagonal(const vector& w) {
    matrix D(w.size(), w.size());
    for (int i = 0; i < w.size(); ++i) D[i,i] = w[i];
    return D;
}

matrix random_matrix(int n, int m, double scale, unsigned seed) {
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dist(-scale, scale);
    matrix A(n,m);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < m; ++j)
            A[i,j] = dist(gen);
    return A;
}

matrix random_symmetric_matrix(int n, double scale, unsigned seed) {
    matrix R = random_matrix(n,n,scale,seed);
    matrix A(n,n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            A[i,j] = 0.5*(R[i,j] + R[j,i]);
        }
    }
    return A;
}

bool approx(const matrix& A, const matrix& B, double acc, double eps) {
    if (A.size1() != B.size1() || A.size2() != B.size2()) return false;
    for (int i = 0; i < A.size1(); ++i)
        for (int j = 0; j < A.size2(); ++j)
            if (!approx(A[i,j], B[i,j], acc, eps)) return false;
    return true;
}

} // namespace pp
