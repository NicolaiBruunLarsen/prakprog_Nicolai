#pragma once
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>
#include <stdexcept>
#include <random>

namespace pp {

// vector
struct vector {
    std::vector<double> data;

    vector() = default;
    explicit vector(int n) : data(n) {}
    vector(std::initializer_list<double> init) : data(init) {}

    int size() const { return static_cast<int>(data.size()); }

    double& operator[](int i) { return data[i]; }
    const double& operator[](int i) const { return data[i]; }

    void print(const std::string& s = "") const {
        if (!s.empty()) std::cout << s << '\n';
        for (double x : data) std::cout << std::setw(12) << x << ' ';
        std::cout << '\n';
    }

    double dot(const vector& other) const {
        assert(size() == other.size());
        double s = 0.0;
        for (int i = 0; i < size(); ++i) s += (*this)[i] * other[i];
        return s;
    }

    double norm() const {
        return std::sqrt(dot(*this));
    }

    vector& operator+=(const vector& other) {
        assert(size() == other.size());
        for (int i = 0; i < size(); ++i) (*this)[i] += other[i];
        return *this;
    }

    vector& operator-=(const vector& other) {
        assert(size() == other.size());
        for (int i = 0; i < size(); ++i) (*this)[i] -= other[i];
        return *this;
    }

    vector& operator*=(double c) {
        for (int i = 0; i < size(); ++i) (*this)[i] *= c;
        return *this;
    }

    vector& operator/=(double c) {
        for (int i = 0; i < size(); ++i) (*this)[i] /= c;
        return *this;
    }
};

inline vector operator+(vector a, const vector& b) { a += b; return a; }
inline vector operator-(vector a, const vector& b) { a -= b; return a; }
inline vector operator-(vector a) { a *= -1; return a; }
inline vector operator*(vector a, double c) { a *= c; return a; }
inline vector operator*(double c, vector a) { a *= c; return a; }
inline vector operator/(vector a, double c) { a /= c; return a; }

inline bool approx(double x, double y, double acc = 1e-9, double eps = 1e-9) {
    double d = std::abs(x - y);
    if (d < acc) return true;
    return d < eps * std::max(std::abs(x), std::abs(y));
}

inline bool approx(const vector& a, const vector& b,
                   double acc = 1e-9, double eps = 1e-9) {
    if (a.size() != b.size()) return false;
    for (int i = 0; i < a.size(); ++i) {
        if (!approx(a[i], b[i], acc, eps)) return false;
    }
    return true;
}


// matrix
// Stored by columns: cols[j][i] = A(i,j)
struct matrix {
    std::vector<vector> cols;

    matrix() = default;
    matrix(int n, int m) : cols(m, vector(n)) {}

    int size1() const { return cols.empty() ? 0 : cols[0].size(); } // rows
    int size2() const { return static_cast<int>(cols.size()); }      // cols

    double& operator()(int i, int j) { return cols[j][i]; }
    const double& operator()(int i, int j) const { return cols[j][i]; }

    vector& operator[](int j) { return cols[j]; }              // j-th column
    const vector& operator[](int j) const { return cols[j]; }

    static matrix id(int n) {
        matrix I(n, n);
        for (int i = 0; i < n; ++i) I(i, i) = 1.0;
        return I;
    }

    matrix transpose() const {
        matrix T(size2(), size1());
        for (int i = 0; i < size1(); ++i)
            for (int j = 0; j < size2(); ++j)
                T(j, i) = (*this)(i, j);
        return T;
    }

    void print(const std::string& s = "") const {
        if (!s.empty()) std::cout << s << '\n';
        for (int i = 0; i < size1(); ++i) {
            for (int j = 0; j < size2(); ++j)
                std::cout << std::setw(12) << (*this)(i, j) << ' ';
            std::cout << '\n';
        }
    }

    matrix& operator+=(const matrix& B) {
        assert(size1() == B.size1() && size2() == B.size2());
        for (int j = 0; j < size2(); ++j) cols[j] += B[j];
        return *this;
    }

    matrix& operator-=(const matrix& B) {
        assert(size1() == B.size1() && size2() == B.size2());
        for (int j = 0; j < size2(); ++j) cols[j] -= B[j];
        return *this;
    }

    matrix& operator*=(double c) {
        for (int j = 0; j < size2(); ++j) cols[j] *= c;
        return *this;
    }

    matrix& operator/=(double c) {
        for (int j = 0; j < size2(); ++j) cols[j] /= c;
        return *this;
    }
};

inline matrix operator+(matrix A, const matrix& B) { A += B; return A; }
inline matrix operator-(matrix A, const matrix& B) { A -= B; return A; }
inline matrix operator*(matrix A, double c) { A *= c; return A; }
inline matrix operator*(double c, matrix A) { A *= c; return A; }
inline matrix operator/(matrix A, double c) { A /= c; return A; }

inline vector operator*(const matrix& A, const vector& x) {
    assert(A.size2() == x.size());
    vector y(A.size1());
    for (int j = 0; j < A.size2(); ++j) {
        double xj = x[j];
        for (int i = 0; i < A.size1(); ++i)
            y[i] += A(i, j) * xj;
    }
    return y;
}

inline matrix operator*(const matrix& A, const matrix& B) {
    assert(A.size2() == B.size1());
    matrix C(A.size1(), B.size2());
    for (int j = 0; j < B.size2(); ++j)
        C[j] = A * B[j];
    return C;
}

inline bool approx(const matrix& A, const matrix& B,
                   double acc = 1e-9, double eps = 1e-9) {
    if (A.size1() != B.size1() || A.size2() != B.size2()) return false;
    for (int i = 0; i < A.size1(); ++i)
        for (int j = 0; j < A.size2(); ++j)
            if (!approx(A(i, j), B(i, j), acc, eps)) return false;
    return true;
}

// qr decomposition by modified Gram-Schmidt
// A is n x m, n >= m
// Q is n x m, R is m x m
struct qr {
    matrix Q;
    matrix R;

    explicit qr(const matrix& A) : Q(A), R(A.size2(), A.size2()) {
        int n = A.size1();
        int m = A.size2();
        assert(n >= m);

        for (int i = 0; i < m; ++i) {
            R(i, i) = Q[i].norm();
            if (approx(R(i, i), 0.0)) {
                throw std::runtime_error("qr: linearly dependent columns / singular matrix");
            }

            Q[i] /= R(i, i);

            for (int j = i + 1; j < m; ++j) {
                R(i, j) = Q[i].dot(Q[j]);
                Q[j] -= Q[i] * R(i, j);
            }
        }
    }

    vector qt_mul(const vector& b) const {
        assert(Q.size1() == b.size());
        vector c(Q.size2());
        for (int i = 0; i < Q.size2(); ++i)
            c[i] = Q[i].dot(b);
        return c;
    }

    vector solve(const vector& b) const {
        assert(Q.size1() == b.size());
        assert(R.size1() == R.size2());

        vector c = qt_mul(b);
        int m = R.size1();
        vector x(m);

        for (int i = m - 1; i >= 0; --i) {
            double s = 0.0;
            for (int k = i + 1; k < m; ++k)
                s += R(i, k) * x[k];
            x[i] = (c[i] - s) / R(i, i);
        }
        return x;
    }

    double det() const {
        assert(R.size1() == R.size2());
        double d = 1.0;
        for (int i = 0; i < R.size1(); ++i) d *= R(i, i);
        return d;
    }

    matrix inverse() const {
        assert(Q.size1() == Q.size2()); // only square matrices are invertible
        int n = Q.size1();
        matrix B(n, n);

        for (int j = 0; j < n; ++j) {
            vector e(n);
            e[j] = 1.0;
            B[j] = solve(e); // j-th column of inverse
        }
        return B;
    }
};

// helpers for testing
inline matrix random_matrix(int n, int m, double scale = 1.0, unsigned seed = 1) {
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dist(-scale, scale);
    matrix A(n, m);
    for (int j = 0; j < m; ++j)
        for (int i = 0; i < n; ++i)
            A(i, j) = dist(gen);
    return A;
}

inline vector random_vector(int n, double scale = 1.0, unsigned seed = 2) {
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dist(-scale, scale);
    vector v(n);
    for (int i = 0; i < n; ++i) v[i] = dist(gen);
    return v;
}

inline bool is_upper_triangular(const matrix& R, double tol = 1e-10) {
    for (int i = 0; i < R.size1(); ++i)
        for (int j = 0; j < i && j < R.size2(); ++j)
            if (std::abs(R(i, j)) > tol) return false;
    return true;
}

}