#include "QR.h"
#include <cmath>
#include <stdexcept>

namespace pp {

QR::QR(const matrix& A) : Q(A), R(A.size2(), A.size2()) {
    int n = A.size1();
    int m = A.size2();

    if (n < m) {
        throw std::runtime_error("QR: expected tall or square matrix, size1 >= size2");
    }

    // Modified Gram-Schmidt. Works for tall matrices.
    for (int i = 0; i < m; ++i) {
        vector qi = Q.col(i);
        R[i,i] = qi.norm();

        if (approx(R[i,i], 0.0)) {
            throw std::runtime_error("QR: linearly dependent columns");
        }

        qi /= R[i,i];
        Q.set_col(i, qi);

        for (int j = i+1; j < m; ++j) {
            vector qj = Q.col(j);
            R[i,j] = qi.dot(qj);
            qj -= qi*R[i,j];
            Q.set_col(j, qj);
        }
    }
}

vector QR::Qt_times(const vector& b) const {
    if (Q.size1() != b.size()) {
        throw std::runtime_error("QR::Qt_times: size mismatch");
    }

    vector c(Q.size2());

    for (int i = 0; i < Q.size2(); ++i) {
        c[i] = Q.col(i).dot(b);
    }

    return c;
}

vector QR::solve(const vector& b) const {
    if (Q.size1() != b.size()) {
        throw std::runtime_error("QR::solve: size mismatch");
    }

    vector c = Qt_times(b);
    int m = R.size1();
    vector x(m);

    // Back-substitution: R*x = Q^T*b.
    for (int i = m-1; i >= 0; --i) {
        double sum = 0;
        for (int k = i+1; k < m; ++k) {
            sum += R[i,k]*x[k];
        }
        x[i] = (c[i] - sum)/R[i,i];
    }

    return x;
}

matrix QR::inverse() const {
    if (Q.size1() != Q.size2()) {
        throw std::runtime_error("QR::inverse: only square matrices can be inverted");
    }

    int n = Q.size1();
    matrix B(n,n);

    for (int j = 0; j < n; ++j) {
        vector e(n);
        e[j] = 1.0;
        B.set_col(j, solve(e));
    }

    return B;
}

double QR::det() const {
    if (R.size1() != R.size2()) {
        throw std::runtime_error("QR::det: R is not square");
    }

    double d = 1.0;
    for (int i = 0; i < R.size1(); ++i) {
        d *= R[i,i];
    }

    return d;
}

} // namespace pp
