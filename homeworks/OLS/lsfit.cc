#include "lsfit.h"
#include "QR.h"
#include <cmath>
#include <stdexcept>

namespace pp {

fit_result lsfit(const Funcs& fs, const vector& x, const vector& y, const vector& dy) {
    int n = x.size();
    int m = static_cast<int>(fs.size());

    if (y.size() != n || dy.size() != n) {
        throw std::runtime_error("lsfit: x, y, dy must have the same size");
    }

    if (n < m) {
        throw std::runtime_error("lsfit: need at least as many data points as fit functions");
    }

    matrix A(n,m);
    vector b(n);

    for (int i = 0; i < n; ++i) {
        if (dy[i] <= 0) {
            throw std::runtime_error("lsfit: all uncertainties dy must be positive");
        }

        b[i] = y[i]/dy[i];

        for (int k = 0; k < m; ++k) {
            A[i,k] = fs[k](x[i])/dy[i];
        }
    }

    QR qr(A);
    vector c = qr.solve(b);

    // Covariance matrix: cov = (A^T A)^(-1) = R^(-1) * (R^(-1))^T.
    // Since R is square upper triangular, invert it by solving R*x=e_j.
    matrix Rinv(m,m);

    for (int j = 0; j < m; ++j) {
        vector e(m);
        e[j] = 1.0;
        vector col(m);

        for (int i = m-1; i >= 0; --i) {
            double sum = 0;
            for (int k = i+1; k < m; ++k) {
                sum += qr.R[i,k]*col[k];
            }
            col[i] = (e[i] - sum)/qr.R[i,i];
        }

        Rinv.set_col(j, col);
    }

    matrix cov = Rinv*Rinv.T();

    vector dc(m);
    for (int k = 0; k < m; ++k) {
        dc[k] = std::sqrt(cov[k,k]);
    }

    double chi2 = 0;
    for (int i = 0; i < n; ++i) {
        double residual = (y[i] - evaluate(fs,c,x[i]))/dy[i];
        chi2 += residual*residual;
    }

    fit_result result{c, cov, dc, chi2, n-m};
    return result;
}

std::tuple<vector,matrix> lsfit_tuple(const Funcs& fs, const vector& x, const vector& y, const vector& dy) {
    fit_result result = lsfit(fs,x,y,dy);
    return std::make_tuple(result.c, result.cov);
}

double evaluate(const Funcs& fs, const vector& c, double z) {
    if (static_cast<int>(fs.size()) != c.size()) {
        throw std::runtime_error("evaluate: number of functions and coefficients mismatch");
    }

    double sum = 0;
    for (int k = 0; k < c.size(); ++k) {
        sum += c[k]*fs[k](z);
    }

    return sum;
}

} // namespace pp
