#pragma once
#include "matrix.h"
#include <functional>
#include <tuple>
#include <vector>

namespace pp {

using Func = std::function<double(double)>;
using Funcs = std::vector<Func>;

struct fit_result {
    vector c;      // best-fit coefficients
    matrix cov;    // covariance matrix
    vector dc;     // uncertainties sqrt(cov_ii)
    double chi2;   // sum of squared residuals
    int dof;       // degrees of freedom
};

fit_result lsfit(const Funcs& fs, const vector& x, const vector& y, const vector& dy);
std::tuple<vector,matrix> lsfit_tuple(const Funcs& fs, const vector& x, const vector& y, const vector& dy);

double evaluate(const Funcs& fs, const vector& c, double z);

} // namespace pp
