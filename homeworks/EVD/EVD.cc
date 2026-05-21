#include "EVD.h"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <stdexcept>

namespace pp {

void EVD::timesJ(pp::matrix& A, int p, int q, double theta) {
    double c = std::cos(theta), s = std::sin(theta);
    for (int i = 0; i < A.size1(); ++i) {
        double aip = A[i,p], aiq = A[i,q];
        A[i,p] = c*aip - s*aiq;
        A[i,q] = s*aip + c*aiq;
    }
}

void EVD::Jtimes(pp::matrix& A, int p, int q, double theta) {
    double c = std::cos(theta), s = std::sin(theta);
    for (int j = 0; j < A.size2(); ++j) {
        double apj = A[p,j], aqj = A[q,j];
        A[p,j] =  c*apj + s*aqj;
        A[q,j] = -s*apj + c*aqj;
    }
}

EVD::EVD(pp::matrix A) : w(A.size1()), V(A.size1(), A.size1()) {
    if (A.size1() != A.size2()) throw std::runtime_error("EVD: matrix must be square");

    int n = A.size1();
    V.setid();

    bool changed;
    do {
        changed = false;
        ++sweeps;

        for (int p = 0; p < n-1; ++p) {
            for (int q = p+1; q < n; ++q) {
                double apq = A[p,q], app = A[p,p], aqq = A[q,q];
                double theta = 0.5*std::atan2(2*apq, aqq-app);
                double c = std::cos(theta), s = std::sin(theta);
                double new_app = c*c*app - 2*s*c*apq + s*s*aqq;
                double new_aqq = s*s*app + 2*s*c*apq + c*c*aqq;

                if (new_app != app || new_aqq != aqq) {
                    changed = true;
                    ++rotations;
                    timesJ(A,p,q, theta);  // A ← A*J
                    Jtimes(A,p,q,-theta);  // A ← J^T*A
                    timesJ(V,p,q, theta);  // V ← V*J
                }
            }
        }
    } while (changed);

    for (int i = 0; i < n; ++i) w[i] = A[i,i];

    sort_eigenpairs();
}

void EVD::sort_eigenpairs() {
    int n = w.size();
    std::vector<int> index(n);
    std::iota(index.begin(), index.end(), 0);

    std::sort(index.begin(), index.end(), [&](int a, int b) {
        return w[a] < w[b];
    });

    vector w_sorted(n);
    matrix V_sorted(V.size1(), V.size2());

    for (int new_col = 0; new_col < n; ++new_col) {
        int old_col = index[new_col];
        w_sorted[new_col] = w[old_col];
        for (int i = 0; i < n; ++i) V_sorted[i,new_col] = V[i,old_col];
    }

    w = w_sorted;
    V = V_sorted;
}

std::tuple<vector,matrix> jacobi(matrix A) {
    EVD evd(A);
    return std::make_tuple(evd.w, evd.V);
}

void timesJ(matrix& A, int p, int q, double theta) {
    EVD::timesJ(A,p,q,theta);
}

void Jtimes(matrix& A, int p, int q, double theta) {
    EVD::Jtimes(A,p,q,theta);
}

} // namespace pp
