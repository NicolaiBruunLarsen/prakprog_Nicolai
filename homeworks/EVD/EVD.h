#pragma once
#include "matrix.h"
#include <tuple>

namespace pp {

struct EVD {
    pp::vector w;  // eigenvalues
    pp::matrix V;  // eigenvectors as columns

    int sweeps = 0;
    int rotations = 0;

    static void timesJ(pp::matrix& A, int p, int q, double theta);
    static void Jtimes(pp::matrix& A, int p, int q, double theta);

    explicit EVD(pp::matrix A); // A is a copy; the original matrix is not changed

private:
    void sort_eigenpairs();
};

// Optional numpy-like interface mentioned on the homework page.
std::tuple<vector,matrix> jacobi(matrix A);

// Free-function wrappers, if you want to call pp::timesJ and pp::Jtimes directly.
void timesJ(matrix& A, int p, int q, double theta);
void Jtimes(matrix& A, int p, int q, double theta);

} // namespace pp
