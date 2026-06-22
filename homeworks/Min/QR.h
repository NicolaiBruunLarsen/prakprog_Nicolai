#pragma once
#include "matrix.h"

namespace pp {

struct QR {
    matrix Q; // reduced Q, size n x m
    matrix R; // upper triangular R, size m x m

    explicit QR(const matrix& A);

    vector Qt_times(const vector& b) const;
    vector solve(const vector& b) const;
    matrix inverse() const;
    double det() const;
};

} // namespace pp
