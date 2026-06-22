#pragma once
#include"matrix.h"

namespace pp{

// QR decomposition by modified Gram-Schmidt.
// A is n by m with n >= m. Q is n by m and R is m by m.
struct QR{
	matrix Q;
	matrix R;

	QR(const matrix& A);

	vector qt_mul(const vector& b) const;
	vector solve(const vector& b) const;
	double det() const;
	matrix inverse() const;
};

} // namespace pp
