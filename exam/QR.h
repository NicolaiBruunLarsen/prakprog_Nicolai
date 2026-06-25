#pragma once
#include "matrix.h"

namespace pp{

struct QR{
	matrix Q,R;
	QR(const matrix& A);
	vector solve(const vector& b) const;
	matrix inverse() const;
};

} // namespace pp
