#pragma once
#include "matrix.h"
#include <functional>

namespace pp{

using objective = std::function<double(vector)>;

struct minimum{
	vector x;
	double fx;
	vector grad;
	int steps;
};

vector gradient(objective phi, vector x);
matrix hessian(objective phi, vector x);

vector gradient_central(objective phi, vector x);
matrix hessian_central(objective phi, vector x);

minimum newton(objective phi, vector x, double acc=1e-3, int max_iter=1000);
minimum newton_central(objective phi, vector x, double acc=1e-3, int max_iter=1000);

} // namespace pp
