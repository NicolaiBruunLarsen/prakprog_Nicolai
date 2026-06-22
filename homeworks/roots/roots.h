#pragma once
#include"matrix.h"
#include<functional>

namespace pp{

using root_function = std::function<vector(vector)>;

matrix jacobian(root_function f, vector x, vector fx, vector dx=vector());
vector newton(root_function f, vector x, double acc=1e-6, double alpha_min=1.0/128, int max_iter=100);
vector newton_quadratic(root_function f, vector x, double acc=1e-6, double alpha_min=1.0/128, int max_iter=100);

} // namespace pp
