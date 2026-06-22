#pragma once
#include "vector.h"
#include <functional>

namespace pp{

struct ann{
	int n;       // number of hidden neurons
	vector p;    // p = {a0..an-1,b0..bn-1,w0..wn-1}

	ann(int n);

	double f(double x) const;       // gaussian wavelet
	double df(double x) const;
	double ddf(double x) const;
	double F(double x) const;       // anti-derivative of f

	double response(double x) const;
	double derivative(double x) const;
	double second_derivative(double x) const;
	double antiderivative(double x) const;

	void initialize(double xmin, double xmax);
	void fit_weights(vector xs, vector ys);
	double cost(vector xs, vector ys) const;
	int train(vector xs, vector ys, double acc=1e-4, int max_iter=30);

	double ode_cost(double a, double b, double c, double yc, double ypc,
	                double alpha=100.0, double beta=100.0, int m=40) const;
	int train_ode(double a, double b, double c, double yc, double ypc,
	              double alpha=100.0, double beta=100.0, double acc=1e-4, int max_iter=20);
};

} // namespace pp
