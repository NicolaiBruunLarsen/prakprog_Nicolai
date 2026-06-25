#pragma once
#include "vector.h"
#include <functional>
#include <vector>

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

	// Ordinary training: optimize all parameters at once.
	int train(vector xs, vector ys, double acc=1e-4, int max_iter=30);

	// SGD-like/block training: optimize selected parameters, one block at a time.
	int train_indices(vector xs, vector ys, std::vector<int> indices, double acc=1e-4, int max_iter=10);
	int train_neurons(vector xs, vector ys, int sweeps=20, double acc=1e-4, int max_iter=10);
	int train_parameter_sets(vector xs, vector ys, int sweeps=20, double acc=1e-4, int max_iter=10);
	int train_random(vector xs, vector ys, int sweeps=100, int nsubset=3, double acc=1e-4, int max_iter=10, unsigned seed=1);
};

} // namespace pp
