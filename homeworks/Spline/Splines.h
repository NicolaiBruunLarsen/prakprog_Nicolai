#pragma once
#include"vector.h"

namespace pp{

using vec = vector;

int binsearch(const vec& x, double z);
double linterp(const vec& x, const vec& y, double z);
double linterpInteg(const vec& x, const vec& y, double z);

struct qspline{
	const int n;
	vec x,y,b,c;
	qspline(const vec& x, const vec& y);
	int binsearch(double z) const;
	double eval(double z) const;
	double deriv(double z) const;
	double integ(double z) const;
};

}//pp
