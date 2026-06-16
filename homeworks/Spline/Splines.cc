#include"Splines.h"
#include<cassert>
#include<cmath>

namespace pp{

int binsearch(const vec& x, double z){/* locates the interval for z by bisection */
	assert(x.size()>1);
	assert(z>=x[0] && z<=x[x.size()-1]);
	int i=0, j=x.size()-1;
	while(j-i>1){
		int mid=(i+j)/2;
		if(z>x[mid]) i=mid; else j=mid;
		}
	return i;
}

double linterp(const vec& x, const vec& y, double z){
	assert(x.size()==y.size());
	int i=binsearch(x,z);
	double dx=x[i+1]-x[i];
	assert(dx>0);
	double dy=y[i+1]-y[i];
	return y[i]+dy/dx*(z-x[i]);
}

double linterpInteg(const vec& x, const vec& y, double z){
	assert(x.size()==y.size());
	int s=binsearch(x,z);
	double sum=0;
	for(int i=0;i<s;i++){
		double dx=x[i+1]-x[i];
		assert(dx>0);
		double p=(y[i+1]-y[i])/dx;
		sum += y[i]*dx + p*dx*dx/2;
		}
	double h=z-x[s];
	double dx=x[s+1]-x[s];
	assert(dx>0);
	double p=(y[s+1]-y[s])/dx;
	sum += y[s]*h + p*h*h/2;
	return sum;
}

qspline::qspline(const vec& xs, const vec& ys)
	: n(xs.size()), x(xs), y(ys), b(n-1), c(n-1)
{// builds quadratic spline
	assert(n>1 && x.size()==y.size());
	vec h(n-1), p(n-1);
	for(int i=0;i<n-1;i++){
		h[i]=x[i+1]-x[i];
		assert(h[i]>0);
		p[i]=(y[i+1]-y[i])/h[i];
		}
	c[0]=0; // recursion up
	for(int i=0;i<n-2;i++)
		c[i+1]=(p[i+1]-p[i]-c[i]*h[i])/h[i+1];
	c[n-2]/=2; // recursion down
	for(int i=n-3;i>=0;i--)
		c[i]=(p[i+1]-p[i]-c[i+1]*h[i+1])/h[i];
	for(int i=0;i<n-1;i++)
		b[i]=p[i]-c[i]*h[i];
}

int qspline::binsearch(double z) const{
	return pp::binsearch(x,z);
}

double qspline::eval(double z) const{
	int i=binsearch(z);
	double h=z-x[i];
	return y[i]+h*(b[i]+h*c[i]);
}

double qspline::deriv(double z) const{
	int i=binsearch(z);
	double h=z-x[i];
	return b[i]+2*c[i]*h;
}

double qspline::integ(double z) const{
	int s=binsearch(z);
	double sum=0;
	for(int i=0;i<s;i++){
		double h=x[i+1]-x[i];
		sum += y[i]*h + b[i]*h*h/2 + c[i]*h*h*h/3;
		}
	double h=z-x[s];
	sum += y[s]*h + b[s]*h*h/2 + c[s]*h*h*h/3;
	return sum;
}

}//pp
