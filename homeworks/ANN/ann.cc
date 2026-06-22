#include "ann.h"
#include "QR.h"
#include "min.h"
#include <cmath>
#include <stdexcept>
#include <iostream>

// copilot has helped with some of the functions for the ann struct

namespace pp{

ann::ann(int n) : n(n), p(3*n) {}

// Gaussian wavelet and its derivatives.
double ann::f(double x) const{
	return x*std::exp(-x*x);
}

double ann::df(double x) const{
	return (1-2*x*x)*std::exp(-x*x);
}

double ann::ddf(double x) const{
	return (4*x*x*x-6*x)*std::exp(-x*x);
}

// Integral of x*exp(-x*x) is -1/2*exp(-x*x).
double ann::F(double x) const{
	return -0.5*std::exp(-x*x);
}

void ann::initialize(double xmin, double xmax){
	for(int i=0;i<n;i++){
		double t = (n==1) ? 0.5 : double(i)/(n-1);
		p[i]     = xmin + t*(xmax-xmin);        // a_i
		p[n+i]   = (xmax-xmin)/n;               // b_i
		p[2*n+i] = 0;                           // w_i
	}
}

double ann::response(double x) const{
	double sum=0;
	for(int i=0;i<n;i++){
		double a=p[i], b=p[n+i], w=p[2*n+i];
		double z=(x-a)/b;
		sum += w*f(z);
	}
	return sum;
}

double ann::derivative(double x) const{
	double sum=0;
	for(int i=0;i<n;i++){
		double a=p[i], b=p[n+i], w=p[2*n+i];
		double z=(x-a)/b;
		sum += w*df(z)/b;
	}
	return sum;
}

double ann::second_derivative(double x) const{
	double sum=0;
	for(int i=0;i<n;i++){
		double a=p[i], b=p[n+i], w=p[2*n+i];
		double z=(x-a)/b;
		sum += w*ddf(z)/(b*b);
	}
	return sum;
}

double ann::antiderivative(double x) const{
	double sum=0;
	for(int i=0;i<n;i++){
		double a=p[i], b=p[n+i], w=p[2*n+i];
		double z=(x-a)/b;
		sum += w*b*F(z);
	}
	return sum;
}

void ann::fit_weights(vector xs, vector ys){
	matrix A(xs.size(),n);
	for(int k=0;k<xs.size();k++){
		for(int i=0;i<n;i++){
			double z=(xs[k]-p[i])/p[n+i];
			A(k,i)=f(z);
		}
	}
	vector w=QR(A).solve(ys);
	for(int i=0;i<n;i++) p[2*n+i]=w[i];
}

double ann::cost(vector xs, vector ys) const{
	for(int i=0;i<n;i++){
		if(p[n+i]<=1e-4 || !std::isfinite(p[n+i])) return 1e99;
	}
	double sum=0;
	for(int k=0;k<xs.size();k++){
		double r=response(xs[k])-ys[k];
		sum += r*r;
	}
	return sum/xs.size();
}

int ann::train(vector xs, vector ys, double acc, int max_iter){
	initialize(xs[0],xs[xs.size()-1]);
	fit_weights(xs,ys);

	vector p0=p;
	double c0=cost(xs,ys);

	objective C = [this,xs,ys](vector q){
		ann net=*this;
		net.p=q;
		return net.cost(xs,ys);
	};

	try{
		minimum min = newton_central(C,p,acc,max_iter);
		ann test=*this;
		test.p=min.x;
		if(std::isfinite(min.fx) && test.cost(xs,ys) < c0){
			p=min.x;
			return min.steps;
		}
	}
	catch(std::exception& e){
		// Keep the least-squares starting point if Newton fails.
	}

	p=p0;
	return 0;
}

double ann::ode_cost(double a, double b, double c, double yc, double ypc,
                     double alpha, double beta, int m) const{
	for(int i=0;i<n;i++){
		if(p[n+i]<=1e-4 || !std::isfinite(p[n+i])) return 1e99;
	}

	double dx=(b-a)/(m-1);
	double sum=0;
	for(int k=0;k<m;k++){
		double x=a+k*dx;
		double phi=second_derivative(x)+response(x); // example: y''+y=0
		double weight = (k==0 || k==m-1) ? 0.5 : 1.0;
		sum += weight*phi*phi*dx;
	}

	double bc1=response(c)-yc;
	double bc2=derivative(c)-ypc;
	return sum + alpha*bc1*bc1 + beta*bc2*bc2;
}

int ann::train_ode(double a, double b, double c, double yc, double ypc,
                   double alpha, double beta, double acc, int max_iter){
	objective C = [this,a,b,c,yc,ypc,alpha,beta](vector q){
		ann net=*this;
		net.p=q;
		return net.ode_cost(a,b,c,yc,ypc,alpha,beta);
	};

	vector p0=p;
	double c0=C(p);
	try{
		minimum min = newton_central(C,p,acc,max_iter);
		if(std::isfinite(min.fx) && min.fx<c0){
			p=min.x;
			return min.steps;
		}
	}
	catch(std::exception& e){ }
	p=p0;
	return 0;
}

} // namespace pp
