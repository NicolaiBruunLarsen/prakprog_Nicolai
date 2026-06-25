#include "ann.h"
#include "QR.h"
#include "min.h"
#include <cmath>
#include <random>
#include <algorithm>
#include <iostream>

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
	catch(std::exception& e){ }

	p=p0;
	return 0;
}

int ann::train_indices(vector xs, vector ys, std::vector<int> indices, double acc, int max_iter){
	vector oldp=p;
	double oldcost=cost(xs,ys);

	vector q0(indices.size());
	for(int i=0;i<q0.size();i++) q0[i]=p[indices[i]];

	objective C = [this,xs,ys,indices,oldp](vector q){
		ann net=*this;
		net.p=oldp;
		for(int i=0;i<q.size();i++) net.p[indices[i]]=q[i];
		return net.cost(xs,ys);
	};

	try{
		minimum min = newton_central(C,q0,acc,max_iter);
		ann test=*this;
		test.p=oldp;
		for(int i=0;i<min.x.size();i++) test.p[indices[i]]=min.x[i];
		double newcost=test.cost(xs,ys);
		if(std::isfinite(newcost) && newcost<oldcost){
			p=test.p;
			return min.steps;
		}
	}
	catch(std::exception& e){ }

	p=oldp;
	return 0;
}

int ann::train_neurons(vector xs, vector ys, int sweeps, double acc, int max_iter){
	initialize(xs[0],xs[xs.size()-1]);
	fit_weights(xs,ys);

	int steps=0;
	for(int sweep=0;sweep<sweeps;sweep++){
		double oldcost=cost(xs,ys);

		for(int i=0;i<n;i++){
			std::vector<int> index = {i,n+i,2*n+i}; // a_i,b_i,w_i
			steps += train_indices(xs,ys,index,acc,max_iter);
		}

		double newcost=cost(xs,ys);
		if(oldcost-newcost < acc*std::max(1.0,std::abs(oldcost))) break;
	}
	return steps;
}

int ann::train_parameter_sets(vector xs, vector ys, int sweeps, double acc, int max_iter){
	initialize(xs[0],xs[xs.size()-1]);
	fit_weights(xs,ys);

	std::vector<int> as,bs,ws;
	for(int i=0;i<n;i++){
		as.push_back(i);
		bs.push_back(n+i);
		ws.push_back(2*n+i);
	}

	int steps=0;
	for(int sweep=0;sweep<sweeps;sweep++){
		double oldcost=cost(xs,ys);

		steps += train_indices(xs,ys,ws,acc,max_iter); // all weights
		steps += train_indices(xs,ys,bs,acc,max_iter); // all widths
		steps += train_indices(xs,ys,as,acc,max_iter); // all centers

		double newcost=cost(xs,ys);
		if(oldcost-newcost < acc*std::max(1.0,std::abs(oldcost))) break;
	}
	return steps;
}

int ann::train_random(vector xs, vector ys, int sweeps, int nsubset, double acc, int max_iter, unsigned seed){
	initialize(xs[0],xs[xs.size()-1]);
	fit_weights(xs,ys);

	std::mt19937 gen(seed);
	std::vector<int> all(p.size());
	for(int i=0;i<p.size();i++) all[i]=i;

	int steps=0;
	double oldcost=cost(xs,ys);
	for(int sweep=0;sweep<sweeps;sweep++){
		std::shuffle(all.begin(),all.end(),gen);

		std::vector<int> index;
		for(int i=0;i<nsubset && i<(int)all.size();i++) index.push_back(all[i]);

		steps += train_indices(xs,ys,index,acc,max_iter);

		double newcost=cost(xs,ys);
		if(std::abs(oldcost-newcost) < acc*std::max(1.0,std::abs(oldcost))) break;
		oldcost=newcost;
	}
	return steps;
}

} // namespace pp
