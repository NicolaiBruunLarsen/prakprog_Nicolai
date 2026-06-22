#include "min.h"
#include "QR.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <limits>

namespace pp{

vector gradient(objective phi, vector x){
	double phix=phi(x);
	vector g(x.size());
	for(int i=0;i<x.size();i++){
		double dxi=(1+std::abs(x[i]))*std::pow(2.0,-26);
		x[i]+=dxi;
		g[i]=(phi(x)-phix)/dxi;
		x[i]-=dxi;
	}
	return g;
}

matrix hessian(objective phi, vector x){
	int n=x.size();
	matrix H(n,n);
	vector gx=gradient(phi,x);
	for(int j=0;j<n;j++){
		double dxj=(1+std::abs(x[j]))*std::pow(2.0,-13);
		x[j]+=dxj;
		vector dg=gradient(phi,x)-gx;
		for(int i=0;i<n;i++) H(i,j)=dg[i]/dxj;
		x[j]-=dxj;
	}
	return H;
}

vector gradient_central(objective phi, vector x){
	vector g(x.size());
	for(int i=0;i<x.size();i++){
		double dxi=(1+std::abs(x[i]))*std::cbrt(std::numeric_limits<double>::epsilon());
		x[i]+=dxi;
		double fp=phi(x);
		x[i]-=2*dxi;
		double fm=phi(x);
		x[i]+=dxi;
		g[i]=(fp-fm)/(2*dxi);
	}
	return g;
}

matrix hessian_central(objective phi, vector x){
	int n=x.size();
	matrix H(n,n);
	vector dx(n);
	for(int i=0;i<n;i++) dx[i]=(1+std::abs(x[i]))*std::cbrt(std::numeric_limits<double>::epsilon());
	double f0=phi(x);

	for(int i=0;i<n;i++){
		// diagonal element
		x[i]+=dx[i];
		double fp=phi(x);
		x[i]-=2*dx[i];
		double fm=phi(x);
		x[i]+=dx[i];
		H(i,i)=(fp-2*f0+fm)/(dx[i]*dx[i]);

		for(int j=i+1;j<n;j++){
			x[i]+=dx[i]; x[j]+=dx[j];
			double fpp=phi(x);
			x[j]-=2*dx[j];
			double fpm=phi(x);
			x[i]-=2*dx[i];
			double fmm=phi(x);
			x[j]+=2*dx[j];
			double fmp=phi(x);
			x[i]+=dx[i]; x[j]-=dx[j];

			double hij=(fpp-fpm-fmp+fmm)/(4*dx[i]*dx[j]);
			H(i,j)=hij;
			H(j,i)=hij;
		}
	}
	return H;
}

static minimum newton_template(objective phi, vector x, double acc, int max_iter, bool central){
	int steps=0;
	vector g = central ? gradient_central(phi,x) : gradient(phi,x);

	while(g.norm()>acc && steps<max_iter){
		matrix H = central ? hessian_central(phi,x) : hessian(phi,x);
		for(int i=0;i<H.size1();i++) H(i,i)+=1e-6; // Levenberg regularization

		vector dx;
		try{
			dx=QR(H).solve(-g);
		}
		catch(std::exception& e){
			dx=-g; // fallback if Hessian is singular
		}

		if(dot(dx,g)>0) dx=-g; // make sure it is a downhill direction

		double phix=phi(x);
		double lambda=1;
		vector z=x+lambda*dx;

		while(lambda>=1.0/1024){
			z=x+lambda*dx;
			if(phi(z)<phix) break;
			lambda/=2;
		}

		x=z;
		g = central ? gradient_central(phi,x) : gradient(phi,x);
		steps++;
	}

	return {x,phi(x),g,steps};
}

// If the last argument is true, use central finite differences for gradient and Hessian. Otherwise, use forward finite differences.

minimum newton(objective phi, vector x, double acc, int max_iter){
	return newton_template(phi,x,acc,max_iter,false);
}

minimum newton_central(objective phi, vector x, double acc, int max_iter){
	return newton_template(phi,x,acc,max_iter,true);
}

} // namespace pp
