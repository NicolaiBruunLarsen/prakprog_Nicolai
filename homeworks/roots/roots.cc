#include"roots.h"
#include"QR.h"
#include<cmath>
#include<algorithm>
#include<stdexcept>

namespace pp{

matrix jacobian(root_function f, vector x, vector fx, vector dx){
	int n=x.size();
	if(dx.size()==0){
		dx=vector(n);
		for(int i=0;i<n;i++) dx[i]=std::max(std::abs(x[i]),1.0)*std::pow(2.0,-26);
	}
	matrix J(n,n);
	for(int j=0;j<n;j++){
		x[j]+=dx[j];
		vector df=f(x)-fx;
		for(int i=0;i<n;i++) J(i,j)=df[i]/dx[j];
		x[j]-=dx[j];
	}
	return J;
}

vector newton(root_function f, vector x, double acc, double alpha_min, int max_iter){
	vector fx=f(x);
	for(int iter=0; iter<max_iter; iter++){
		if(fx.norm()<acc) return x;
		matrix J=jacobian(f,x,fx);
		QR QRJ(J);
		vector Dx=QRJ.solve(-fx);

		if(Dx.norm()<std::pow(2.0,-26)*(x.norm()+1)) return x;

		double alpha=1;
		vector z=x+alpha*Dx;
		vector fz=f(z);
		while(fz.norm()>fx.norm() && alpha>alpha_min){
			alpha/=2;
			z=x+alpha*Dx;
			fz=f(z);
		}
		x=z;
		fx=fz;
	}
	return x;
}

vector newton_quadratic(root_function f, vector x, double acc, double alpha_min, int max_iter){
	vector fx=f(x);
	for(int iter=0; iter<max_iter; iter++){
		if(fx.norm()<acc) return x;
		matrix J=jacobian(f,x,fx);
		QR QRJ(J);
		vector Dx=QRJ.solve(-fx);

		if(Dx.norm()<std::pow(2.0,-26)*(x.norm()+1)) return x;

		double phi0=0.5*fx.dot(fx);
		double dphi0=-fx.dot(fx);
		double alpha=1;
		vector z=x+alpha*Dx;
		vector fz=f(z);
		double phia=0.5*fz.dot(fz);

		while(phia>phi0 && alpha>alpha_min){
			double c=(phia-phi0-dphi0*alpha)/(alpha*alpha);
			double alpha_new=-dphi0/(2*c);
			if(alpha_new<alpha/4 || alpha_new>alpha/2 || !std::isfinite(alpha_new)) alpha_new=alpha/2;
			alpha=alpha_new;
			z=x+alpha*Dx;
			fz=f(z);
			phia=0.5*fz.dot(fz);
		}
		x=z;
		fx=fz;
	}
	return x;
}

} // namespace pp
