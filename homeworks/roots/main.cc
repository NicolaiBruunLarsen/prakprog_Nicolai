#include"roots.h"
#include"ODE.h"
#include<cmath>
#include<fstream>
#include<iomanip>
#include<iostream>

const double pi = 3.141592653589793238462643383279502884;

pp::vector rosen_grad(pp::vector x){
	double X=x[0], Y=x[1];
	return pp::vector{
		-2*(1-X)-400*X*(Y-X*X),
		200*(Y-X*X)
	};
}

pp::vector himmel_grad(pp::vector x){
	double X=x[0], Y=x[1];
	double A=X*X+Y-11;
	double B=X+Y*Y-7;
	return pp::vector{
		4*X*A+2*B,
		2*A+4*Y*B
	};
}

double hydrogen_endpoint(double E, double rmin, double rmax, double acc, double eps){
	std::function<pp::vector(double,pp::vector)> F = [=](double r, pp::vector y){ // c++ lambda function 
		return pp::vector{y[1],-2*(E+1/r)*y[0]};
	};
	pp::vector y0{rmin-rmin*rmin,1-2*rmin};
	auto [rs,ys] = pp::driver(F,rmin,rmax,y0,0.05,acc,eps);
	return ys.back()[0];
}

double find_energy(double rmin, double rmax, double acc, double eps){
	std::function<pp::vector(pp::vector)> M = [=](pp::vector x){ 
		return pp::vector{hydrogen_endpoint(x[0],rmin,rmax,acc,eps)};
	};
	pp::vector E = pp::newton(M,pp::vector{-0.7},1e-5);
	return E[0];
}

void write_wavefunction(double E, double rmin, double rmax, double acc, double eps){
	std::function<pp::vector(double,pp::vector)> F = [=](double r, pp::vector y){
		return pp::vector{y[1],-2*(E+1/r)*y[0]};
	};
	pp::vector y0{rmin-rmin*rmin,1-2*rmin};
	auto [rs,ys] = pp::driver(F,rmin,rmax,y0,0.05,acc,eps);
	std::ofstream out("out_wavefunction.txt");
	for(size_t i=0;i<rs.size();i++){
		double r=rs[i];
		out << r << " " << ys[i][0] << " " << r*std::exp(-r) << "\n";
	}
}

void write_convergence(){
	double rmin=1e-3;
	double rmax=8;
	std::ofstream out_rmax("out_convergence_rmax.txt");
	for(double R : {4.0,5.0,6.0,7.0,8.0,10.0,12.0}){
		double E=find_energy(rmin,R,1e-6,1e-6);
		out_rmax << R << " " << E << " " << std::abs(E+0.5) << "\n";
	}

	std::ofstream out_rmin("out_convergence_rmin.txt");
	for(double r0 : {1e-1,5e-2,2e-2,1e-2,5e-3,1e-3}){
		double E=find_energy(r0,rmax,1e-6,1e-6);
		out_rmin << r0 << " " << E << " " << std::abs(E+0.5) << "\n";
	}

	std::ofstream out_acc("out_convergence_acc.txt");
	for(double acc : {1e-3,1e-4,1e-5,1e-6,1e-7}){
		double E=find_energy(rmin,rmax,acc,acc);
		out_acc << acc << " " << E << " " << std::abs(E+0.5) << "\n";
	}
}

int main(){
	std::cout << std::setprecision(16);

	std::cout << "Part A: Newton method\n";

	std::function<pp::vector(pp::vector)> f1 = [](pp::vector x){
		return pp::vector{x[0]*x[0]-2};
	};
	pp::vector r1=pp::newton(f1,pp::vector{1.0},1e-10);
	std::cout << "root of x^2-2 starting at 1: " << r1 << "  f(x)=" << f1(r1).norm() << "\n";

	std::function<pp::vector(pp::vector)> f2 = [](pp::vector x){
		return pp::vector{x[0]*x[0]+x[1]*x[1]-1, x[0]-x[1]};
	};
	pp::vector r2=pp::newton(f2,pp::vector{0.5,0.9},1e-10);
	std::cout << "root of circle/line system:      " << r2 << "  f(x)=" << f2(r2).norm() << "\n\n";

	pp::vector rosen=pp::newton(rosen_grad,pp::vector{0.0,0.0},1e-8);
	std::cout << "Rosenbrock minimum from gradient: " << rosen << "  |grad|=" << rosen_grad(rosen).norm() << "\n";

	std::cout << "Himmelblau minima from different starts:\n";
	for(pp::vector start : {pp::vector{3,2}, pp::vector{-3,3}, pp::vector{-4,-3}, pp::vector{4,-2}}){
		pp::vector root=pp::newton(himmel_grad,start,1e-8);
		std::cout << "start " << start << " -> " << root << "  |grad|=" << himmel_grad(root).norm() << "\n";
	}

	std::cout << "\nPart B: Hydrogen shooting method\n";
	double rmin=1e-3, rmax=8;
	double E=find_energy(rmin,rmax,1e-6,1e-6);
	std::cout << "rmin=" << rmin << " rmax=" << rmax << " E0=" << E << " exact -0.5, error=" << std::abs(E+0.5) << "\n";
	write_wavefunction(E,rmin,rmax,1e-6,1e-6);
	write_convergence();
	std::cout << "wrote out.wavefunction.txt and convergence files\n";

	std::cout << "\nPart C: quadratic interpolation line-search\n";
	pp::vector rosen_q=pp::newton_quadratic(rosen_grad,pp::vector{0.0,0.0},1e-8,1.0/128,1000);
	std::cout << "Rosenbrock with quadratic line-search: " << rosen_q << "  |grad|=" << rosen_grad(rosen_q).norm() << "\n";

	return 0;
}
