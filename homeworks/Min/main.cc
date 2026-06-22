#include "min.h"
#include "QR.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

const double pi = 3.141592653589793238462643383279502884;

using pp::vector;

double rosenbrock(vector x){
	return std::pow(1-x[0],2) + 100*std::pow(x[1]-x[0]*x[0],2);
}

double himmelblau(vector x){
	return std::pow(x[0]*x[0]+x[1]-11,2) + std::pow(x[0]+x[1]*x[1]-7,2);
}

double breit_wigner(double E, double m, double Gamma, double A){
	return A/(std::pow(E-m,2) + Gamma*Gamma/4);
}

int main(){
	std::cout << std::setprecision(12);

	std::cout << "Part A: Newton minimization with numerical gradient and Hessian\n";

	pp::minimum R = pp::newton(rosenbrock,vector{5,10},1e-5);
	std::cout << "Rosenbrock minimum: start at " << vector{5,10} << " -> x = " << R.x << "  f(x) = " << R.fx
	          << "  |grad| = " << R.grad.norm() << "  steps = " << R.steps << "\n";

	std::vector<vector> starts = { {20,20}, {-2,3}, {-4,-3}, {4,-2} };
	std::cout << "Himmelblau minima from different starting points:\n";
	for(auto start : starts){
		pp::minimum H = pp::newton(himmelblau,start,1e-5);
		std::cout << "start " << start << " -> x = " << H.x << "  f(x) = " << H.fx
		          << "  |grad| = " << H.grad.norm() << "  steps = " << H.steps << "\n";
	}

	std::cout << "\nPart B: Higgs Breit-Wigner fit\n";
	std::vector<double> energy,signal,error;
	std::string line;
	double E,s,ds;
	while(std::getline(std::cin,line)){
		if(line.size()==0 || line[0]=='#') continue;
		std::istringstream iss(line);
		if(iss >> E >> s >> ds){
			energy.push_back(E);
			signal.push_back(s);
			error.push_back(ds);
		}
	}

	if(energy.size()==0){
		std::cerr << "No Higgs data on stdin. Run as: ./main < higgs.data.txt\n";
		return 1;
	}

	auto deviation = [&](vector c){
		double m=c[0], Gamma=c[1], A=c[2];
		if(Gamma<=0 || A<=0) return 1e30;
		double sum=0;
		for(size_t i=0;i<energy.size();i++){
			double y=breit_wigner(energy[i],m,Gamma,A);
			sum += std::pow((y-signal[i])/error[i],2);
		}
		return sum;
	};

	// Parameters are c={m,Gamma,A}. These are just reasonable starting values.
	pp::minimum fit = pp::newton(deviation,vector{125,5,50},1e-5,1000);
	double m=fit.x[0], Gamma=std::abs(fit.x[1]), A=fit.x[2];
	std::cout << "m     = " << m     << " GeV\n";
	std::cout << "Gamma = " << Gamma << " GeV\n";
	std::cout << "A     = " << A     << "\n";
	std::cout << "D     = " << fit.fx << "\n";
	std::cout << "steps = " << fit.steps << "\n";

	std::ofstream data("out.higgs.data.txt");
	for(size_t i=0;i<energy.size();i++) data << energy[i] << " " << signal[i] << " " << error[i] << "\n";

	std::ofstream curve("out.higgs.fit.txt");
	for(double x=100;x<=160;x+=0.1) curve << x << " " << breit_wigner(x,m,Gamma,A) << "\n";

	std::cout << "\nPart C: forward versus central finite differences\n";
	pp::minimum Rf = pp::newton(rosenbrock,vector{20,20},1e-5);
	pp::minimum Rc = pp::newton_central(rosenbrock,vector{20,20},1e-5);
	std::cout << "start at " << vector{20,20} << "\n";
	std::cout << "Rosenbrock forward: x = " << Rf.x << "  f = " << Rf.fx << "  steps = " << Rf.steps << "\n";
	std::cout << "Rosenbrock central: x = " << Rc.x << "  f = " << Rc.fx << "  steps = " << Rc.steps << "\n";

	pp::minimum Hf = pp::newton(himmelblau,vector{20,20},1e-5);
	pp::minimum Hc = pp::newton_central(himmelblau,vector{20,20},1e-5);
	std::cout << "start at " << vector{20,20} << "\n";
	std::cout << "Himmelblau forward: x = " << Hf.x << "  f = " << Hf.fx << "  steps = " << Hf.steps << "\n";
	std::cout << "Himmelblau central: x = " << Hc.x << "  f = " << Hc.fx << "  steps = " << Hc.steps << "\n";

	return 0;
}
