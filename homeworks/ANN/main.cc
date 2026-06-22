#include "ann.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

using pp::vector;

const double pi = 3.141592653589793238462643383279502884;

double exact_fun(double x){
	return std::cos(5*x-1)*std::exp(-x*x);
}

double exact_derivative(double x){
	return std::exp(-x*x)*(-5*std::sin(5*x-1)-2*x*std::cos(5*x-1));
}

double exact_second(double x){
	return std::exp(-x*x)*((4*x*x-27)*std::cos(5*x-1)+20*x*std::sin(5*x-1));
}

int main(){
	std::cout << std::setprecision(12);

	int N=25;
	vector xs(N),ys(N);
	for(int i=0;i<N;i++){
		xs[i]=-1.0 + 2.0*i/(N-1);
		ys[i]=exact_fun(xs[i]);
	}

	pp::ann net(6);
	int steps = net.train(xs,ys,1e-5,25);

	std::cout << "Part A: interpolation with artificial neural network\n";
	std::cout << "number of hidden neurons = " << net.n << "\n";
	std::cout << "training steps = " << steps << "\n";
	std::cout << "cost after training = " << net.cost(xs,ys) << "\n\n";

	std::ofstream train("out.train.txt");
	for(int i=0;i<N;i++) train << xs[i] << " " << ys[i] << "\n";

	std::ofstream out("out.ann.txt");
	double I_exact=0, I_net0=net.antiderivative(-1.5);
	double xold=-1.5, yold=exact_fun(xold);
	for(int i=0;i<=400;i++){
		double x=-1.5+3.0*i/400; // 400 points from -1.5 to 1.5
		double y=exact_fun(x);
		if(i>0) I_exact += 0.5*(y+yold)*(x-xold);
		double I_net=net.antiderivative(x)-I_net0;
		out << x << " "
		    << y << " " << net.response(x) << " "
		    << exact_derivative(x) << " " << net.derivative(x) << " "
		    << exact_second(x) << " " << net.second_derivative(x) << " "
		    << I_exact << " " << I_net << "\n";
		xold=x; yold=y;
	}

	std::cout << "Part C: ANN solution of y''+y=0, y(0)=0, y'(0)=1\n";
	int M=25;
	vector xo(M),yo(M);
	for(int i=0;i<M;i++){
		xo[i]=pi*i/(M-1);
		yo[i]=std::sin(xo[i]); // only used as a starting guess
	}
	pp::ann ode(5);
	ode.train(xo,yo,1e-4,15);
	int odesteps=ode.train_ode(0,pi,0,0,1,100,100,1e-4,15);
	std::cout << "ODE training steps = " << odesteps << "\n";
	std::cout << "ODE cost = " << ode.ode_cost(0,pi,0,0,1,100,100) << "\n";

	std::ofstream odeout("out.ode.txt");
	for(int i=0;i<=400;i++){
		double x=pi*i/400;
		odeout << x << " " << std::sin(x) << " " << ode.response(x)
		       << " " << std::cos(x) << " " << ode.derivative(x)
		       << " " << -std::sin(x) << " " << ode.second_derivative(x) << "\n";
	}

	return 0;
}
