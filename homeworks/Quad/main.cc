#include "QUAD.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>

const double pi = 3.141592653589793238462643383279502884;

int ncalls;

double erf_integral(double z, double acc=0.001, double eps=0.001)
{
	if(z<0) return -erf_integral(-z,acc,eps);

	if(z<=1){
		std::function<double(double)> f = [](double x){return std::exp(-x*x);};
		return 2/std::sqrt(pi)*integrate(f,0,z,acc,eps);
	}
	else{
		std::function<double(double)> f = [=](double t){
			if(t==0) return 0.0;
			double x = z + (1-t)/t;
			return std::exp(-x*x)/(t*t);
		};
		return 1 - 2/std::sqrt(pi)*integrate(f,0,1,acc,eps);
	}
}

void write_erf_data()
{
	std::ofstream data("out_erf.txt");
	for(double z=-3; z<=3.0001; z+=0.05){
		data << z << " " << erf_integral(z,1e-5,1e-5) << " " << std::erf(z) << "\n";
	}
}

void write_erf_convergence()
{
	std::ofstream data("out_erf_convergence.txt");
	double exact = 0.84270079294971486934;
	for(double acc=1e-1; acc>=1e-7; acc/=10){
		double y = erf_integral(1,acc,0);
		data << acc << " " << std::abs(y-exact) << "\n";
	}
}


int main()
{
	std::cout << std::setprecision(16);

	std::cout << "Part A:\n";
	std::cout << "integral                         result             exact              ncalls\n";

	ncalls = 0;
	std::function<double(double)> f_sqrt = [](double x){ncalls++; return std::sqrt(x);};
	double Q_sqrt = integrate(f_sqrt,0,1,1e-6,1e-6);
	std::cout << "int_0^1 sqrt(x) dx              " << Q_sqrt << "  " << 2.0/3.0 << "  " << ncalls << "\n";

	ncalls = 0;
	std::function<double(double)> f_inv_sqrt = [](double x){ncalls++; return 1/std::sqrt(x);};
	double Q_inv_sqrt = integrate(f_inv_sqrt,0,1,1e-4,1e-4);
	std::cout << "int_0^1 1/sqrt(x) dx            " << Q_inv_sqrt << "  " << 2.0 << "  " << ncalls << "\n";

	ncalls = 0;
	std::function<double(double)> f_circle = [](double x){ncalls++; return std::sqrt(1-x*x);};
	double Q_circle = integrate(f_circle,0,1,1e-6,1e-6);
	std::cout << "int_0^1 sqrt(1-x^2) dx          " << Q_circle << "  " << pi/4 << "  " << ncalls << "\n";

	ncalls = 0;
	std::function<double(double)> f_log = [](double x){ncalls++; return std::log(x)/std::sqrt(x);};
	double Q_log = integrate(f_log,0,1,1e-4,1e-4);
	std::cout << "int_0^1 log(x)/sqrt(x) dx       " << Q_log << "  " << -4.0 << "  " << ncalls << "\n\n";

	std::cout << "erf(1) = " << erf_integral(1,1e-6,0) << "\n";
	std::cout << "tabulated erf(1) = 0.8427007929497149\n\n";

	std::cout << "Part B:\n";
	ncalls = 0;
	std::function<double(double)> g1 = [](double x){ncalls++; return 1/std::sqrt(x);};
	double ordinary = integrate(g1,0,1,1e-4,0);
	int ordinary_calls = ncalls;

	ncalls = 0;
	std::function<double(double)> g2 = [](double x){ncalls++; return 1/std::sqrt(x);};
	double cc = integrate_CC(g2,0,1,1e-4,0);
	int cc_calls = ncalls;
	std::cout << "1/sqrt(x), ordinary:       " << ordinary << "   calls = " << ordinary_calls << "\n";
	std::cout << "1/sqrt(x), Clenshaw-Curtis:" << cc       << "   calls = " << cc_calls << "\n";

	ncalls = 0;
	std::function<double(double)> h1 = [](double x){ncalls++; return std::log(x)/std::sqrt(x);};
	ordinary = integrate(h1,0,1,1e-4,0);
	ordinary_calls = ncalls;

	ncalls = 0;
	std::function<double(double)> h2 = [](double x){ncalls++; return std::log(x)/std::sqrt(x);};
	cc = integrate_CC(h2,0,1,1e-4,0);
	cc_calls = ncalls;
	std::cout << "log(x)/sqrt(x), ordinary:       " << ordinary << "   calls = " << ordinary_calls << "\n";
	std::cout << "log(x)/sqrt(x), Clenshaw-Curtis:" << cc       << "   calls = " << cc_calls << "\n\n";

	std::cout << "Infinite limits using transformation + Clenshaw-Curtis\n";
	double inf = std::numeric_limits<double>::infinity();

	ncalls = 0;
	std::function<double(double)> e1 = [](double x){ncalls++; return std::exp(-x);};
	double Q_inf1 = integrate_CC(e1,0,inf,1e-5,1e-5);
	std::cout << "int_0^inf exp(-x) dx       " << Q_inf1 << "  exact 1, calls = " << ncalls << "\n";

	ncalls = 0;
	std::function<double(double)> e2 = [](double x){ncalls++; return std::exp(-x*x);};
	double Q_inf2 = integrate_CC(e2,-inf,inf,1e-5,1e-5);
	std::cout << "int_-inf^inf exp(-x^2) dx  " << Q_inf2 << "  exact sqrt(pi), calls = " << ncalls << "\n\n";

	std::cout << "\nComparison with Python/scipy\n";

	int ncalls = 0;
	std::function<double(double)> f1 = [&](double x){
		ncalls++;
		return 1/std::sqrt(x);
	};

	ncalls = 0;
	double Qordinary1 = integrate(f1,0,1,1e-4,0);
	int calls_ordinary1 = ncalls;

	ncalls = 0;
	double Qcc1 = integrate_CC(f1,0,1,1e-4,0);
	int calls_cc1 = ncalls;

	std::cout << "Integral of 1/sqrt(x) from 0 to 1\n";
	std::cout << "ordinary integrator: " << Qordinary1 << "  ncalls = " << calls_ordinary1 << "\n";
	std::cout << "Clenshaw-Curtis:     " << Qcc1       << "  ncalls = " << calls_cc1 << "\n";
	std::cout << "Python/scipy:        " << 1.9999999999999984 << "  ncalls = " << 231 << "\n\n";


	std::function<double(double)> f2 = [&](double x){
		ncalls++;
		return std::log(x)/std::sqrt(x);
	};

	ncalls = 0;
	double Qordinary2 = integrate(f2,0,1,1e-4,0);
	int calls_ordinary2 = ncalls;

	ncalls = 0;
	double Qcc2 = integrate_CC(f2,0,1,1e-4,0);
	int calls_cc2 = ncalls;

	std::cout << "Integral of ln(x)/sqrt(x) from 0 to 1\n";
	std::cout << "ordinary integrator: " << Qordinary2 << "  ncalls = " << calls_ordinary2 << "\n";
	std::cout << "Clenshaw-Curtis:     " << Qcc2       << "  ncalls = " << calls_cc2 << "\n";
	std::cout << "Python/scipy:        " << -4.000000000000085 << "  ncalls = " << 315 << "\n\n";



	std::cout << "Part C:\n";
	std::cout << "Integrate log(x)/sqrt(x) from 0 to 1\n";
	double Qe,err;
	std::function<double(double)> ferr = [](double x){return std::log(x)/std::sqrt(x);};
	std::tie(Qe,err) = integrate_with_err(ferr,0,1,1e-4,0);
	std::cout << "result = " << Qe << "  estimated error = " << err << "  true error = " << std::abs(Qe+4) << "\n";

	write_erf_data();
	write_erf_convergence();

	return 0;
}
