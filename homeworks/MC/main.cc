#include"MC.h"
#include<cmath>
#include<fstream>
#include<iomanip>
#include<iostream>
#include<vector>

const double pi = 3.141592653589793238462643383279502884;

void write_circle_errors()
{
	std::ofstream data("out.circle.txt");
	vec a = {-1,-1};
	vec b = { 1, 1};
	func circle = [](vec x){return x[0]*x[0]+x[1]*x[1] <= 1 ? 1.0 : 0.0;};

	for(int N: {100,300,1000,3000,10000,30000,100000}){
		double Q,err;
		std::tie(Q,err) = plainmc_lcg(circle,a,b,N,17);
		data << N << " " << err << " " << std::abs(Q-pi) << " " << Q << "\n";
	}
}

void write_quasi_errors()
{
	std::ofstream data("out.quasi.txt");
	vec a = {-1,-1};
	vec b = { 1, 1};
	func circle = [](vec x){return x[0]*x[0]+x[1]*x[1] <= 1 ? 1.0 : 0.0;};

	for(int N: {100,300,1000,3000,10000,30000,100000}){
		double Qp,errp,Qq,errq;
		std::tie(Qp,errp) = plainmc_lcg(circle,a,b,N,19);
		std::tie(Qq,errq) = quasymc(circle,a,b,N);
		data << N << " " << std::abs(Qp-pi) << " " << errp << " " << std::abs(Qq-pi) << " " << errq << "\n";
	}
}

int main()
{
	std::cout << std::setprecision(16);

	std::cout << "Part A: plain Monte Carlo\n";

	vec a2 = {-1,-1};
	vec b2 = { 1, 1};
	func circle = [](vec x){return x[0]*x[0]+x[1]*x[1] <= 1 ? 1.0 : 0.0;};

	double Q,err;
	std::tie(Q,err) = plainmc_lcg(circle,a2,b2,100000,1);
	std::cout << "area of unit circle = " << Q << "  error estimate = " << err << "  exact = " << pi << "\n";

	vec ae = {-1,-2,-3};
	vec be = { 1, 2, 3};
	func ellipsoid = [](vec x){
		double r2 = x[0]*x[0]/(1*1) + x[1]*x[1]/(2*2) + x[2]*x[2]/(3*3);
		return r2 <= 1 ? 1.0 : 0.0;
	};

	std::tie(Q,err) = plainmc_lcg(ellipsoid,ae,be,200000,2);
	std::cout << "volume of ellipsoid = " << Q << "  error estimate = " << err << "  exact = " << 4.0/3*pi*1*2*3 << "\n\n";

	std::cout << "Part B: quasi-random sequences\n";

	vec a3 = {0,0,0};
	vec b3 = {1,1,1};
	func difficult = [](vec x){
		return 1.0/(1-std::cos(pi*x[0])*std::cos(pi*x[1])*std::cos(pi*x[2]));
	};

	double exact = 1.3932039296856768591842462603255;
	double Q_lcg,err_lcg,Q_std,err_std,Q_quasi,err_quasi;
	int N = 100000;

	std::tie(Q_lcg,err_lcg) = plainmc_lcg(difficult,a3,b3,N,3);
	std::tie(Q_std,err_std) = plainmc_std(difficult,a3,b3,N,3);
	std::tie(Q_quasi,err_quasi) = quasymc(difficult,a3,b3,N);

	std::cout << "Difficult integral, exact = " << exact << "\n";
	std::cout << "LCG:          " << Q_lcg   << "  estimated error = " << err_lcg   << "  actual error = " << std::abs(Q_lcg-exact) << "\n";
	std::cout << "std::mt19937: " << Q_std   << "  estimated error = " << err_std   << "  actual error = " << std::abs(Q_std-exact) << "\n";
	std::cout << "Halton:       " << Q_quasi << "  estimated error = " << err_quasi << "  actual error = " << std::abs(Q_quasi-exact) << "\n\n";

	std::cout << "Part C: recursive stratified sampling\n";
	double Q_plain,err_plain,Q_strat,err_strat;
	std::tie(Q_plain,err_plain) = plainmc_lcg(circle,a2,b2,50000,4);
	std::tie(Q_strat,err_strat) = stratifiedmc(circle,a2,b2,50000,64,4);
	std::cout << "circle plain:       " << Q_plain << "  error = " << err_plain << "  actual error = " << std::abs(Q_plain-pi) << "\n";
	std::cout << "circle stratified:  " << Q_strat << "  error = " << err_strat << "  actual error = " << std::abs(Q_strat-pi) << "\n";

	write_circle_errors();
	write_quasi_errors();

	return 0;
}
