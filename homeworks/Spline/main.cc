#include"Splines.h"
#include<iostream>
#include<fstream>
#include<cmath>
#include<cassert>
#include<string>

pp::vec make_grid(double a, double b, double dx){
	int n=(int)((b-a)/dx+0.5)+1;
	pp::vec x(n);
	for(int i=0;i<n;i++) x[i]=a+i*dx;
	return x;
}

pp::vec cos_values(const pp::vec& x){
	pp::vec y(x.size());
	for(int i=0;i<x.size();i++) y[i]=std::cos(x[i]);
	return y;
}

void test(){
	pp::vec x{0,1,2,3};
	pp::vec y{0,1,2,3};
	assert(pp::approx(pp::linterp(x,y,0.5),0.5));
	assert(pp::approx(pp::linterp(x,y,2.5),2.5));
	assert(pp::approx(pp::linterpInteg(x,y,3.0),4.5));

	pp::vec xs{1,2,3,4,5};
	pp::vec y0{1,1,1,1,1};
	pp::vec y1{1,2,3,4,5};
	pp::vec y2{1,4,9,16,25};

	pp::qspline s0(xs,y0), s1(xs,y1), s2(xs,y2);
	for(int i=0;i<s0.b.size();i++){
		assert(pp::approx(s0.b[i],0.0));
		assert(pp::approx(s0.c[i],0.0));
		assert(pp::approx(s1.b[i],1.0));
		assert(pp::approx(s1.c[i],0.0));
		assert(pp::approx(s2.b[i],2.0*(i+1)));
		assert(pp::approx(s2.c[i],1.0));
		}
	assert(pp::approx(s2.eval(2.5),2.5*2.5));
	assert(pp::approx(s2.deriv(2.5),5.0));
	assert(pp::approx(s2.integ(5.0),(125.0-1.0)/3.0));

	std::cout << "all tests passed\n";
	std::cout << "\ny=x^2 debug coefficients\n";
	std::cout << "i b[i] c[i]\n";
	for(int i=0;i<s2.b.size();i++)
		std::cout << i << " " << s2.b[i] << " " << s2.c[i] << "\n";
}

void make_data(){
	pp::vec x=make_grid(0,9,0.5);
	pp::vec y=cos_values(x);
	pp::qspline qs(x,y);

	std::ofstream datapoints("data.txt");
	datapoints << "# x y\n";
	for(int i=0;i<x.size();i++) datapoints << x[i] << " " << y[i] << "\n";

	std::ofstream lin("data_lin.txt");
	lin << "# z cos(z) sin(z) linterp linterpInteg\n";
	for(double z=x[0];z<=x[x.size()-1]+1e-12;z+=0.02)
		lin << z << " "
		    << std::cos(z) << " "
		    << std::sin(z)-std::sin(x[0]) << " "
		    << pp::linterp(x,y,z) << " "
		    << pp::linterpInteg(x,y,z) << "\n";

	std::ofstream quad("data_quad.txt");
	quad << "# z cos(z) sin(z) qspline qsplineInteg qsplineDeriv -sin(z)\n";
	for(double z=x[0];z<=x[x.size()-1]+1e-12;z+=0.02)
		quad << z << " "
		     << std::cos(z) << " "
		     << std::sin(z)-std::sin(x[0]) << " "
		     << qs.eval(z) << " "
		     << qs.integ(z) << " "
		     << qs.deriv(z) << " "
		     << -std::sin(z) << "\n";

	std::cout << "Part A: linear spline\n";
	std::cout << "xs = {";
	for(int i=0;i<x.size();i++){
	    std::cout << x[i];
	    if(i+1<x.size()) std::cout << ", ";
	}
	std::cout << "}\n";
	
	double z = 0.5;
	double exact_integral = std::sin(z)-std::sin(x[0]);
	double lin_integral = pp::linterpInteg(x,y,z);
	
	std::cout << "Integral from " << x[0] << " to " << z << "\n";
	std::cout << "Exact integral from sin(x): " << exact_integral << "\n";
	std::cout << "Calculated integral using linterpInteg: " << lin_integral << "\n\n";
	
	std::cout << "Part B: quadratic spline\n";
	double quad_integral = qs.integ(z);
	std::cout << "Integral from " << x[0] << " to " << z << "\n";
	std::cout << "Exact integral from sin(x): " << exact_integral << "\n";
	std::cout << "Calculated integral using qspline.integ: " << quad_integral << "\n\n";
	
	std::cout << "Quadratic spline coefficients for y=cos(x)\n";
	std::cout << "i x[i] b[i] c[i]\n";
	for(int i=0;i<qs.b.size();i++)
	    std::cout << i << " " << x[i] << " " << qs.b[i] << " " << qs.c[i] << "\n";
	
	std::cout << "\nData files written:\n";
	std::cout << "data.txt\n";
	std::cout << "data_lin.txt\n";
	std::cout << "data_quad.txt\n";
}

int main(int argc, char** argv){
	if(argc>1){
		std::string arg=argv[1];
		if(arg=="--test"){test(); return 0;}
		}
	make_data();
	return 0;
}
