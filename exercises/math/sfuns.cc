#include <cmath>
#include <numbers>
#include"sfuns.h"



namespace sfuns {

constexpr double PI = std::numbers::pi;



//  for task 2
double fgamma(double x){
  if(x<0) return PI/std::sin(PI*x)/fgamma(1-x);
	if(x<9) return fgamma(x+1)/x;
	double lnfgamma=x*std::log(x+1/(12*x-1/x/10))-x+std::log(2*PI/x)/2;
	return std::exp(lnfgamma);
	}


// for task 3
#include <limits>


static double stirling_lngamma_core(double x) {
  const double t = x + 1.0 / (12.0 * x - 1.0 / (10.0 * x));
  const double lnf = x * std::log(t) - x + 0.5 * std::log(2.0 * PI / x);
  return lnf;
}


double lngamma(double x) {
  constexpr double NaN = std::numeric_limits<double>::quiet_NaN();

  if (x <= 0.0) return NaN;

  if (x < 9.0) {
    return lngamma(x + 1.0) - std::log(x);
  }

  return stirling_lngamma_core(x);
}

}
