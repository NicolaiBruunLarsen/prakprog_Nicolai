#pragma once

namespace sfuns {
  double erf(double x);   // Error function
  double sgamma(double x);   // Gamma(x) approximation (Stirling)
  double lngamma(double x);  // log(Gamma(x)) approximation (positive x only)
}
