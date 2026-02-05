#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <numbers>

#include "sfuns.h"

int main() {
  using cd = std::complex<double>;
  const cd I(0.0, 1.0);

  constexpr double PI = std::numbers::pi;
  constexpr double E  = std::numbers::e;

  std::cout << std::setprecision(15);

  // Task 1
  std::cout << "== Task 1 ==\n";
  std::cout << "sqrt(2)      \t = " << std::sqrt(2.0) << "\n";
  std::cout << "2^(1/5)      \t = " << std::pow(2.0, 1.0 / 5.0) << "\n";
  std::cout << "e^pi         \t = " << std::exp(PI) << "\n";
  std::cout << "e^i          \t = " << std::exp(I) << "\n";
  std::cout << "pi^e         \t = " << std::pow(PI, E) << "\n";
  std::cout << "pi^i         \t = " << std::pow(PI, I) << "\n";
  std::cout << "e^i (pow)    \t = " << std::pow(E, I) << "\n";
  std::cout << "i^i          \t = " << std::pow(I, I) << "\n";
  std::cout << "log(i)       \t = " << std::log(I) << "\n\n";

  // Task 2
  std::cout << "== Task 2: Gamma via Stirling approximation ==\n";
  for (int n = 1; n <= 10; ++n) {
    const double x = static_cast<double>(n);
    std::cout << "Gamma(" << n << ") \t = " << sfuns::fgamma(x) << "\n";
  }
  std::cout << "\n";

  // Task 3
  std::cout << "== Task 3: log-Gamma (positive x only) ==\n";
  for (int n = 1; n <= 10; ++n) {
    const double x = static_cast<double>(n);
    std::cout << "lngamma(" << n << ") \t = " << sfuns::lngamma(x) << "\n";
  }

  return 0;
}
