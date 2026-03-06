#include <iostream>
#include <iomanip>
#include <cmath>
#include "sfuns.h"

int main() {
    std::cout << std::scientific;

    // Dataset 0: erf(x)
    for (double x = -3; x <= 3; x += 1.0 / 8) {
        std::cout << x << " " << sfuns::erf(x) << "\n";
    }

    std::cout << "\n\n";

    // Dataset 1: gamma(x)
    double dx = 1.0 / 64;
    for (double x = -5 + dx / 2; x <= 8; x += dx) {
        std::cout << x << " " << sfuns::sgamma(x) << "\n";
    }

    std::cout << "\n\n";

    // Dataset 2: ln(gamma(x))
    for (double x = dx; x <= 10; x += dx) {
        std::cout << x << " " << sfuns::lngamma(x) << "\n";
    }

    std::cout << "\n\n";

    // Dataset 3: factorial points, since Gamma(n+1)=n!
    double f = 1;
    std::cout << 1 << " " << f << "\n"; // Gamma(1)=0!=1
    for (int n = 1; n <= 10; n++) {
        f *= n; // now f = n!
        std::cout << n + 1 << " " << f << "\n";
    }

    return 0;
}