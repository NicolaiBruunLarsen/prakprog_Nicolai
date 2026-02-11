#include<iostream>
#include<limits>
#include<iomanip>

#include <cmath>      // std::abs
#include <algorithm>  // std::max

bool approx(double a, double b, double acc = 1e-9, double eps = 1e-9) {
    double diff = std::abs(a - b);

    if (diff <= acc) return true; // absolute tolerance

    double max_ab = std::max(std::abs(a), std::abs(b));

    if (max_ab == 0.0) return false;

    return diff <= eps * max_ab;
}


int main() {

    float       f=1.0f; while((float)      (1.0f+f) != 1.0f){f/=2.0f;} f*=2.0f;
    double      d=1.0d; while((double)     (1.0d+d) != 1.0d){d/=2.0d;} d*=2.0d;
    long double l=1.0L; while((long double)(1.0L+l) != 1.0L){l/=2.0L;} l*=2.0L;
    std::printf("      float eps=%g\n",f);
    std::printf("     double eps=%g\n",d);
    std::printf("long double eps=%Lg\n",l);


    // using limit header
        
    std::cout << "Using limts" << "\n";
    std::cout << std::numeric_limits<float>::epsilon() << "\n";
    std::cout << std::numeric_limits<double>::epsilon() << "\n";
    std::cout << std::numeric_limits<long double>::epsilon() << "\n";


    // 2. Non-commutativity of addition
    std::cout << "2. Non-commutativity of addition" << "\n";

    double tiny = std::numeric_limits<double>::epsilon()/2.0;

    double a=1+tiny+tiny; // should be 1+epsilon, that is, larger than 1, no?
    double b=tiny+tiny+1; // should be the same mathematically, no?

    std::cout << "a=" << a << "\n";
    std::cout << "b=" << b << "\n";
    std::cout << "a==b? " << (a==b ? "yes":"no") << "\n";


    // use iomanip
    std::cout << "Using iomanip" << "\n";
    std::cout << std::fixed << std::setprecision(17);
    std::cout << "       tiny=" << tiny << "\n";
    std::cout << "1+tiny+tiny=" << a << "\n";
    std::cout << "tiny+tiny+1=" << b << "\n";

    // 3. Comparing doubles: introduction

    double d1 = 0.1+0.1+0.1+0.1+0.1+0.1+0.1+0.1;
    double d2 = 8*0.1;

    std::cout << "d1==d2? " << (d1==d2 ? "true":"false") << "\n"; 

    std::cout << std::fixed << std::setprecision(17);
    std::cout << "d1=" << d1 << "\n";
    std::cout << "d2=" << d2 << "\n";

    // Comparing doubles: the task

    std::cout << "Comparing doubles task" << "\n";

    std::cout << "d1 = " << d1 << "\n";
    std::cout << "d2 = " << d2 << "\n";

    std::cout << "d1 == d2? " << (d1 == d2 ? "yes" : "no") << "\n";
    std::cout << "approx(d1,d2)? " << (approx(d1, d2) ? "yes" : "no") << "\n";


    return 0;
}