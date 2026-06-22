#include "ODE.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace pp {

std::tuple<vector,vector> rkstep12(
    std::function<vector(double,vector)> f,
    double x,
    vector y,
    double h)
{
    vector k0 = f(x,y);                  // Euler, lower order
    vector k1 = f(x+h/2, y+k0*(h/2));    // midpoint, higher order
    vector yh = y + k1*h;                // estimate y(x+h)
    vector dy = (k1-k0)*h;               // error estimate
    return std::make_tuple(yh,dy);
}

std::tuple<std::vector<double>,std::vector<vector>> driver(
    std::function<vector(double,vector)> f,
    double a,
    double b,
    vector yinit,
    double h,
    double acc,
    double eps)
{
    double x = a;
    vector y = yinit;

    std::vector<double> xlist;
    std::vector<vector> ylist;
    xlist.push_back(x);
    ylist.push_back(y);

    if(b<a) h = -std::abs(h);
    else    h =  std::abs(h);

    int steps = 0;
    int max_steps = 1000000;

    while(true){
        if((b>=a && x>=b) || (b<a && x<=b))
            return std::make_tuple(xlist,ylist);

        if((b>=a && x+h>b) || (b<a && x+h<b)) h = b-x;

        auto [yh,dy] = rkstep12(f,x,y,h);

        double tol = (acc + eps*yh.norm())*std::sqrt(std::abs(h/(b-a)));
        double err = dy.norm();

        if(err <= tol){
            x += h;
            y = yh;
            xlist.push_back(x);
            ylist.push_back(y);
        }

        if(err > 0){
            double factor = std::pow(tol/err,0.25)*0.95;
            factor = std::min(factor,2.0);
            factor = std::max(factor,0.1);
            h *= factor;
        }
        else h *= 2;

        steps++;
        if(steps > max_steps) throw std::runtime_error("driver: too many steps");
        if(std::abs(h) < 1e-14) throw std::runtime_error("driver: step-size underflow");
    }
}

} // namespace pp
