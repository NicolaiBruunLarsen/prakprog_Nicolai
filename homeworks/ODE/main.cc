#include "ODE.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <cassert>
#include <string>

void write_solution(const std::string& filename,
                    const std::vector<double>& xs,
                    const std::vector<pp::vector>& ys)
{
    std::ofstream out(filename);
    for(int i=0;i<(int)xs.size();i++){
        out << xs[i];
        for(int k=0;k<ys[i].size();k++) out << "\t" << ys[i][k];
        out << "\n";
    }
}

void write_oscillator(const std::string& filename,
                      const std::vector<double>& xs,
                      const std::vector<pp::vector>& ys)
{
    std::ofstream out(filename);
    out << "# x u uprime sin(x) cos(x)\n";
    for(int i=0;i<(int)xs.size();i++){
        double x=xs[i];
        out << x << "\t" << ys[i][0] << "\t" << ys[i][1]
            << "\t" << std::sin(x) << "\t" << std::cos(x) << "\n";
    }
}

void write_orbit(const std::string& filename,
                 const std::vector<double>& phis,
                 const std::vector<pp::vector>& ys)
{
    std::ofstream out(filename);
    out << "# phi u uprime x y\n";
    for(int i=0;i<(int)phis.size();i++){
        double phi=phis[i];
        double u=ys[i][0];
        double r=1/u;
        out << phi << "\t" << u << "\t" << ys[i][1]
            << "\t" << r*std::cos(phi) << "\t" << r*std::sin(phi) << "\n";
    }
}

void part_a(){
    std::cout << "Part a: embedded Runge-Kutta stepper and adaptive driver\n";

    const double pi = std::acos(-1.0);

    auto oscillator = [](double x, pp::vector y){
        (void)x;
        pp::vector dydx(2);
        dydx[0] = y[1];
        dydx[1] = -y[0];
        return dydx;
    };

    pp::vector yinit{0,1};
    auto [xs,ys] = pp::driver(oscillator,0,2*pi,yinit,0.125,1e-5,1e-5);
    write_oscillator("out.oscillator.txt",xs,ys);

    std::cout << "u''=-u test at x=2*pi:\n";
    std::cout << "u(2*pi)  = " << ys.back()[0] << " exact " << 0 << "\n";
    std::cout << "u'(2*pi) = " << ys.back()[1] << " exact " << 1 << "\n";
    assert(pp::approx(ys.back()[0],0,1e-3,1e-3));
    assert(pp::approx(ys.back()[1],1,1e-3,1e-3));

    auto pendulum = [](double t, pp::vector y){
        (void)t;
        double b=0.25, c=5.0;
        pp::vector dydt(2);
        dydt[0] = y[1];
        dydt[1] = -b*y[1] - c*std::sin(y[0]);
        return dydt;
    };

    pp::vector pend_init{pi-0.1,0};
    auto [tp,yp] = pp::driver(pendulum,0,10,pend_init,0.05,1e-4,1e-4);
    write_solution("out.pendulum.txt",tp,yp);

    auto lotka = [](double t, pp::vector y){
        (void)t;
        double alpha=1.5, beta=1.0, gamma=3.0, delta=1.0;
        pp::vector dydt(2);
        dydt[0] = alpha*y[0] - beta*y[0]*y[1];
        dydt[1] = delta*y[0]*y[1] - gamma*y[1];
        return dydt;
    };

    pp::vector lv_init{10,5};
    auto [tlv,ylv] = pp::driver(lotka,0,15,lv_init,0.05,1e-4,1e-4);
    write_solution("out.lotka.txt",tlv,ylv);

    std::cout << "wrote out.oscillator.txt, out.pendulum.txt, out.lotka.txt\n";
}

void part_b(){
    std::cout << "\nPart b: relativistic precession of planetary orbit\n";

    const double pi = std::acos(-1.0);
    double a=0;
    double b=12*pi;

    auto orbit = [](double epsilon){
        return [epsilon](double phi, pp::vector y){
            (void)phi;
            pp::vector dydphi(2);
            dydphi[0] = y[1];
            dydphi[1] = 1 - y[0] + epsilon*y[0]*y[0];
            return dydphi;
        };
    };

    pp::vector y_circular{1,0};
    auto [pc,yc] = pp::driver(orbit(0.0),a,b,y_circular,0.05,1e-5,1e-5);
    write_orbit("out.orbit.circular.txt",pc,yc);

    pp::vector y_ellipse{1,-0.5};
    auto [pe,ye] = pp::driver(orbit(0.0),a,b,y_ellipse,0.05,1e-5,1e-5);
    write_orbit("out.orbit.ellipse.txt",pe,ye);

    pp::vector y_relativistic{1,-0.5};
    auto [pr,yr] = pp::driver(orbit(0.01),a,b,y_relativistic,0.05,1e-5,1e-5);
    write_orbit("out.orbit.relativistic.txt",pr,yr);

    std::cout << "wrote out.orbit.circular.txt, out.orbit.ellipse.txt, out.orbit.relativistic.txt\n";
}

int main(){
    part_a();
    part_b();
    std::cout << "\nall done\n";
    return 0;
}
