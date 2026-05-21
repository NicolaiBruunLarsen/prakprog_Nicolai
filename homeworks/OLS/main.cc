#include "matrix.h"
#include "QR.h"
#include "lsfit.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

pp::vector thx_time() {
    return pp::vector{1,2,3,4,6,9,10,13,15};
}

pp::vector thx_activity() {
    return pp::vector{117,100,88,72,53,29.5,25.2,15.2,11.1};
}

pp::vector thx_uncertainty() {
    return pp::vector{6,5,4,4,4,3,3,2,2};
}

bool test_qr_tall() {
    pp::matrix A(5,3);
    A[0,0]=1; A[0,1]=2; A[0,2]=3;
    A[1,0]=4; A[1,1]=5; A[1,2]=6;
    A[2,0]=7; A[2,1]=8; A[2,2]=10;
    A[3,0]=2; A[3,1]=1; A[3,2]=3;
    A[4,0]=3; A[4,1]=4; A[4,2]=5;

    pp::QR qr(A);
    pp::matrix QTQ = qr.Q.T()*qr.Q;
    pp::matrix I = pp::identity(3);
    bool orthogonal = pp::approx(QTQ, I, 1e-10, 1e-10);
    bool reconstruct = pp::approx(qr.Q*qr.R, A, 1e-10, 1e-10);

    std::cout << "Q.T()*Q == 1      " << orthogonal << "\n";
    std::cout << "Q*R == A          " << reconstruct << "\n";

    return orthogonal && reconstruct;
}

bool test_lsfit_line() {
    pp::vector x{0,1,2,3,4};
    pp::vector y{2,5,8,11,14};
    pp::vector dy{1,1,1,1,1};

    pp::Funcs fs = {
        [](double) { return 1.0; },
        [](double z) { return z; }
    };

    pp::fit_result fit = pp::lsfit(fs,x,y,dy);

    bool ok0 = pp::approx(fit.c[0], 2.0, 1e-10, 1e-10);
    bool ok1 = pp::approx(fit.c[1], 3.0, 1e-10, 1e-10);

    std::cout << "linear fit c0=2   " << ok0 << "\n";
    std::cout << "linear fit c1=3   " << ok1 << "\n";

    return ok0 && ok1;
}

void write_data_files(const pp::fit_result& fit, const pp::Funcs& fs,
                      const pp::vector& t, const pp::vector& y, const pp::vector& dy) {
    std::ofstream data("out.data.txt");
    data << "# t y dy ln_y dln_y\n";
    for (int i = 0; i < t.size(); ++i) {
        data << t[i] << " "
             << y[i] << " "
             << dy[i] << " "
             << std::log(y[i]) << " "
             << dy[i]/y[i] << "\n";
    }

    std::ofstream curve("out.fit.txt");
    curve << "# t best plus_plus plus_minus minus_plus minus_minus\n";

    double tmin = 0.0;
    double tmax = 16.0;
    int npoints = 200;

    for (int i = 0; i < npoints; ++i) {
        double z = tmin + (tmax-tmin)*i/(npoints-1);

        pp::vector cpp = fit.c;
        pp::vector cpm = fit.c;
        pp::vector cmp = fit.c;
        pp::vector cmm = fit.c;

        cpp[0] += fit.dc[0]; cpp[1] += fit.dc[1];
        cpm[0] += fit.dc[0]; cpm[1] -= fit.dc[1];
        cmp[0] -= fit.dc[0]; cmp[1] += fit.dc[1];
        cmm[0] -= fit.dc[0]; cmm[1] -= fit.dc[1];

        curve << z << " "
              << std::exp(pp::evaluate(fs, fit.c, z)) << " "
              << std::exp(pp::evaluate(fs, cpp, z)) << " "
              << std::exp(pp::evaluate(fs, cpm, z)) << " "
              << std::exp(pp::evaluate(fs, cmp, z)) << " "
              << std::exp(pp::evaluate(fs, cmm, z)) << "\n";
    }

    std::ofstream logcurve("out.logfit.txt");
    logcurve << "# t best plus_plus plus_minus minus_plus minus_minus\n";

    for (int i = 0; i < npoints; ++i) {
        double z = tmin + (tmax-tmin)*i/(npoints-1);

        pp::vector cpp = fit.c;
        pp::vector cpm = fit.c;
        pp::vector cmp = fit.c;
        pp::vector cmm = fit.c;

        cpp[0] += fit.dc[0]; cpp[1] += fit.dc[1];
        cpm[0] += fit.dc[0]; cpm[1] -= fit.dc[1];
        cmp[0] -= fit.dc[0]; cmp[1] += fit.dc[1];
        cmm[0] -= fit.dc[0]; cmm[1] -= fit.dc[1];

        logcurve << z << " "
                 << pp::evaluate(fs, fit.c, z) << " "
                 << pp::evaluate(fs, cpp, z) << " "
                 << pp::evaluate(fs, cpm, z) << " "
                 << pp::evaluate(fs, cmp, z) << " "
                 << pp::evaluate(fs, cmm, z) << "\n";
    }
}

void run_fit() {
    pp::vector t = thx_time();
    pp::vector y = thx_activity();
    pp::vector dy = thx_uncertainty();

    pp::vector logy(t.size());
    pp::vector dlogy(t.size());

    for (int i = 0; i < t.size(); ++i) {
        logy[i] = std::log(y[i]);
        dlogy[i] = dy[i]/y[i];
    }

    pp::Funcs fs = {
        [](double) { return 1.0; },
        [](double z) { return z; }
    };

    pp::fit_result fit = pp::lsfit(fs,t,logy,dlogy);

    double ln_a = fit.c[0];
    double lambda = -fit.c[1];
    double d_ln_a = fit.dc[0];
    double d_lambda = fit.dc[1];
    double a = std::exp(ln_a);
    double d_a = a*d_ln_a;

    double half_life = std::log(2.0)/lambda;
    double d_half_life = std::log(2.0)/(lambda*lambda)*d_lambda;

    double modern_half_life = 3.6319; // days, modern value for 224Ra.

    std::cout << std::setprecision(10);
    std::cout << "# Fit of ln(y)=ln(a)-lambda*t\n";
    std::cout << "ln(a)       = " << ln_a << " +/- " << d_ln_a << "\n";
    std::cout << "a           = " << a << " +/- " << d_a << "\n";
    std::cout << "lambda      = " << lambda << " +/- " << d_lambda << " 1/day\n";
    std::cout << "T_1/2       = " << half_life << " +/- " << d_half_life << " days\n";
    std::cout << "modern 224Ra T_1/2 = " << modern_half_life << " days\n";
    std::cout << "difference  = " << half_life - modern_half_life << " days\n";
    std::cout << "chi2        = " << fit.chi2 << "\n";
    std::cout << "dof         = " << fit.dof << "\n";
    std::cout << "chi2/dof    = " << fit.chi2/fit.dof << "\n";
    std::cout << "covariance matrix:\n";
    fit.cov.print();

    write_data_files(fit, fs, t, y, dy);

    std::cerr << "wrote out.data.txt, out.fit.txt, out.logfit.txt\n";
}

void print_help(const char* prog) {
    std::cout
        << "Usage:\n"
        << "  " << prog << " --test\n"
        << "  " << prog << " --fit\n"
        << "  " << prog << " --help\n";
}

int main(int argc, char** argv) {
    bool do_test = false;
    bool do_fit = false;

    if (argc == 1) {
        do_fit = true;
    }

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "--test") {
            do_test = true;
        }
        else if (arg == "--fit") {
            do_fit = true;
        }
        else if (arg == "--help" || arg == "-h") {
            print_help(argv[0]);
            return 0;
        }
        else {
            throw std::runtime_error("unknown command-line argument: " + arg);
        }
    }

    if (do_test) {
        bool ok = test_qr_tall() && test_lsfit_line();
        return ok ? 0 : 1;
    }

    if (do_fit) {
        run_fit();
        return 0;
    }

    print_help(argv[0]);
    return 0;
}
