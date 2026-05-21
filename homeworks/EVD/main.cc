#include "matrix.h"
#include "EVD.h"

#include <random>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

// Random symmetric matrix for testing and timing
pp::matrix random_symmetric_matrix(int n, double scale = 1.0, unsigned seed = 1) {
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dist(-scale, scale);

    pp::matrix A(n,n);

    for (int i = 0; i < n; ++i) {
        for (int j = i; j < n; ++j) {
            double x = dist(gen);
            A[i,j] = x;
            A[j,i] = x;
        }
    }

    return A;
}

pp::matrix hydrogen_hamiltonian(double rmax, double dr) {
    int npoints = static_cast<int>(rmax/dr) - 1;

    if (npoints < 2) {
        throw std::runtime_error("too few grid points: increase rmax or decrease dr");
    }

    pp::vector r(npoints);

    for (int i = 0; i < npoints; ++i) {
        r[i] = dr*(i+1);
    }

    pp::matrix H(npoints,npoints);

    for (int i = 0; i < npoints-1; ++i) {
        H[i,i]   = -2*(-0.5/dr/dr);
        H[i,i+1] =  1*(-0.5/dr/dr);
        H[i+1,i] =  1*(-0.5/dr/dr);
    }

    H[npoints-1,npoints-1] = -2*(-0.5/dr/dr);

    for (int i = 0; i < npoints; ++i) {
        H[i,i] += -1/r[i];
    }

    return H;
}

double exact_energy(int k) {
    double n = k + 1.0;
    return -0.5/(n*n);
}

double exact_s_wave_reduced(int k, double r) {
    // First three exact reduced radial s-wave functions f(r)=r*R_{n0}(r), in Bohr units.
    if (k == 0) {
        return 2*r*std::exp(-r);
    }

    if (k == 1) {
        return (1/std::sqrt(2.0))*r*(1-r/2)*std::exp(-r/2);
    }

    if (k == 2) {
        return (2.0/(81*std::sqrt(3.0)))*r*(27 - 18*r + 2*r*r)*std::exp(-r/3);
    }

    return 0;
}

bool test_jacobi() {
    pp::matrix A = random_symmetric_matrix(6, 1.0, 1234);
    pp::EVD evd(A);

    pp::matrix D  = pp::diagonal(evd.w);
    pp::matrix VT = evd.V.T();
    pp::matrix I  = pp::identity(A.size1());

    bool test1 = pp::approx(VT*A*evd.V, D, 1e-6, 1e-6);
    bool test2 = pp::approx(evd.V*D*VT, A, 1e-6, 1e-6);
    bool test3 = pp::approx(VT*evd.V, I, 1e-6, 1e-6);
    bool test4 = pp::approx(evd.V*VT, I, 1e-6, 1e-6);

    std::cout << "V.T()*A*V == D    " << test1 << "\n";
    std::cout << "V*D*V.T() == A    " << test2 << "\n";
    std::cout << "V.T()*V == 1      " << test3 << "\n";
    std::cout << "V*V.T() == 1      " << test4 << "\n";
    std::cout << "sweeps             " << evd.sweeps << "\n";
    std::cout << "rotations          " << evd.rotations << "\n";

    return test1 && test2 && test3 && test4;
}

void write_wavefunctions(const pp::EVD& evd, double dr, const std::string& filename) {
    std::ofstream out(filename);

    int npoints = evd.V.size1();
    int nwaves = std::min(3, evd.V.size2());

    out << "# r";

    for (int k = 0; k < nwaves; ++k) {
        out << " numeric_f" << k << " exact_f" << k;
    }

    out << "\n";

    for (int i = 0; i < npoints; ++i) {
        double r = dr*(i+1);
        out << r;

        for (int k = 0; k < nwaves; ++k) {
            double sign = evd.V[0,k] < 0 ? -1.0 : 1.0;
            double numeric = sign*evd.V[i,k]/std::sqrt(dr);
            double exact = exact_s_wave_reduced(k,r);

            out << " " << numeric << " " << exact;
        }

        out << "\n";
    }
}

void run_hydrogen(double rmax, double dr) {
    pp::matrix H = hydrogen_hamiltonian(rmax, dr);
    pp::EVD evd(H);

    std::cout << "# rmax = " << rmax
              << ", dr = " << dr
              << ", npoints = " << H.size1() << "\n";

    std::cout << "# k numeric_energy exact_energy error\n";

    int nprint = std::min(5, evd.w.size());

    for (int k = 0; k < nprint; ++k) {
        double exact = exact_energy(k);

        std::cout << k << " "
                  << evd.w[k] << " "
                  << exact << " "
                  << evd.w[k] - exact << "\n";
    }

    write_wavefunctions(evd, dr, "out.wavefunctions.txt");

    std::cerr << "wrote out.wavefunctions.txt\n";
}

void convergence_dr(double rmax) {
    std::ofstream out("out.convergence.dr.txt");

    out << "# dr epsilon0 exact error npoints\n";

    double exact = exact_energy(0);

    for (double dr : {0.50, 0.40, 0.30, 0.25, 0.20, 0.15, 0.10}) {
        pp::matrix H = hydrogen_hamiltonian(rmax, dr);
        pp::EVD evd(H);

        out << dr << " "
            << evd.w[0] << " "
            << exact << " "
            << evd.w[0] - exact << " "
            << H.size1() << "\n";

        std::cerr << "dr = " << dr << " done\n";
    }
}

void convergence_rmax(double dr) {
    std::ofstream out("out.convergence.rmax.txt");

    out << "# rmax epsilon0 exact error npoints\n";

    double exact = exact_energy(0);

    for (double rmax : {3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0}) {
        pp::matrix H = hydrogen_hamiltonian(rmax, dr);
        pp::EVD evd(H);

        out << rmax << " "
            << evd.w[0] << " "
            << exact << " "
            << evd.w[0] - exact << " "
            << H.size1() << "\n";

        std::cerr << "rmax = " << rmax << " done\n";
    }
}

void timing_test(int N) {
    pp::matrix A = random_symmetric_matrix(N, 1.0, 42);

    auto start = std::chrono::high_resolution_clock::now();

    pp::EVD evd(A);

    auto stop = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> seconds = stop - start;

    std::cout << N << " "
              << seconds.count() << " "
              << evd.sweeps << " "
              << evd.rotations << "\n";
}

void print_help(const char* prog) {
    std::cout
        << "Usage:\n"
        << "  " << prog << " --test\n"
        << "  " << prog << " -rmax 10 -dr 0.3\n"
        << "  " << prog << " --convergence-dr [-rmax 10]\n"
        << "  " << prog << " --convergence-rmax [-dr 0.2]\n"
        << "  " << prog << " --time N\n";
}

int main(int argc, char** argv) {
    double rmax = 10.0;
    double dr = 0.3;

    bool do_test = false;
    bool do_conv_dr = false;
    bool do_conv_rmax = false;

    int time_N = 0;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "-rmax" && i+1 < argc) {
            rmax = std::stod(argv[++i]);
        }
        else if (arg == "-dr" && i+1 < argc) {
            dr = std::stod(argv[++i]);
        }
        else if (arg == "--test") {
            do_test = true;
        }
        else if (arg == "--convergence-dr") {
            do_conv_dr = true;
        }
        else if (arg == "--convergence-rmax") {
            do_conv_rmax = true;
        }
        else if (arg == "--time" && i+1 < argc) {
            time_N = std::stoi(argv[++i]);
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
        return test_jacobi() ? 0 : 1;
    }

    if (do_conv_dr) {
        convergence_dr(rmax);
        return 0;
    }

    if (do_conv_rmax) {
        convergence_rmax(dr);
        return 0;
    }

    if (time_N > 0) {
        timing_test(time_N);
        return 0;
    }

    run_hydrogen(rmax, dr);

    return 0;
}
