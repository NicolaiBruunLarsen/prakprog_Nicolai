#include "linalg.hpp"
#include <iostream>

int main(int argc, char** argv) {
    int N = 5;
    if (argc == 3 && std::string(argv[1]) == "-size") {
        N = std::stoi(argv[2]);
    }


    
    {
        std::cout << "-----Task A-----\n";
        std::cout << "-----Part 2-----\n";

        int n = 6, m = 4;
        pp::matrix A = pp::random_matrix(n, m, 1.0, 11);

        A.print("Tall matrix A:");

        pp::qr decomp(A);

        pp::matrix QTQ = decomp.Q.transpose() * decomp.Q;
        pp::matrix I = pp::matrix::id(m);
        pp::matrix QR = decomp.Q * decomp.R;

        decomp.Q.print("Q:");

        decomp.R.print("R:");

        QR.print("Q * R = A:");

        QTQ.print("Q^T * Q = I:");

        std::cout << "Check list:\n";
        std::cout << "R upper triangular: " << pp::is_upper_triangular(decomp.R) << "\n";
        std::cout << "Q^T Q = I:          " << pp::approx(QTQ, I, 1e-8, 1e-8) << "\n";
        std::cout << "Q R = A:            " << pp::approx(QR, A, 1e-8, 1e-8) << "\n\n";
    }

    //  Solve Ax=b
    {
        std::cout << "-----Part 3-----\n";
        std::cout << "Solve Ax=b\n";

        int n = 5;
        pp::matrix A = pp::random_matrix(n, n, 1.0, 21);
        pp::vector b = pp::random_vector(n, 1.0, 22);

        pp::qr decomp(A);
        pp::vector x = decomp.solve(b);
        pp::vector Ax = A * x;

        b.print("b:");
        x.print("x = A^{-1} b:");

        Ax.print("A * x: ");

        std::cout << "approx(A*x, b) = " << pp::approx(Ax, b) << "\n";
        std::cout << "det(A) from R:    " << decomp.det() << "\n\n";
    }

    // Part B: Inverse
    {
        std::cout << "-----Task B-----\n";
        std::cout << "Inverse of A\n";

        int n = 5;
        pp::matrix A = pp::random_matrix(n, n, 1.0, 31);

        A.print("Matrix A:");

        pp::qr decomp(A);
        pp::matrix A_inv = decomp.inverse();
        pp::matrix AB = A * A_inv;
        pp::matrix I = pp::matrix::id(n);

        A_inv.print("Inverse A^{-1}:");
        AB.print("A * A^{-1}:");

        std::cout << "Inverse checks:\n";
        std::cout << "approx(A*A^{-1}, I) = " << pp::approx(AB, I) << "\n\n";
    }

    // Part C:Timing 
    {
        pp::matrix A = pp::random_matrix(N, N, 1.0, 123);
        pp::qr decomp(A);
    }

    return 0;
}