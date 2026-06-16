#include "matrix.h"

namespace pp {

void matrix::setid() {
    if(nrows!=ncols) throw std::runtime_error("setid: matrix is not square");
    for(int i=0;i<nrows;i++)
        for(int j=0;j<ncols;j++)
            (*this)[i,j] = (i==j ? 1.0 : 0.0);
}

matrix matrix::T() const {
    matrix B(ncols,nrows);
    for(int i=0;i<nrows;i++)
        for(int j=0;j<ncols;j++)
            B[j,i] = (*this)[i,j];
    return B;
}

void matrix::print(const std::string& s) const {
    if(!s.empty()) std::cout << s << '\n';
    for(int i=0;i<nrows;i++){
        for(int j=0;j<ncols;j++) std::cout << std::setw(12) << (*this)[i,j] << ' ';
        std::cout << '\n';
    }
}

matrix identity(int n) {
    matrix I(n,n);
    I.setid();
    return I;
}

matrix diagonal(const vector& d) {
    int n=d.size();
    matrix D(n,n);
    for(int i=0;i<n;i++) D[i,i]=d[i];
    return D;
}

matrix operator+(const matrix& A, const matrix& B) {
    if(A.size1()!=B.size1() || A.size2()!=B.size2()) throw std::runtime_error("matrix +: size mismatch");
    matrix C(A.size1(),A.size2());
    for(int i=0;i<A.size1();i++) for(int j=0;j<A.size2();j++) C[i,j]=A[i,j]+B[i,j];
    return C;
}

matrix operator-(const matrix& A, const matrix& B) {
    if(A.size1()!=B.size1() || A.size2()!=B.size2()) throw std::runtime_error("matrix -: size mismatch");
    matrix C(A.size1(),A.size2());
    for(int i=0;i<A.size1();i++) for(int j=0;j<A.size2();j++) C[i,j]=A[i,j]-B[i,j];
    return C;
}

matrix operator*(const matrix& A, double c) {
    matrix B(A.size1(),A.size2());
    for(int i=0;i<A.size1();i++) for(int j=0;j<A.size2();j++) B[i,j]=A[i,j]*c;
    return B;
}

matrix operator*(double c, const matrix& A) { return A*c; }

matrix operator/(const matrix& A, double c) { return A*(1.0/c); }

vector operator*(const matrix& A, const vector& x) {
    if(A.size2()!=x.size()) throw std::runtime_error("matrix*vector: size mismatch");
    vector y(A.size1());
    for(int i=0;i<A.size1();i++){
        double sum=0;
        for(int j=0;j<A.size2();j++) sum += A[i,j]*x[j];
        y[i]=sum;
    }
    return y;
}

matrix operator*(const matrix& A, const matrix& B) {
    if(A.size2()!=B.size1()) throw std::runtime_error("matrix*matrix: size mismatch");
    matrix C(A.size1(),B.size2());
    for(int i=0;i<A.size1();i++)
        for(int j=0;j<B.size2();j++){
            double sum=0;
            for(int k=0;k<A.size2();k++) sum += A[i,k]*B[k,j];
            C[i,j]=sum;
        }
    return C;
}

bool approx(const matrix& A, const matrix& B, double acc, double eps) {
    if(A.size1()!=B.size1() || A.size2()!=B.size2()) return false;
    for(int i=0;i<A.size1();i++)
        for(int j=0;j<A.size2();j++)
            if(!approx(A[i,j],B[i,j],acc,eps)) return false;
    return true;
}

} // namespace pp
