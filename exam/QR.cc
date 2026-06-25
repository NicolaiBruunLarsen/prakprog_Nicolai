#include "QR.h"
#include <stdexcept>
#include <cmath>

namespace pp{

QR::QR(const matrix& A) : Q(A), R(A.size2(),A.size2()){
	int n=A.size1(), m=A.size2();
	if(n<m) throw std::runtime_error("QR: n<m");

	for(int i=0;i<m;i++){
		vector qi=Q.col(i);
		double rii=qi.norm();
		if(approx(rii,0.0)) throw std::runtime_error("QR: singular matrix");
		R(i,i)=rii;
		qi/=rii;
		Q.set_col(i,qi);

		for(int j=i+1;j<m;j++){
			vector qj=Q.col(j);
			double rij=dot(qi,qj);
			R(i,j)=rij;
			qj-=qi*rij;
			Q.set_col(j,qj);
		}
	}
}

vector QR::solve(const vector& b) const{
	if(b.size()!=Q.size1()) throw std::runtime_error("QR solve: bad size");
	int n=R.size1();
	vector c(n),x(n);
	for(int i=0;i<n;i++) c[i]=dot(Q.col(i),b);

	for(int i=n-1;i>=0;i--){
		double s=c[i];
		for(int k=i+1;k<n;k++) s-=R(i,k)*x[k];
		x[i]=s/R(i,i);
	}
	return x;
}

matrix QR::inverse() const{
	if(Q.size1()!=Q.size2()) throw std::runtime_error("QR inverse: not square");
	int n=Q.size1();
	matrix B(n,n);
	for(int i=0;i<n;i++){
		vector e(n);
		e[i]=1;
		B.set_col(i,solve(e));
	}
	return B;
}

} // namespace pp
