#include"QR.h"
#include<stdexcept>
#include<cmath>

namespace pp{

QR::QR(const matrix& A) : Q(A), R(A.size2(),A.size2()){
	int n=A.size1();
	int m=A.size2();
	if(n<m) throw std::runtime_error("QR: expected n>=m");

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

vector QR::qt_mul(const vector& b) const{
	if(Q.size1()!=b.size()) throw std::runtime_error("QR qt_mul: size mismatch");
	vector c(Q.size2());
	for(int i=0;i<Q.size2();i++) c[i]=dot(Q.col(i),b);
	return c;
}

vector QR::solve(const vector& b) const{
	if(Q.size1()!=b.size()) throw std::runtime_error("QR solve: size mismatch");
	int n=R.size1();
	vector c=qt_mul(b);
	vector x(n);

	for(int i=n-1;i>=0;i--){
		double sum=c[i];
		for(int k=i+1;k<n;k++) sum-=R(i,k)*x[k];
		x[i]=sum/R(i,i);
	}
	return x;
}

double QR::det() const{
	if(R.size1()!=R.size2()) throw std::runtime_error("QR det: R not square");
	double d=1.0;
	for(int i=0;i<R.size1();i++) d*=R(i,i);
	return d;
}

matrix QR::inverse() const{
	if(Q.size1()!=Q.size2()) throw std::runtime_error("QR inverse: matrix not square");
	int n=Q.size1();
	matrix B(n,n);
	for(int j=0;j<n;j++){
		vector e(n);
		e[j]=1.0;
		B.set_col(j,solve(e));
	}
	return B;
}

} // namespace pp
