#include "vector.h"

namespace pp {

void vector::print(const std::string& s) const {
    if(!s.empty()) std::cout << s << '\n';
    for(double x : data) std::cout << std::setw(12) << x << ' ';
    std::cout << '\n';
}

double vector::dot(const vector& other) const {
    if(size()!=other.size()) throw std::runtime_error("vector dot: size mismatch");
    double sum=0;
    for(int i=0;i<size();i++) sum += (*this)[i]*other[i];
    return sum;
}

double vector::norm() const {
    return std::sqrt(dot(*this));
}

vector& vector::operator+=(const vector& other) {
    if(size()!=other.size()) throw std::runtime_error("vector +=: size mismatch");
    for(int i=0;i<size();i++) (*this)[i] += other[i];
    return *this;
}

vector& vector::operator-=(const vector& other) {
    if(size()!=other.size()) throw std::runtime_error("vector -=: size mismatch");
    for(int i=0;i<size();i++) (*this)[i] -= other[i];
    return *this;
}

vector& vector::operator*=(double c) {
    for(int i=0;i<size();i++) (*this)[i] *= c;
    return *this;
}

vector& vector::operator/=(double c) {
    for(int i=0;i<size();i++) (*this)[i] /= c;
    return *this;
}

vector operator+(vector a, const vector& b) { a+=b; return a; }
vector operator-(vector a, const vector& b) { a-=b; return a; }
vector operator-(vector a) { a*=-1; return a; }
vector operator*(vector a, double c) { a*=c; return a; }
vector operator*(double c, vector a) { a*=c; return a; }
vector operator/(vector a, double c) { a/=c; return a; }

bool approx(double x, double y, double acc, double eps) {
    double d = std::abs(x-y);
    if(d <= acc) return true;
    return d <= eps*std::max(std::abs(x),std::abs(y));
}

bool approx(const vector& a, const vector& b, double acc, double eps) {
    if(a.size()!=b.size()) return false;
    for(int i=0;i<a.size();i++) if(!approx(a[i],b[i],acc,eps)) return false;
    return true;
}

} // namespace pp
