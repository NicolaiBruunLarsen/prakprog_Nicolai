#include "vector.h"
#include <iomanip>
#include <stdexcept>
#include <algorithm>

namespace pp {

vector::vector(int n) : data(n) {}
vector::vector(std::initializer_list<double> init) : data(init) {}

int vector::size() const {
    return static_cast<int>(data.size());
}

double& vector::operator[](int i) {
    return data.at(i);
}

const double& vector::operator[](int i) const {
    return data.at(i);
}

vector& vector::operator+=(const vector& other) {
    if (size() != other.size()) throw std::runtime_error("vector +=: size mismatch");
    for (int i = 0; i < size(); ++i) data[i] += other[i];
    return *this;
}

vector& vector::operator-=(const vector& other) {
    if (size() != other.size()) throw std::runtime_error("vector -=: size mismatch");
    for (int i = 0; i < size(); ++i) data[i] -= other[i];
    return *this;
}

vector& vector::operator*=(double c) {
    for (double& x : data) x *= c;
    return *this;
}

vector& vector::operator/=(double c) {
    for (double& x : data) x /= c;
    return *this;
}

double vector::dot(const vector& other) const {
    if (size() != other.size()) throw std::runtime_error("dot: size mismatch");
    double sum = 0;
    for (int i = 0; i < size(); ++i) sum += data[i]*other[i];
    return sum;
}

double vector::norm() const {
    return std::sqrt(dot(*this));
}

void vector::print(const std::string& s) const {
    if (!s.empty()) std::cout << s << '\n';
    for (int i = 0; i < size(); ++i) {
        std::cout << std::setw(16) << data[i];
    }
    std::cout << '\n';
}

vector operator+(vector a, const vector& b) { return a += b; }
vector operator-(vector a, const vector& b) { return a -= b; }
vector operator-(vector a) { return a *= -1; }
vector operator*(vector a, double c) { return a *= c; }
vector operator*(double c, vector a) { return a *= c; }
vector operator/(vector a, double c) { return a /= c; }

std::ostream& operator<<(std::ostream& os, const vector& v) {
    os << "{";
    for (int i = 0; i < v.size(); ++i) {
        os << v[i];
        if (i + 1 < v.size()) os << ", ";
    }
    os << "}";
    return os;
}

double dot(const vector& a, const vector& b) {
    return a.dot(b);
}

double norm(const vector& v) {
    return v.norm();
}

bool approx(double x, double y, double acc, double eps) {
    double diff = std::abs(x-y);
    if (diff <= acc) return true;
    return diff <= eps*std::max(std::abs(x), std::abs(y));
}

bool approx(const vector& a, const vector& b, double acc, double eps) {
    if (a.size() != b.size()) return false;
    for (int i = 0; i < a.size(); ++i) {
        if (!approx(a[i], b[i], acc, eps)) return false;
    }
    return true;
}

} // namespace pp
