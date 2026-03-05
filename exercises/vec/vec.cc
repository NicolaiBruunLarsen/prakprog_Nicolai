#include "vec.h"
#include <cmath>        // std::abs, std::sqrt

vec& vec::operator+=(const vec& v) {
    x += v.x;
    y += v.y;
    z += v.z;
    return *this;
}

vec& vec::operator-=(const vec& v) {
    x -= v.x;
    y -= v.y;
    z -= v.z;
    return *this;
}

vec& vec::operator*=(double a) {
    x *= a;
    y *= a;
    z *= a;
    return *this;
}

vec& vec::operator/=(double a) {
    // Optional safety check. You can remove it if you prefer silent inf/nan.
    if (a == 0.0) throw std::runtime_error("vec::operator/=: division by zero");
    x /= a;
    y /= a;
    z /= a;
    return *this;
}

void vec::print(const std::string& s) const {
    std::cout << s << "[" << x << " " << y << " " << z << "]\n";
}


std::ostream& operator<<(std::ostream& os, const vec& v) {
    os << "{ " << v.x << ", " << v.y << ", " << v.z << " }";
    return os;
}


vec operator-(const vec& v) {
    return vec(-v.x, -v.y, -v.z);
}

vec operator-(const vec& a, const vec& b) {
    vec r = a;
    r -= b;
    return r;
}

// this is the old form 
vec operator+(const vec& a, const vec& b) {
    vec r = a;
    r += b;
    return r;
}

vec operator*(const vec& v, double a) {
    vec r = v;
    r *= a;
    return r;
}

vec operator*(double a, const vec& v) {
    return v * a; // reuse
}

vec operator/(const vec& v, double a) {
    vec r = v;
    r /= a;
    return r;
}


// dot product
double operator%(const vec& a, const vec& b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

// cross product
vec operator^(const vec& a, const vec& b) {
    return vec(
        a.y*b.z - a.z*b.y,
        a.z*b.x - a.x*b.z,
        a.x*b.y - a.y*b.x
    );
}

// norm
double norm(const vec& v) {
    return std::sqrt(v % v);
}


static bool approx_double(double a, double b, double acc, double eps) {
    // True if |a-b| <= acc OR |a-b| <= max(|a|,|b|)*eps
    double diff = std::abs(a - b);
    if (diff <= acc) return true;
    double scale = std::max(std::abs(a), std::abs(b));
    return diff <= scale * eps;
}

bool approx(const vec& a, const vec& b, double acc, double eps) {
    if (!approx_double(a.x, b.x, acc, eps)) return false;
    if (!approx_double(a.y, b.y, acc, eps)) return false;
    if (!approx_double(a.z, b.z, acc, eps)) return false;
    return true;
}