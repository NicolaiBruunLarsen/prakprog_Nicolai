#pragma once
#include <vector>
#include <string>
#include <iostream>
#include <cmath>

namespace pp {

struct vector {
    std::vector<double> data;

    vector() = default;
    explicit vector(int n);
    vector(std::initializer_list<double> init);

    int size() const;

    double& operator[](int i);
    const double& operator[](int i) const;

    vector& operator+=(const vector& other);
    vector& operator-=(const vector& other);
    vector& operator*=(double c);
    vector& operator/=(double c);

    double dot(const vector& other) const;
    double norm() const;

    void print(const std::string& s = "") const;
};

vector operator+(vector a, const vector& b);
vector operator-(vector a, const vector& b);
vector operator-(vector a);
vector operator*(vector a, double c);
vector operator*(double c, vector a);
vector operator/(vector a, double c);

std::ostream& operator<<(std::ostream& os, const vector& v);

double dot(const vector& a, const vector& b);
double norm(const vector& v);

bool approx(double x, double y, double acc = 1e-9, double eps = 1e-9);
bool approx(const vector& a, const vector& b, double acc = 1e-9, double eps = 1e-9);

} // namespace pp
