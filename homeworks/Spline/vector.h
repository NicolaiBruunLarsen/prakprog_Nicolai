#pragma once
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>
#include <stdexcept>

namespace pp {

struct vector {
    std::vector<double> data;

    vector() = default;
    explicit vector(int n) : data(n) {}
    vector(int n, double value) : data(n,value) {}
    vector(std::initializer_list<double> init) : data(init) {}

    int size() const { return static_cast<int>(data.size()); }
    void resize(int n) { data.resize(n); }

    double& operator[](int i) { return data.at(i); }
    const double& operator[](int i) const { return data.at(i); }

    void print(const std::string& s="") const;
    double dot(const vector& other) const;
    double norm() const;

    vector& operator+=(const vector& other);
    vector& operator-=(const vector& other);
    vector& operator*=(double c);
    vector& operator/=(double c);
};

vector operator+(vector a, const vector& b);
vector operator-(vector a, const vector& b);
vector operator-(vector a);
vector operator*(vector a, double c);
vector operator*(double c, vector a);
vector operator/(vector a, double c);

bool approx(double x, double y, double acc=1e-9, double eps=1e-9);
bool approx(const vector& a, const vector& b, double acc=1e-9, double eps=1e-9);

} // namespace pp
