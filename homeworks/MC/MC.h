#pragma once
#include<functional>
#include<vector>
#include<tuple>

using vec = std::vector<double>;
using func = std::function<double(vec)>;

struct lcg{
	unsigned int seed;
	unsigned int a = 1664525;
	unsigned int c = 1013904223;
	unsigned long long m = 4294967296ULL;
	lcg(unsigned int s=1):seed(s){}
	double operator()();
};

std::tuple<double,double> plainmc(func f, vec a, vec b, int N, std::function<double()> random_double);
std::tuple<double,double> plainmc_lcg(func f, vec a, vec b, int N, unsigned int seed=1);
std::tuple<double,double> plainmc_std(func f, vec a, vec b, int N, unsigned int seed=1);

double corput(int n, int b);
vec halton(int n, int dim, int startbase=0);
std::tuple<double,double> quasymc(func f, vec a, vec b, int N);

std::tuple<double,double> stratifiedmc(func f, vec a, vec b, int N, int nmin=64, unsigned int seed=1);
