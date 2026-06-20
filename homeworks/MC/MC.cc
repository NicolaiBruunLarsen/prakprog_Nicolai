#include"MC.h"
#include<cmath>
#include<random>
#include<algorithm>


double lcg::operator()()
{
	seed = (a*seed + c)%m;
	return (seed+1.0)/(m+1.0);
}

std::tuple<double,double> plainmc(func f, vec a, vec b, int N, std::function<double()> random_double)
{
	int dim = a.size();
	double V = 1;
	for(int i=0;i<dim;i++) V *= b[i]-a[i];

	double sum = 0, sum2 = 0;
	vec x(dim);

	for(int i=0;i<N;i++){
		for(int k=0;k<dim;k++) x[k] = a[k] + random_double()*(b[k]-a[k]);
		double fx = f(x);
		sum += fx;
		sum2 += fx*fx;
	}

	double mean = sum/N;
	double sigma = std::sqrt(std::max(0.0,sum2/N - mean*mean));

	return std::make_tuple(mean*V, sigma*V/std::sqrt(N));
}

std::tuple<double,double> plainmc_lcg(func f, vec a, vec b, int N, unsigned int seed)
{
	lcg rnd(seed);
	return plainmc(f,a,b,N,[&](){return rnd();});
}

std::tuple<double,double> plainmc_std(func f, vec a, vec b, int N, unsigned int seed)
{
	std::mt19937 gen(seed);
	std::uniform_real_distribution<double> dist(0.0,1.0);
	return plainmc(f,a,b,N,[&](){return dist(gen);});
}

std::vector<int> prime_numbers(int n)
{
	std::vector<int> primes;
	int candidate = 2;
	while((int)primes.size()<n){
		bool prime = true;
		for(int p:primes){
			if(p*p>candidate) break;
			if(candidate%p==0){prime=false; break;}
		}
		if(prime) primes.push_back(candidate);
		candidate++;
	}
	return primes;
}

double corput(int n, int b)
{
	double q = 0;
	double bk = 1.0/b;
	while(n>0){
		q += (n%b)*bk;
		n /= b;
		bk /= b;
	}
	return q;
}

vec halton(int n, int dim, int startbase)
{
	std::vector<int> bases = prime_numbers(dim+startbase);
	vec x(dim);
	for(int k=0;k<dim;k++) x[k] = corput(n,bases[k+startbase]);
	return x;
}

std::tuple<double,double> quasymc(func f, vec a, vec b, int N)
{
	int dim = a.size();
	double V = 1;
	for(int i=0;i<dim;i++) V *= b[i]-a[i];

	double sum1 = 0, sum2 = 0;
	vec x(dim), y(dim);

	for(int i=1;i<=N;i++){
		vec u = halton(i,dim,0);
		vec v = halton(i,dim,dim);
		for(int k=0;k<dim;k++){
			x[k] = a[k] + u[k]*(b[k]-a[k]);
			y[k] = a[k] + v[k]*(b[k]-a[k]);
		}
		sum1 += f(x);
		sum2 += f(y);
	}

	double Q1 = sum1/N*V;
	double Q2 = sum2/N*V;

	return std::make_tuple((Q1+Q2)/2,std::abs(Q1-Q2));
}

static std::tuple<double,double> stratified(func f, vec a, vec b, int N, int nmin, std::function<double()> random_double)
{
	if(N<=2*nmin) return plainmc(f,a,b,N,random_double);

	int dim = a.size();
	int Npilot = nmin;

	vec sum_left(dim,0), sum_right(dim,0), sum2_left(dim,0), sum2_right(dim,0);
	std::vector<int> n_left(dim,0), n_right(dim,0);
	vec x(dim);

	for(int i=0;i<Npilot;i++){
		for(int k=0;k<dim;k++) x[k] = a[k] + random_double()*(b[k]-a[k]);
		double fx = f(x);
		for(int k=0;k<dim;k++){
			if(x[k] < (a[k]+b[k])/2){
				n_left[k]++;
				sum_left[k] += fx;
				sum2_left[k] += fx*fx;
			}
			else{
				n_right[k]++;
				sum_right[k] += fx;
				sum2_right[k] += fx*fx;
			}
		}
	}

	int kdiv = 0;
	double maxvar = -1;
	double var_left_best = 1, var_right_best = 1;

	for(int k=0;k<dim;k++){
		if(n_left[k]==0 || n_right[k]==0) continue;
		double mean_left = sum_left[k]/n_left[k];
		double mean_right = sum_right[k]/n_right[k];
		double var_left = std::max(0.0,sum2_left[k]/n_left[k]-mean_left*mean_left);
		double var_right = std::max(0.0,sum2_right[k]/n_right[k]-mean_right*mean_right);
		double var = var_left + var_right;
		if(var>maxvar){
			maxvar = var;
			kdiv = k;
			var_left_best = var_left;
			var_right_best = var_right;
		}
	}

	int Nleft, Nright;
	int Nrest = N - Npilot;
	if(var_left_best+var_right_best==0){
		Nleft = Nrest/2;
		Nright = Nrest-Nleft;
	}
	else{
		Nleft = int(Nrest*var_left_best/(var_left_best+var_right_best));
		if(Nleft<nmin) Nleft = nmin;
		if(Nleft>Nrest-nmin) Nleft = Nrest-nmin;
		Nright = Nrest-Nleft;
	}

	vec a_left=a, b_left=b, a_right=a, b_right=b;
	double mid = (a[kdiv]+b[kdiv])/2;
	b_left[kdiv] = mid;
	a_right[kdiv] = mid;

	double Q1,err1,Q2,err2;
	std::tie(Q1,err1) = stratified(f,a_left,b_left,Nleft,nmin,random_double);
	std::tie(Q2,err2) = stratified(f,a_right,b_right,Nright,nmin,random_double);

	return std::make_tuple(Q1+Q2,std::sqrt(err1*err1+err2*err2));
}

std::tuple<double,double> stratifiedmc(func f, vec a, vec b, int N, int nmin, unsigned int seed)
{
	lcg rnd(seed);
	return stratified(f,a,b,N,nmin,[&](){return rnd();});
}
