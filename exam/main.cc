#include "ann.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

using pp::vector;

const double pi = 3.141592653589793238462643383279502884;

double exact_fun(double x){
	return std::cos(5*x-1)*std::exp(-x*x);
}

double exact_derivative(double x){
	return std::exp(-x*x)*(-5*std::sin(5*x-1)-2*x*std::cos(5*x-1));
}

double exact_second(double x){
	return std::exp(-x*x)*((4*x*x-27)*std::cos(5*x-1)+20*x*std::sin(5*x-1));
}

int main(){
	std::cout << std::setprecision(12);

	int N=25;
	vector xs(N),ys(N);
	for(int i=0;i<N;i++){
		xs[i]=-1.0 + 2.0*i/(N-1);
		ys[i]=exact_fun(xs[i]);
	}

	int n_hidden=6;

	pp::ann all(n_hidden);
	int steps_all=all.train(xs,ys,1e-5,25);

	pp::ann by_neuron(n_hidden);
	int steps_neuron=by_neuron.train_neurons(xs,ys,20,1e-5,8);

	pp::ann by_set(n_hidden);
	int steps_set=by_set.train_parameter_sets(xs,ys,20,1e-5,8);

	pp::ann random_subset(n_hidden);
	int steps_random=random_subset.train_random(xs,ys,100,4,1e-5,8,3);

	std::cout << "\"Stochastic gradient descent\" training of a network\n";
	std::cout << "hidden neurons = " << n_hidden << "\n";
	std::cout << "training points = " << N << "\n\n";

	std::cout << "method                         cost                minimizer steps\n";
	std::cout << "all parameters at once        " << all.cost(xs,ys)           << "   " << steps_all    << "\n";
	std::cout << "one neuron at a time          " << by_neuron.cost(xs,ys)     << "   " << steps_neuron << "\n";
	std::cout << "one parameter set at a time   " << by_set.cost(xs,ys)        << "   " << steps_set    << "\n";
	std::cout << "random parameter subset       " << random_subset.cost(xs,ys) << "   " << steps_random << "\n\n";



	std::ofstream costfile("out.costs.txt");
	costfile << "0 " << all.cost(xs,ys) << "\n";
	costfile << "1 " << by_neuron.cost(xs,ys) << "\n";
	costfile << "2 " << by_set.cost(xs,ys) << "\n";
	costfile << "3 " << random_subset.cost(xs,ys) << "\n";


	std::ofstream methods("out.methods.txt");
	methods << "# x exact all by_neuron by_set random_subset\n";

	for(int i=0;i<=400;i++){
		double x=-1.0+2.0*i/400;
		methods << x << " "
		        << exact_fun(x) << " "
		        << all.response(x) << " "
		        << by_neuron.response(x) << " "
		        << by_set.response(x) << " "
		        << random_subset.response(x) << "\n";
	}

	std::cout << "The values of the cost function for the different training methods are shown in fig.costs.svg.\n";
	std::cout << "I can conclude that for this run it is better to train all parameters at once, but the other methods are still able to reach a reasonable cost.\n";
	std::cout << "It could look like the single neuron training is good, but is also has the biggets amount of minimizer steps.\n";

	return 0;
}
