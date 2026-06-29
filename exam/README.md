This project is an extension of the ANN homework.

Most of the code is reused from the Homerworks.

The network is

    F(x) = sum_i w_i f((x-a_i)/b_i)

where a_i, b_i, and w_i are the parameters of neuron number i. In my implementation the parameters are stored in one vector

    p = {a0, a1, ..., b0, b1, ..., w0, w1, ...}

and the activation function is the Gaussian wavelet f(x)=x exp(-x^2).

The network is trained by minimizing the cost/loss function

    C(p) = sum_k (F(x_k)-y_k)^2

where (x_k,y_k) are the tabulated training points. The minimization is done using my own minimization from the previous minimization homework.

It compares four training approaches:

0. all parameters at once, (this is to compare)
1. one neuron at a time (a_i,b_i,w_i),
2. one parameter set at a time: all w_i, then all b_i, then all a_i,
3. randomly chosen parameter subsets.

the class therfore has these for methods for training the network:

int train(vector xs, vector ys, double acc, int max_iter); 
int train_neurons(vector xs, vector ys, int sweeps, double acc, int max_iter); 
int train_parameter_sets(vector xs, vector ys, int sweeps, double acc, int max_iter); 
int train_random(vector xs, vector ys, int sweeps, int nsubset, double acc, int max_iter, unsigned seed);

After implementing the methods, I tested them on the function

g(x) = cos(5x-1) exp(-x^2)

sampled on the interval [-1,1].

The plots are:

fig.methods.svg  compares the exact function with the four training methods
fig.costs.svg    compares the final cost/loss for each method

The cost plot shows which method gives the smallest final loss. The number of minimizer steps is also included in Out.txt because the different methods do not necessarily use the same amount of training work. A method with a lower cost is not automatically better if it used many more minimizer steps.

I have completed all the tasks therefore i think this could be a 10/10.