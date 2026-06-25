This project is an extension of the ANN homework.

Most of the code is reused from the Homerworks.

The network is

    F(x) = sum_i w_i f((x-a_i)/b_i)

with Gaussian wavelets f(x)=x exp(-x^2).

It compares four training approaches:

0. all parameters at once, (this is to compare)
1. one neuron at a time (a_i,b_i,w_i),
2. one parameter set at a time: all w_i, then all b_i, then all a_i,
3. randomly chosen parameter subsets.

The block methods use the same minimizer as before, but only on selected parameters while the rest are kept fixed.
