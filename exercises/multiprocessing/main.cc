#include<thread>
#include<string>
#include<iostream>
#include<vector>
#include<functional>

struct data { int start,end; double sum; };

void harm(data& arg){
	int start=arg.start, end=arg.end; // you better work
	double sum=0;                     // with thread-local variables
	for(int i=start;i<end;i++) sum+=1.0/i;
	arg.sum=sum; // only write once into memory
}

int main(int argc, char** argv){
	int nthreads = 1, nterms = (int)1e9; /* default values */
	for(int i=0;i<argc;i++) {
		std::string arg = argv[i];
		if(arg=="-terms" && i+1<argc) nterms=(int)std::stod(argv[++i]);
		if(arg=="-threads" && i+1<argc) nthreads=std::stoi(argv[++i]);
	}

	std::cout << "terms: "<<nterms<<"\n";
	std::cout << "threads: "<<nthreads<<"\n";

	std::vector<data> params(nthreads);
	for(int i=0;i<nthreads;i++) {
params[i].start = 1 + nterms/nthreads*i;
	params[i].end = 1 + nterms/nthreads*(i+1);
}
params[params.size()-1].end = nterms+1; /* the endpoint might need adjustment */

	std::vector<std::thread> threads;
	threads.reserve(nthreads); // reserve empty slots for our threads

	for(int i=0;i<nthreads;i++) {
		threads.emplace_back(harm,std::ref(params[i]));
	}

	for(auto &thread : threads) thread.join();

	double total = 0;
	for (auto &p : params) {
		total += p.sum;
	}

	std::cout << total << "\n";
	return 0;
}