#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>


// int main (int argc, char* argv[]) {
// 	std::vector<double> numbers;
// 	for(int i=0;i<argc;++i){
// 		std::string arg=argv[i];
// 		if(arg=="-n" && i+1<argc)
// 			numbers.push_back(std::stod(argv[i+1]));
// 	}
    
// // exit(EXIT_SUCCESS);

// // for(auto n: numbers)
// // 	std::cout << n <<" "<< std::sin(n) <<" "<< std::cos(n) <<std::endl;



// double x;
// while( std::cin >> x ){
// 	std::cout << x <<" "<< std::sin(x) <<" "<< std::cos(x) << std::endl;
// 	}

// exit(EXIT_SUCCESS);


// }

// file streams

int main (int argc, char *argv[]) {
	std::string infile="", outfile="";
	for(int i=0;i<argc;i++){
		std::string arg=argv[i];
		if(arg=="--input" && i+1 < argc) infile=argv[i+1];
		if(arg=="--output" && i+1 < argc) outfile=argv[i+1];
	}

	// default to the provided input/output files when no args are given (not in the exercise description, but useful i think)
	if(infile.empty()) infile = "input.txt";
	if(outfile.empty()) outfile = "Out.txt";

	std::ifstream myinput(infile);
	std::ofstream myoutput(outfile);
double x;
if( myinput.is_open() && myoutput.is_open() ){
	while( myinput >> x ){
		myoutput << x <<" "<<std::sin(x)<<" "<<std::cos(x)<<std::endl;
		}
	}
else{
	std::cerr << "Error opening files: " << infile << outfile << std::endl;
	return EXIT_FAILURE;
    }
myinput.close();
myoutput.close();
exit(EXIT_SUCCESS);
}