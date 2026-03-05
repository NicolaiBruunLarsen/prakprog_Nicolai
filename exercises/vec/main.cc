#include <iostream>
#include "vec.h"

int main() {

    vec a(1,2,3), b(4,5,6);

    a.print("a = ");
    b.print("b = "); 

    vec c = a + b;
    c.print("c = a + b = "); // c = 5 7 9
    // std::cout << c << "\n"; // { 5, 7, 9 }


    vec d = -a;
    d.print("d = -a = "); // d = -1 -2 -3

    vec e = 2.0 * a;
    e.print("e = 2.0 * a =");

    vec f = e / 3.0;
    f.print("f = e / 3.0 = ");

    // dot product wih %
    double dp = a % b;
    std::cout << "a % b = " << dp << "\n"; // a

    // cross product with ^
    vec cp = a ^ b;
    cp.print("a ^ b = "); 

    // norm 
    double n = norm(a);
    std::cout << "norm(a) = " << n << "\n";

    std::cout << a;








    return 0;
}