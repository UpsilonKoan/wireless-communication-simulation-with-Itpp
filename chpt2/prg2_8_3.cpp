#include "stdafx.h"
#include "itpp/base/random.h"
using std::cin;
using std::cout;
using std::endl;
using namespace itpp;

int main() {

  Exponential_RNG gen(2.0);
  gen.setup(1.0);

  cout << gen( ) << endl; // prints a random integer
  cout << gen(3) << endl; // prints 3 random integers
  cout << gen(5,4)<<endl; //prints 5*4 random matrix integers
  
  return 0;
  }
