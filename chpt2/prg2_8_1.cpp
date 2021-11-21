#include "stdafx.h"
#include "itpp/base/random.h"
using std::cin;
using std::cout;
using std::endl;
using namespace itpp;

double min,max;
vec vout1;
mat mout1;
  int main() {

  Uniform_RNG gen(0, 10);
  gen.get_setup(min,max);
  cout <<min<<endl;//prints sets value
  cout <<max<<endl;
  cout << gen() << endl; // prints a random integer
  cout << gen(10) << endl; // prints 10 random integers
  cout << gen(5,4)<<endl; //prints 5*4 random integers
  gen.sample_vector(5,vout1);
  cout <<vout1<<endl;
  gen.sample_matrix(4,3,mout1);
  cout <<mout1<<endl;
  return 0;
  }
