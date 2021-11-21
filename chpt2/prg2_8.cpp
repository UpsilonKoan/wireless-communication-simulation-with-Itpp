#include "stdafx.h"
#include "itpp/base/random.h"
using std::cin;
using std::cout;
using std::endl;
using namespace itpp;

int min,max;
  int main() {

  I_Uniform_RNG gen(0, 10);
  gen.get_setup(min,max);
  cout <<min<<endl;//prints sets value
  cout <<max<<endl;
  cout << gen() << endl; // prints a random integer
  cout << gen(10) << endl; // prints 10 random integers
  cout << gen(5,4)<<endl; //prints 5*4 random integers
  return 0;
  }
