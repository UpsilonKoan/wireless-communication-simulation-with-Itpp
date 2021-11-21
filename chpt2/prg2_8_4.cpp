#include "stdafx.h"
#include "itpp/base/random.h"
using std::cin;
using std::cout;
using std::endl;
using namespace itpp;

double mean1,var1,cho1;
int main() {

  AR1_Normal_RNG  gen(0,2.0,0.1);

  gen.get_setup(mean1,var1,cho1);
  cout<<"mean1="<<mean1<<endl;
  cout<<"var1="<<var1<<endl;
  cout<<"cho1="<<cho1<<endl;

  cout << gen( ) << endl; // prints a random integer
  cout << gen(3) << endl; // prints 10 random integers
  cout << gen(5,4)<<endl; //prints 5*4 random integers
  
  return 0;
  }
