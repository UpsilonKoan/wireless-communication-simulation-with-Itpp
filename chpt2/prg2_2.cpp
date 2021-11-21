#include "stdafx.h"
#include "itpp/itbase.h"

using std::cin;
using std::cout;
using std::endl;
using namespace itpp;

int main(int argc, char* argv[])
{
  vec a = randn(2);
  vec b = randn(3);
  vec c = randn(4);
  Array<vec> my_array(3);
  my_array(0) = a;
  my_array(1) = b;
  my_array(2) = c;
  cout<<endl;
  cout<<"a="<<a<<endl;
  cout<<"b="<<b<<endl;
  cout<<"c="<<c<<endl;
  cout<<"my_array="<<my_array<<endl;
 

  return 0;
}
