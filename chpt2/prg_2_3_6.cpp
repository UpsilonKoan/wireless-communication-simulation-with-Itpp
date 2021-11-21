#include "itpp/itbase.h"
#include "stdafx.h"

using std::cin;
using std::cout;
using std::endl;
using namespace itpp;

double f(const double x)
  {
    return x*log(x);
  }

  int main()
  {
    double res = quad( f, 1.5, 3.5);
    cout << "res = " << res << endl;

    return 0;
  }
