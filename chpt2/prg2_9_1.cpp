#include "stdafx.h"
#include <itpp/stat/histogram.h>
#include <itpp/itbase.h>

using std::cin;
using std::cout;
using std::endl;
using namespace itpp;

int main() {
  
	Histogram<double>  hist(0, 99, 100);
	//hist.setup(0, 99, 100);
 
    hist.update(randn(100));

  // Get position of bin number 5
  double bin5_center = hist.get_bin_center(5);
  cout<<"bin5_center="<<bin5_center<<endl;
  // Get corresponding bin counter
  int bin5_counter = hist.get_bin(5);
  cout<<"bin5_counter="<<bin5_counter<<endl;

  // Get bin 5 left boundary:
  double bin5_left = hist.get_bin_left(5);
  cout<<"bin5_left="<<bin5_left<<endl;

  // compute PDF & CDF of experimental data
  vec my_data_pdf = hist.get_pdf();
  vec my_data_cdf = hist.get_cdf();
  cout<<"my_data_pdf="<<my_data_pdf<<endl;
  cout<<"my_data_cdf="<<my_data_cdf<<endl;

  return 0;
  }
