#include "itpp/itcomm.h"
#include "stdio.h"
using std::cin;
using std::cout;
using std::endl;
using namespace itpp;

int main()
{
	//dec2bin(int,int)
	bvec my_bvec;
	my_bvec = dec2bin(8,7);
	cout << my_bvec << endl;

	//dec2bin(int,bool)
	bvec my_bvec_2,my_bvec_3;
	my_bvec_2 = dec2bin(16,false);
	my_bvec_3 = dec2bin(16,true);

	cout << my_bvec_2 << endl;
    cout << my_bvec_3 << endl;

	//round(double)
	double round_x = 4.567;
	double my_x;
	my_x = itpp::round(round_x);
	cout << my_x <<endl;

	//ceil(double)
	double ceil_x = 3.4567;
	double my_ceil;
	my_ceil = ceil(ceil_x);
	cout << my_ceil << endl;
	pause(-1);
	return 0;
}
