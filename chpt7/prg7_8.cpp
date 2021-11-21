#include "stdafx.h"
#include <itpp/itbase.h>
#include <itpp/itcomm.h>
using namespace itpp;

//These lines are needed for use of cout and endl
using std::cin;
using std::cout;
using std::endl;

int main()
{
cout<<"**************spreading code************"<<endl;
	smat SpreadingCode(8,8);
	SpreadingCode = wcdma_spreading_codes(8);
	cout<<"wcdmaÀ©ÆµÂë"<<endl;
cout<<SpreadingCode<<endl;
return 0;
}
