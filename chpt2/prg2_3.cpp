#include "stdafx.h"
#include "itpp/itbase.h"

using std::cin;
using std::cout;
using std::endl;
using namespace itpp;

int main(int argc, char* argv[])
{
	//Declarations of scalars and vectors: 
	double x=0.9;
	double y=2;
 
	cout<<"x="<<x<<endl;
	cout<<"erf(x)="<<erf(x)<<endl;
	cout<<"asinh(x)="<<asinh(x)<<endl;

    ivec  i = "4:16";
	vec  e;e = "3:2.5:13";
	cout<<"i="<<i<<endl;
	cout<<"e="<<e<<endl;

	cout<<"pow (e, 2.1)="<<pow (e, 2.1)<<endl;
	cout<<"int2bits (i)="<<int2bits (i)<<endl;
    cout<<"levels2bits  (i)="<<levels2bits(i)<<endl;

	bvec bin_list="1 0 1 1 1 0 0 0 0";
	
	cout<<"bin_list="<<bin_list<<endl;
	cout<<"find(bin_list)="<<find(bin_list)<<endl;
	
	cout<<"x="<<x<<endl;
    cout<<"y="<<y<<endl;
    cout<<"e="<<e<<endl;
	cout<<"rem(x,y)="<<rem(x,y)<<endl;
	cout<<"rem(e,y)="<<rem(e,y)<<endl;
	cout<<"ren(y,e)="<<rem(y,e)<<endl;
	cout<<"fact(3)="<<fact(3)<<endl;
	cout<<"binom(5,2)="<<binom(5,2)<<endl; 

  	return 0;
}
