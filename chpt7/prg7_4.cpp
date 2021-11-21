#include "stdafx.h"
#include "itpp/itcomm.h"

using std::cin;
using std::cout;
using std::endl;
using namespace itpp;

int main(int argc, char* argv[])
{
	bmat GrayCode;
   GrayCode = graycode(4);
   cout<<GrayCode<<endl;



   bvec Sequence1,Sequence2;
Sequence1.set_length (6);//2个序列长度必须一样，不然会出错
Sequence2.set_length (6);
Sequence1 = randb(6);
Sequence2 = randb(6);
cout<<Sequence1<<endl<<Sequence2<<endl;
int HammingDist;
HammingDist = hamming_distance(Sequence1,Sequence2);
cout<<"hamming distance="<<HammingDist<<endl;

bvec Sequence3;
Sequence3.set_length (8);
Sequence3 = randb(8);
cout<<"code: "<<Sequence3<<endl;
int Weight;
Weight = weight(Sequence3);
cout<<"Weight="<<Weight<<endl;


	return 0;
}
