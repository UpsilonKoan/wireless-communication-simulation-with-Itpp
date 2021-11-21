#include "stdafx.h"
#include "itpp/itcomm.h"

using std::cin;
using std::cout;
using std::endl;
using namespace itpp;

int main(int argc, char* argv[])
{
	int L=601;
	double pi=3.14159;
	double interval=0.001;
	vec t(L),x(L),r(L);
	cvec y(L);
	r=randn(L);
	for(int i=0;i<L;i++)
   {
	  t(i)=i*interval;
	  x(i)=sin(2*pi*50*t(i))+sin(2*pi*120*t(i));
	 // y(i)=x(i)+2*r(i);
     y(i)=x(i);

   }
  int N=512;
  cvec Y;
     Y=fft(y,N);

	 //Declarations of classes:
    it_file ff; 
  //Save the results to file:
  ff.open("fft_file.it");
  //ff << Name("num_p") << num_p;
  ff << Name("time_signal") << y;
  ff << Name("freq_signal") <<sqr(Y);

  ff.close();

	return 0;
}
