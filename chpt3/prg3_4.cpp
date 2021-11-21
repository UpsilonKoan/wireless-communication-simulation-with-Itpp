#include "stdafx.h"
#include "itpp/itcomm.h"
#include <itpp/comm/channel.h>
#include <iostream>
#include <fstream>

using namespace std;
using std::cin;
using std::cout;
using std::endl;
using namespace itpp;

int main(int argc, char* argv[])
{
    int n=10000;
	FIR_Fading_Generator fad1(0.5);
	fad1.set_filter_length(100);
	fad1.init();
    cvec tc;
	fad1.generate(n,tc);
   ofstream f1("FIR_Fading_Generator.txt");
   for(int i=0;i<n;i++)
   {
	   f1<<tc(i)<<endl;
   }
   f1.close();		
	system("Pause");
	return 0;
}
