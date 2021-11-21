#include "stdafx.h"
#include <itpp/itbase.h>
#include <itpp/itcomm.h>
#include <complex>
typedef std::complex<float> complex;
using namespace itpp;

//These lines are needed for use of cout and endl
using std::cin;
using std::cout;
using std::endl;

int main()
{
	mat a;
	int finish,i,j;
	int antenna_num;
	mat U,V;//svd分解
	vec s;//奇异值向量
	vec alpha;
	vec power;
	double total_power,power_sum;

	antenna_num = 8;
	total_power = 4.33;
	a.set_size(antenna_num,antenna_num);//产生随机矩阵a
    for(i=0;i<a.rows();i++)
	    for(j=0;j<a.cols();j++)
			a(i,j)=randn();
	U.set_size(antenna_num,antenna_num,false);
	V.set_size(antenna_num,antenna_num,false);
	s.set_length(antenna_num);
	alpha.set_length(antenna_num);
	power.set_length(antenna_num);
	
	svd(a,U,s,V);//svd分解
	alpha = s;
	power_sum=0;
	power = itpp::waterfilling(alpha,total_power);
	for(i=0;i<antenna_num;i++)
		power_sum = power_sum+power(i);
	cout<<"alpha"<<endl<<alpha<<endl<<endl;
	cout<<"power"<<endl<<power<<endl<<endl;
	cout<<"total_power"<<endl<<total_power<<endl<<endl;
    cin >> finish;
    return 0;
}
