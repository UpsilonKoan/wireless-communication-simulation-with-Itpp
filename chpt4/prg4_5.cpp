#include "stdafx.h"
#include "itpp/itcomm.h"
#include <itpp/comm/reedsolomon.h>
using std::cin;
using std::cout;
using std::endl;
using namespace itpp;
int main(int argc, char* argv[])
{
		int m=3,t=1;
		bvec uncoded_bit,coded_bit,decoded_bit;
		Reed_Solomon rs(m,t,false);	//RS码类的实例化
		uncoded_bit="1 0 0 0 0 0 0 0 0 0 0 0 0 0 0"; //输入15位的信息比特
		rs.encode(uncoded_bit,coded_bit);
		rs.decode(coded_bit,decoded_bit);
		cout<<uncoded_bit<<endl;
    	cout<<coded_bit<<endl;
		cout<<decoded_bit<<endl;
		system("Pause");
}
