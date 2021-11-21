#include "stdafx.h"
#include "itpp/itcomm.h"
#include <itpp/comm/hammcode.h>

using std::cin;
using std::cout;
using std::endl;
using namespace itpp;
int main(int argc, char* argv[])
{
	    int m=3;
        bvec uncoded_bit, coded_bit;
		bmat h;
		Hamming_Code hamming(m);  // 类的实例化，初始化参数m=3
		uncoded_bit="1 0 0 0 0 0 0 0 0 0 0 0";  // 输入12位的信息比特
		h=hamming.get_H();  // 获得校验矩阵
	    hamming.encode(uncoded_bit,coded_bit);
		cout<<"uncoded_bit="<<uncoded_bit<<endl;
		cout<<"h="<<h<<endl;        
    	cout<<"coded_bit"<<coded_bit<<endl;
 	    system("Pause");
}
