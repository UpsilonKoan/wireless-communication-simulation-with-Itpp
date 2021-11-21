#include "stdafx.h"
#include "itpp/itcomm.h"
#include "itpp/comm/spread.h"

using std::cin;
using std::cout;
using std::endl;
using namespace itpp;

int main(int argc, char* argv[])
{
   //产生I Q支路扩频码
   vec spreading_code_I = "1 1 1 1";
   vec spreading_code_Q = "1 1 1 1";

   //初始化spread_2d类
   Spread_2d spread_2d(spreading_code_I,spreading_code_Q);

   //产生传输比特
   bvec transmitted_bits = randb(6);

   //QPSK调制
    QPSK qpsk;                

    cvec transmitted_symbols = qpsk.modulate_bits(transmitted_bits);

	//扩频调制
	cvec transmitted_signal = spread_2d.spread(transmitted_symbols);	

    cvec received_signal = transmitted_signal;

	//解扩信号
	cvec received_symbols  = spread_2d.despread(received_signal,0);

	//解调
    bvec received_bits = qpsk.demodulate_bits(received_symbols);

    //显示结果
	cout<< "transmitted_bits = " << transmitted_bits <<endl;
	cout<< "transmitted_symbols = " << transmitted_symbols <<endl;
	cout<< "transmitted_signal = " << transmitted_signal <<endl;
	cout<< "received_symbols = " << received_symbols <<endl;
    cout<< "received_bits = " << received_bits <<endl;
	return 0;
}
