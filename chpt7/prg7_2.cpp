#include "stdafx.h"
#include "itpp/itcomm.h"
#include "itpp/comm/spread.h"

using std::cin;
using std::cout;
using std::endl;
using namespace itpp;

int main(int argc, char* argv[])
{
    //产生扩频码
    vec spreading_code = "-1 1 -1 1";

    //初始化Spread_1d类
    Spread_1d spread_1d(spreading_code);

    //产生传输符号
    bvec transmitted_bits = randb(10);
    BPSK bpsk;
    vec transmitted_symbols = bpsk.modulate_bits(transmitted_bits);

    //对序列进行扩频
    vec transmitted_signal = spread_1d.spread(transmitted_symbols);

    vec received_signal = transmitted_signal;

    //对接收信号解扩
    vec received_symbols  = spread_1d.despread(received_signal,0);

    //解调信号
    bvec received_bits = bpsk.demodulate_bits(received_symbols);

    //显示结果
	cout<< "transmitted_bits = " << transmitted_bits <<endl;
	cout<< "transmitted_symbols = " << transmitted_symbols <<endl;
	cout<< "transmitted_signal = " << transmitted_signal <<endl;
	cout<< "received_symbols = " << received_symbols <<endl;
    cout<< "received_bits = " << received_bits <<endl;
	return 0;
}
