#include "stdafx.h"
#include "itpp/itcomm.h"
#include <itpp/comm/hammcode.h>
using std::cin;
using std::cout;
using std::endl;
using namespace itpp;
int main(int argc, char* argv[])
{
		BPSK bpsk; //BPSK类的实例化
		bvec uncoded_bit,coded_bit,decoded_bit,decoded_bit_bpsk;
		vec tx_signal, rx_signal;
        Extended_Golay golay;
		uncoded_bit="1 0 0 0 0 0 0 0 0 0 0 0";
	
		golay.encode(uncoded_bit, coded_bit);
		tx_signal=bpsk.modulate_bits(coded_bit); //调制
		rx_signal = tx_signal + sqrt(0.5)*randn(tx_signal.size());//加噪声
		bpsk.demodulate_bits(rx_signal,decoded_bit_bpsk);
		golay.decode(decoded_bit_bpsk,decoded_bit); //vec型数据译码
		cout<<" uncoded_bit="<< uncoded_bit<<endl;
		cout<<" coded_bit ="<<coded_bit<<endl;
		cout<<"decoded_bit_bpsk=" <<decoded_bit_bpsk<<endl;
		cout<<"decoded_bit ="<<decoded_bit<<endl;

		system("Pause");
}
