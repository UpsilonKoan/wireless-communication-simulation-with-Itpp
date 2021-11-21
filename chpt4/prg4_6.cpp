#include "stdafx.h"
#include "itpp/itcomm.h"
#include <itpp/comm/reedsolomon.h>
using std::cin;
using std::cout;
using std::endl;
using namespace itpp;
int main(int argc, char* argv[])
{
	BPSK bpsk;
	int m=3,t=1;
	bvec uncoded_bit,coded_bit,decoded_bit,decoded_bit_bpsk;
	vec tx_signal, rx_signal;
	Reed_Solomon rs(m,t,false);	
	uncoded_bit="1 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
	rs.encode(uncoded_bit,coded_bit);
	tx_signal=bpsk.modulate_bits(coded_bit);
	rx_signal = tx_signal + sqrt(0.5)*randn(tx_signal.size());
    bpsk.demodulate_bits(rx_signal,decoded_bit_bpsk);
	decoded_bit=rs.decode(decoded_bit_bpsk);
	cout<<" uncoded_bit ="<<uncoded_bit<<endl;
	cout<<"coded_bit="<<coded_bit<<endl;
    cout<<"decoded_bit_bpsk="<<decoded_bit_bpsk<<endl;
	cout<<"decoded_bit="<<decoded_bit<<endl;

	system("Pause");
}
