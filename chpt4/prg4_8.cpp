#include "stdafx.h"
#include "itpp/itcomm.h"

using std::cin;
using std::cout;
using std::endl;
using namespace itpp;

int main(int argc, char* argv[])
{
	BPSK bpsk;
    Convolutional_Code   code;
	ivec generator(3); 
	generator(0)=0133; 
	generator(1)=0165; 
	generator(2)=0171; 
	code.set_generator_polynomials(generator, 7); // �������ɶ���ʽ

	bvec bits=randb(100), encoded_bits, decoded_bits;
	vec  tx_signal,  rx_signal;

	code.encode_tail(bits, encoded_bits); //   ��������
	tx_signal = bpsk.modulate_bits(encoded_bits);// ����
	rx_signal = tx_signal + sqrt(0.5)*randn(tx_signal.size());    //������
	code.decode_tail(rx_signal, decoded_bits);//����������

	cout<<"uncoded_bits="<<bits<<endl;
	cout<<"encoded_bits="<<encoded_bits<<endl;
	cout<<"tx_signal="<<tx_signal<<endl;
	cout<<"rx_signal="<<rx_signal<<endl;
	cout<<"decoded_bits="<<decoded_bits<<endl;

	return 0;
}
