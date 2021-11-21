#include "stdafx.h"
#include "itpp/itcomm.h"
#include <itpp/comm/interleave.h>

using std::cin;
using std::cout;
using std::endl;
using namespace itpp;

int main(int argc, char* argv[])
{
	BPSK bpsk;
	bvec bits = "0 1 1 0 0 0 1 1 1 1 0 0 1 0 0 1";
	vec symbols = bpsk.modulate_bits(bits);
	vec interleaved_symbols, deinterleaved_symbols;
	Block_Interleaver<double> block_interleaver(4,4);
	block_interleaver.interleave(symbols , interleaved_symbols);
	block_interleaver.interleave(interleaved_symbols , deinterleaved_symbols);
	cout<<bits<<endl;
	cout<<symbols<<endl;
	cout<<interleaved_symbols<<endl;
	cout<<deinterleaved_symbols<<endl;


	return 0;
}
