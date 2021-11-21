#include "stdafx.h"
#include "itpp/itcomm.h"

using std::cin;
using std::cout;
using std::endl;
using namespace itpp;

int main(int argc, char* argv[])
{
	BPSK bpsk;
 	bvec bits = "0 1 1 0 0 0 1 1 1 1 0 0 1 0 0 1";
  	vec symbols = bpsk.modulate_bits(bits);
    ivec interleave_sequence = sort_index(randu(16));

  	Sequence_Interleaver<double> sequence_interleaver(interleave_sequence);
  	//sequence_interleaver.;
  	vec interleaved_symbols = sequence_interleaver.interleave(symbols);

	cout<<bits<<endl;
	cout<<symbols<<endl;
	cout<<interleave_sequence<<endl;
	//cout<<sequence_interleaver.get_interleaver_sequence()<<endl;
	cout<<interleaved_symbols<<endl;

	

	return 0;
}
