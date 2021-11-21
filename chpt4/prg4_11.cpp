#include <itpp/itbase.h>
#include <itpp/itcomm.h>

using namespace itpp;
using namespace std;

int main()
{
  bvec uncoded_bits = "0 1 1 0 1 0 1", tail_bits(4), decoded_bits(11);
  bmat parity_bits(11, 1);
  vec received_systematic(11), extrinsic_output(11), L(11), symbols;
  mat received_parity(11, 1);

  double Ec = 1.0, N0 = 0.00000000000000000000000000000000000000001, Lc = 4.0 * sqrt(Ec) / N0;
  Normal_RNG noise(0.0, N0 / 2.0);

  BPSK bpsk;
  Rec_Syst_Conv_Code rscc;
 
  rscc.set_generator_polynomials("25 23", 5);//设定生成多项式
  rscc.set_awgn_channel_parameters(Ec, N0);//AWGN信道参数设定

  rscc.encode_tail(uncoded_bits, tail_bits, parity_bits);        //tail 编码
  bpsk.modulate_bits(concat(uncoded_bits, tail_bits), symbols);  //bpsk调制
  received_systematic = symbols + noise(11);                  //添加噪声

  bpsk.modulate_bits(parity_bits.get_col(0), symbols);
  received_parity.set_col(0, symbols + noise(11));

  vec extrinsic_input = zeros(11);
  rscc.map_decode(received_systematic, received_parity, extrinsic_input, 
                  extrinsic_output);         //最大后验概率(MAP) 译码

  L = Lc * received_systematic + extrinsic_output;
  for (int k = 0; k < 11; k++) {
    (L(k) > 0) ? (decoded_bits(k) = bin(0)) : (decoded_bits(k) = bin(1));
  }
  cout << "uncoded_bits = " << uncoded_bits << endl;
  cout << "tail_bits = " << tail_bits << endl;
  cout << "decoded_bits = " << decoded_bits << endl;

  return 0;
}
