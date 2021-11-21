#include <itpp/itcomm.h>

using namespace itpp;
using namespace std;

int main()
{
  cout << "====================================" << endl;
  cout << "    Test of convolutional coders    " << endl;
  cout << "====================================" << endl;

  Array<ivec> spectrum, spectrum_fast;
  spectrum.set_size(2);
  spectrum_fast.set_size(2);
 
  Convolutional_Code code;
  BPSK bpsk;
  BERC berc;

  const int no_bits = 2500;
  const int packet_size = 500;

  int coded_packet_size;
  bvec bits, tail_coded_bits, tail_decoded_bits, tailbite_coded_bits,
  tailbite_decoded_bits, trunc_coded_bits, trunc_decoded_bits;
  vec symbols;
  ivec dist_profile;

  ivec G(2);
  G(0) = 0133;
  G(1) = 0171;
  int L = max(int2bits(G(0)), int2bits(G(1))); // L = 7

  code.set_generator_polynomials(G, L);  //设定生成多项式

  cout << "------------------------------------------------------------------------" << endl;
  cout << "1) Rate 1/2 code" << endl;
  cout << "------------------------------------------------------------------------" << endl;

  cout << "Catastrophic test = " << code.catastrophic() << endl; //检查是否是恶性码，如果是返回1
  cout << "Code rate         = " << code.get_rate() << endl << endl;//获取码率R=1/2

  code.calculate_spectrum(spectrum, 10, 10);   //计算谱
  code.fast(spectrum_fast, 10, 10);//Cederwall快速算法
  code.distance_profile(dist_profile, 10); //距离谱谱

  cout << "Spectrum:" << endl;
  cout << "* Ad = " << spectrum(0) << endl;   //重量谱
  cout << "* Cd = " << spectrum(1) << endl;   //信息重量谱

  cout << "Spectrum, fast:" << endl;
  cout << "* Ad = " << spectrum_fast(0) << endl;
  cout << "* Cd = " << spectrum_fast(1) << endl << endl;

  cout << "Distance profile  = " << dist_profile << endl << endl;

  cout << "Tail method test. Printing 30 bits starting from bit 1400:" << endl;    //tail编译码
  bits = randb(no_bits);
  cout << "* Input bits    = " << bits.mid(1400, 30) << endl;
  tail_coded_bits = code.encode_tail(bits);
  cout << "* Coded bits    = " << tail_coded_bits.mid(1400, 30) << endl;
  bpsk.modulate_bits(tail_coded_bits, symbols);
  tail_decoded_bits = code.decode_tail(symbols);
  cout << "* Decoded bits  = " << tail_decoded_bits.mid(1400, 30) << endl;
  berc.count(bits, tail_decoded_bits);
  cout << "BER = " << berc.get_errorrate() << endl << endl;

  cout << "Tailbite method test. Printing 30 bits starting from bit 1400:"   //tailbite编译码
       << endl;
  cout << "* Input bits    = " << bits.mid(1400, 30) << endl;
  tailbite_coded_bits = code.encode_tailbite(bits);
  cout << "* Coded bits    = " << tailbite_coded_bits.mid(1400, 30) << endl;
  bpsk.modulate_bits(tailbite_coded_bits, symbols);
  tailbite_decoded_bits = code.decode_tailbite(symbols);
  cout << "* Decoded bits  = " << tailbite_decoded_bits.mid(1400, 30) << endl;
  berc.clear();
  berc.count(bits, tailbite_decoded_bits);
  cout << "BER = " << berc.get_errorrate() << endl << endl;

  cout << "Trunc method test. Printing 30 bits starting from bit 1400:" << endl;  //Trunc编译码
  cout << "* Input bits    = " << bits.mid(1400, 30) << endl;
  trunc_coded_bits.set_size(0);
  for (int i = 0; i < no_bits / packet_size; i++) {
    trunc_coded_bits = concat(trunc_coded_bits,
                              code.encode_trunc(bits.mid(i * packet_size,
                                                         packet_size)));
  }
  cout << "* Coded bits    = " << trunc_coded_bits.mid(1400, 30) << endl;
  bpsk.modulate_bits(trunc_coded_bits, symbols);
  trunc_decoded_bits.set_size(0);
  coded_packet_size = round_i(packet_size / code.get_rate());
  for (int i = 0; i < no_bits / packet_size; i++) {
    trunc_decoded_bits =
      concat(trunc_decoded_bits,
             code.decode_trunc(symbols.mid(i * coded_packet_size,
                                           coded_packet_size)));
  }
  cout << "* Decoded bits  = " << trunc_decoded_bits.mid(1400, 30) << endl;
  berc.clear();
  berc.count(bits, trunc_decoded_bits);
  cout << "BER = " << berc.get_errorrate() << endl << endl;
  return 0;
}

