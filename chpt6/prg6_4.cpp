// 8-PSK调制解调
#include <itpp/itcomm.h>
using namespace itpp;
using namespace std;
int main()
{
  RNG_reset(12345);   //设置随机数种子
  const int no_symbols = 5;  //发送符号数
  const double N0 = 0.1;  //噪声功率
{
    cout << endl << "8-PSK" << endl;
    PSK psk(8);  // 8-PSK调制
    int bps = round_i(psk.bits_per_symbol()); //每个符号的比特数
    bvec tx_bits = randb(no_symbols * bps); //发送比特
    ivec tx_sym_numbers = randi(no_symbols, 0, pow2i(bps) - 1); //发送符号
    cvec noise = sqrt(N0) * randn_c(no_symbols); //噪声

    cvec tx_symbols = psk.modulate(tx_sym_numbers);//符号调制
    cvec rx_symbols = tx_symbols + noise; //接收符号
    ivec dec_sym_numbers = psk.demodulate(rx_symbols);//符号解调

    cout << "* modulating symbol numbers:" << endl;
    cout << "  tx_sym_numbers  = " << tx_sym_numbers << endl;
    cout << "  tx_symbols      = " << tx_symbols << endl;
    cout << "  rx_symbols      = " << rx_symbols << endl;
    cout << "  dec_sym_numbers = " << dec_sym_numbers << endl;

    tx_symbols = psk.modulate_bits(tx_bits);//比特调制
    rx_symbols = tx_symbols + noise; //接收符号
    bvec decbits = psk.demodulate_bits(rx_symbols);//比特硬解调
    vec softbits_approx = psk.demodulate_soft_bits(rx_symbols, N0, APPROX); //AWGN信道的软比特近似解调
    vec softbits = psk.demodulate_soft_bits(rx_symbols, N0, LOGMAP);
 //AWGN信道的软比特Log-MAP解调

    cout << "* modulating bits:" << endl;
    cout << "  tx_bits         = " << tx_bits << endl;
    cout << "  tx_symbols      = " << tx_symbols << endl;
    cout << "  rx_symbols      = " << rx_symbols << endl;
    cout << "  decbits         = " << decbits << endl;
    cout << "  softbits        = " << softbits << endl;
    cout << "  softbits_approx = " << softbits_approx << endl << endl;
  }
}
