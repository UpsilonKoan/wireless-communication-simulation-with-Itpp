#include <itpp/itcomm.h>
using namespace itpp;
using namespace std;
int main()
{
  RNG_reset(12345);//随机数种子
  const int no_symbols = 5; //符号个数
  const double N0 = 0.05;//噪声功率谱密度
{
    cout << endl << "4-PAM(real)" << endl;
    PAM pam(4);
    int bps = round_i(pam.bits_per_symbol());//每符号比特个数

    bvec tx_bits = randb(no_symbols * bps);//发送比特
    ivec tx_sym_numbers = randi(no_symbols, 0, pow2i(bps) - 1);//发送符号
    vec noise = sqrt(N0) * randn(no_symbols);//噪声

    vec tx_symbols = pam.modulate(tx_sym_numbers);//符号调制
    vec rx_symbols = tx_symbols + noise;//接收符号
    ivec dec_sym_numbers = pam.demodulate(rx_symbols);//符号解调

    cout << "* modulating symbol numbers:" << endl;
    cout << "  tx_sym_numbers  = " << tx_sym_numbers << endl;
    cout << "  tx_symbols      = " << tx_symbols << endl;
    cout << "  rx_symbols      = " << rx_symbols << endl;
    cout << "  dec_sym_numbers = " << dec_sym_numbers << endl;

    tx_symbols = pam.modulate_bits(tx_bits);//比特调制
    rx_symbols = tx_symbols + noise;//接收符号
    bvec decbits = pam.demodulate_bits(rx_symbols);//比特硬解调
    vec softbits_approx = pam.demodulate_soft_bits(rx_symbols, N0, APPROX); //AWGN信道的软比特近似解调
    vec softbits = pam.demodulate_soft_bits(rx_symbols, N0, LOGMAP);
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
