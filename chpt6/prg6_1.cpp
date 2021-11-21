#include <itpp/itcomm.h>
using namespace itpp;
using namespace std;
int main()
{
  RNG_reset(12345); //设置随机数种子
  const int no_symbols = 5;//5个符号
  const double N0 = 0.1;
  {
    cout << endl << "Modulator_1D (configured as BPSK)" << endl;
    Modulator_1D  mod("1.0 -1.0", "0 1"); //1D实调制器
    int bps = round_i(mod.bits_per_symbol());

    bvec tx_bits = randb(no_symbols * bps);//产生发送bit
    ivec tx_sym_numbers = randi(no_symbols, 0, pow2i(bps) - 1);//发送符号
    vec noise = sqrt(N0) * randn(no_symbols);//产生噪声

    vec tx_symbols = mod.modulate_bits(tx_bits);//比特调制，得到调制后的符号
    vec rx_symbols = tx_symbols + noise;//接收到的有噪声的符号
    bvec decbits = mod.demodulate_bits(rx_symbols);//比特硬解调，得到解调结果

    cout << "* modulating bits:" << endl;
    cout << "  tx_bits         = " << tx_bits << endl;
    cout << "  tx_symbols      = " << tx_symbols << endl;
    cout << "  rx_symbols      = " << rx_symbols << endl;
    cout << "  decbits         = " << decbits << endl;

    tx_symbols = mod.modulate(tx_sym_numbers);//符号调制
    rx_symbols = tx_symbols + noise;
    ivec dec_sym_numbers = mod.demodulate(rx_symbols);//符号解调

    cout << "* modulating symbol numbers:" << endl;
    cout << "  tx_sym_numbers  = " << tx_sym_numbers << endl;
    cout << "  tx_symbols      = " << tx_symbols << endl;
    cout << "  rx_symbols      = " << rx_symbols << endl;
    cout << "  dec_sym_numbers = " << dec_sym_numbers << endl;
}
}
