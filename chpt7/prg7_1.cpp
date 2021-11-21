#include "stdafx.h"
#include "itpp/itcomm.h"

using std::cin;
using std::cout;
using std::endl;
using namespace itpp;

int main(int argc, char* argv[])
{
  //变量声明
  int i, Number_of_bits;                            
  double Ec, Eb;
  vec EbN0dB, EbN0, N0, noise_variance, bit_error_rate; //vec 是浮点型向量
  bvec transmitted_bits, received_bits;                 //bvec是二进制型向量
  cvec transmitted_symbols, received_symbols,ofdm_mod_symbols,ofdm_demod_symbols;   //cvec 是复向量

  //类声明:
  OFDM ofdm;                     //定义OFDM类
  QPSK qpsk;                     //定义QPSK类
  AWGN_Channel awgn_channel;     //定义高斯白噪声信道类
  it_file ff;                    //用于存储结果的文件
  BERC berc;                     //计算误码率
  Real_Timer tt;                 //计时器用于测量运行时间
 
  //计时器开始
  tt.tic();

  //变量初始化：
  Ec = 1.0;                      
  Eb = Ec / 2.0;                 
  EbN0dB = linspace(0.0,9.0,10); 
  EbN0 = inv_dB(EbN0dB);        
  N0 = Eb * pow(EbN0,-1.0);      
  Number_of_bits = 128*2*6;    
  ofdm.set_parameters(128,32,1); //设置OFDM的参数，FFT的点数是，CP循环前缀的长度是，采样点为（默认值）

  //为结果申请空间， "false" 表示不要复制旧的变量到新开辟的空间
  bit_error_rate.set_size(EbN0dB.length(),false);

  //随机数产生器:
  RNG_randomize();

  for (i=0; i<EbN0dB.length(); i++) {

    cout << "Now simulating Eb/N0 value number " << i+1 << " of " << EbN0dB.length() << endl;

    transmitted_bits = randb(Number_of_bits);

    transmitted_symbols = qpsk.modulate_bits(transmitted_bits);

	ofdm.modulate(transmitted_symbols,ofdm_mod_symbols);   //ofdm调制，IFFT变换，加循环前缀
    awgn_channel.set_noise(N0(i));

    received_symbols = awgn_channel(ofdm_mod_symbols);

	ofdm_demod_symbols = ofdm.demodulate(received_symbols); //ofdm解调，FFT变换，去循环前缀
received_bits = qpsk.demodulate_bits(ofdm_demod_symbols);

    //计算误码率:
    berc.clear();                               
    berc.count(transmitted_bits,received_bits);
    bit_error_rate(i) = berc.get_errorrate(); 

  }
  tt.toc_print();
  //显示结果:
  cout << endl;
  cout << "EbN0dB = " << EbN0dB << " [dB]" << endl;
  cout << "BER = " << bit_error_rate << endl;
  cout << "Saving results to ./ofdm_result_file.it" << endl;
  cout << endl;

  //保存结果到文件中:
  ff.open("ofdm_result_file.it");
  ff << Name("EbN0dB") << EbN0dB;
  ff << Name("ber") << bit_error_rate;
  ff.close();

	return 0;
}

