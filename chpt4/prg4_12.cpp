#include "stdafx.h"
#include "itpp/itcomm.h"
#include <itpp/comm/turbo.h>


using std::cin;
using std::cout;
using std::endl;
using namespace itpp;

int main(int argc, char* argv[])
{
 // 类的声明:
QPSK qpsk;                       // QPSK调制
AWGN_Channel awgn_channel;      // AWGN信道
BERC berc;                       // 计算误码率
Turbo_Codec turbo;

// 标量和矢量的声明和初始化: 
int Number_of_bits=400;             // 比特数目
int ml=2;                          // 调制电平数，QPSK ml=2
int i;
int constraint_length ;                // turbo码内交织序列约束长度
double Es, Eb;                      // 符号能量和比特能量
double code_rate;                   // turbo码编码码率
vec EbN0dB, EbN0, N0, noise_variance, bit_error_rate, received_soft_bits; 
// vec 是一个存储double型变量的向量
bvec transmitted_bits, encoded_bits, received_bits;                 
// bvec 是一个存储二进制比特的向量
cvec transmitted_symbols, received_symbols;          
// cvec 是一个存储double_complex型变量的向量
ivec gen(2);                        // Turbo码生成多项式

// 设置Turbo码的参数 ：
gen(0) = 013;
gen(1) = 015;
constraint_length=4;       
code_rate = (double)(Number_of_bits / (3.0 * Number_of_bits + 4.0 *(constraint_length-1))); 

// 设置turbo码内交织序列 
ivec interleaver_sequence = wcdma_turbo_interleaver_sequence( 40 );
turbo.set_parameters(gen, gen, constraint_length,interleaver_sequence);

// 计算噪声单边带功率谱密度：
Es = 1.0;                                   // 每个QPSK符号能量
Eb = Es / (ml * code_rate);                     // 每个信道的比特能量
EbN0dB = linspace(0.0,4.0,5);                  // Eb/N0 的值dB，从0dB到4dB,共5个值
EbN0 = inv_dB(EbN0dB);                     // Eb/N0 线性的值 
N0 = Eb * pow(EbN0,-1.0);                    // N0 复值噪声的方差  

bit_error_rate.set_size(EbN0dB.length(),false);     // 为结果向量分配内存空间. 
RNG_randomize();                           // 产生随机数:

// 开始仿真：
for (i=0; i<EbN0dB.length(); i++) {                    // 遍历所有EbN0dB值:
 
transmitted_bits = randb(Number_of_bits);             // 产生随机二进制数0 1 
turbo.encode(transmitted_bits,encoded_bits);               // turbo编码
transmitted_symbols = qpsk.modulate_bits(encoded_bits);    // QPSK调制

awgn_channel.set_noise(N0(i));                       // 设置AWGN信道噪声
received_symbols = awgn_channel(transmitted_symbols);    // 通过AWGN信道
received_soft_bits = qpsk.demodulate_soft_bits(received_symbols ,N0(i)); // QPSK软解调
turbo.decode(received_soft_bits,received_bits);            // turbo译码

// 计算误比特率：
berc.clear();                                       // 清除误比特率计数器
berc.count(transmitted_bits,received_bits);               // 计算误比特数
bit_error_rate(i) = berc.get_errorrate();         // 将误比特率保存在结果向量中
}

//输出结果:
cout << endl;
cout << "EbN0dB = " << EbN0dB << " [dB]" << endl;
cout << "BER = " << bit_error_rate << endl;
cout << endl;


	return 0;
}
