#include "stdafx.h"
#include "itpp/itcomm.h"
#include <itpp/comm/crc.h>
#include <itpp/itbase.h>


using std::cin;
using std::cout;
using std::endl;
using std::string;


using namespace itpp;

int main(int argc, char* argv[])
{
	//Declarations of scalars and vectors: 
  CRC_Code crc(string("CRC-4"));//运用字符串值CRC-4设定生成CRC
  //标准多项式
  bvec bits = randb(10),coded_bits,decoded_bits,rec_bits;//设定存储变量，
                                               //生成信息比特
  bool error;
  coded_bits = crc.encode(bits);   //计算并且添加检验字段，生成CRC码 
  rec_bits =coded_bits;
  error = crc.decode(rec_bits, decoded_bits);//校验正确性
  cout<<"输入比特="<<bits<<endl;
  cout<<"CRC编码="<<coded_bits<<endl;
  cout<<"译码结果"<<decoded_bits<<endl;

  
  return 0;
}
