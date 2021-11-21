#include "stdafx.h"
#include "itpp/itcomm.h"
#include <itpp/comm/channel.h>
#include <iostream>
#include <fstream>
using namespace std;
using std::cin;
using std::cout;
using std::endl;
using namespace itpp;

int main(int argc, char* argv[])
{
	//信道每个抽样点变化一次
	ofstream f1("staticchannel.txt");
  	Static_Fading_Generator fad1;   //实例化类确定衰落类
   //初始化衰落产生器
	fad1.set_LOS_power(1.0);      //设置直射径相对功率值
	cvec out;                    //设置输出向量
	out.set_length(10);           //设置每次静态信道输出向量长度

	for(int i=0;i<1000;i++)
	{
		fad1.init();  
		fad1.generate(10,out);        //产生信道系数存放在out里面
		f1<<abs(out)<<endl;
	}
	
	f1.close();
	system("Pause");
	return 0;
}
