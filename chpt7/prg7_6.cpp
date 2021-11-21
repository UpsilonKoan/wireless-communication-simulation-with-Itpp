#include "stdafx.h"
#include <itpp/itbase.h>
#include <itpp/itcomm.h>
using namespace itpp;

//These lines are needed for use of cout and endl
using std::cin;
using std::cout;
using std::endl;

int main()
{
    cout<<"******************LFSR*****************"<<endl;
	LFSR lfsr;
	bvec Connection="1 0 1 1";//连接多项式
	bvec State="0 0 1";//初始状态
	lfsr.set_connections(Connection);//设置连接多项式
	lfsr.set_state(State);//设置初始状态
	cout<<lfsr.get_state();
	cout<<"连接多项式："<<Connection<<endl;
	cout<<"寄存器长度"<<lfsr.get_length()<<endl;
	bvec Sequence(7);
	int i;
	for(i = 0;i < Sequence.length();i ++)
	{
		Sequence(i)=lfsr.shift();
		cout<<"step  "<<"寄存器状态      "<<"output  "<<endl;
		cout<<i<<"     "<<lfsr.get_state()<<"    "<<Sequence(i)<<endl;
	}
	cout<<endl<<"序列:"<<Sequence<<endl;
	return 0;
}
