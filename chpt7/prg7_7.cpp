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

    cout<<"***************Gold 序列**************"<<endl;
	Gold gold(5);//degree=5
	bvec GoldSequence(31);//Gold序列
	bvec State1="0 0 0 0 1";//2个寄存器的初始状态
	bvec State2="0 0 0 0 1";
	gold.set_state(State1,State2);
	cout<<"Gold序列长度="<<gold.get_sequence_length()<<endl;
	int i;
	for(i=0;i<gold.get_sequence_length();i++)
	{
		GoldSequence(i) = gold.shift();
	}
	cout<<"Gold序列"<<endl;
	cout<<GoldSequence<<endl;
	bmat GoldFamily(33,31);
	GoldFamily = gold.get_family();
	cout<<"GoldFamily"<<endl;
	cout<<GoldFamily<<endl;
return 0;
}
