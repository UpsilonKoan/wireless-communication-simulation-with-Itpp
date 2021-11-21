#include "stdafx.h"
#include "itpp/itcomm.h"
#include "stdio.h"

using std::cin;
using std::cout;
using std::endl;
using namespace itpp;

int main()
{
	vec my_vec_1("1 0 1"); //must be 0 or 1
	
	//to_bvec
	bvec my_bvec;
	my_bvec = to_bvec(my_vec_1);
    cout<<" my_vec_1="<<my_vec_1<<endl;
	cout<<" my_bvec=" <<my_bvec <<endl;

	//to_svec
	vec	my_vec_2("1 1 1");   //must be 0 or 1
	svec my_svec;
   cout<<"my_vec_2="<<my_vec_2<<endl;
	my_svec = to_svec(my_vec_2);
	cout<<"my_svec=" <<my_svec <<endl;

	//to_ivec
	vec my_vec_3("2.3 4.8 5");
	ivec my_ivec;
	my_ivec = to_ivec(my_vec_3);
	cout<<"my_vec_3="<<my_vec_3<<endl;
	cout<<" my_ivec=" <<my_ivec <<endl; //向零取整

	//to_cvec
	vec my_vec_4("2.4 3.6 8");
	cvec my_cvec;
	my_cvec = to_cvec(my_vec_4);
    cout<<" my_vec_4="<<my_vec_4<<endl;
	cout<<"my_cvec=" <<my_cvec <<endl;

	vec my_vec_real("1 2 3"); 
	vec my_vec_imag("3 2 1");
	cvec my_cvec1;
	my_cvec1 = to_cvec(my_vec_real,my_vec_imag);

	cout<<" my_vec_real="<<my_vec_real<<endl;
    cout<<" my_vec_imag="<<my_vec_imag<<endl;
	cout<<"my_cvec1=" <<my_cvec1 <<endl;

	//to_ivec 将一个整型数转换成为一个ivec
	int int_s = 5;
	ivec my_ivec1;
	my_ivec1 = to_ivec(int_s);
	cout<<"int_s="<<int_s<<endl;
	cout<<" my_ivec1="<<my_ivec1<<endl;
	
	return 0;
}
