#include "stdafx.h"
#include "itpp/itcomm.h"
#include <itpp/base/mat.h> 
#include <itpp/base/vec.h> 
#include <iostream>


using std::cin;
using std::cout;
using std::endl;
using namespace itpp;
using namespace std;
int main()
{
  int  i,j,m; 
  int  Number_of_bits=400000; // total number of bits to be simulated
  int  Kc=4; //Number of antenas 
  int  M=4;  //M-QAM
  int  Rm=round_i(log2(double(M)));//Number of bits per symbol
  vec  N0,bit_error_rate;
  vec  EbN0dB = linspace(0.0,25.0,12); 
  vec  EbN0 = inv_dB(EbN0dB);
  //Transmitted vectors
  bvec Transmit_bit_bvec; //bit vector to be transmitted
  cvec Transmit_symbol_cvec; // symbol vector after QAM  
  cmat Transmit_symbol_cmat; // symbol matrix after QAM 
  cmat Transmit_symbol_colum_cmat;//transmited column vector for each loop
  
  //received vectors (complex)
  cmat Receive_symbol_colum_cmat; //received column vector after fading channel
  cmat Receive_symbol_colum_cmat_1;//receive vector multiplied by detection matrix
  cvec Receive_symbol_cvec;//vector for decision
  bvec Receive_bit_colum_bvec;//restored bit vector for each loop
  bvec Receive_bit_bvec;//restored vector of bits
  //Noise matrix
  mat  Noise_mat_r,Noise_mat_i;
  cmat Noise_cmat;
  //Fading channel matrix
  mat  H_mat_r,H_mat_i;
  cmat H_cmat;//Rayleigh fading channel matrix (complex)

  double spow;//signal power
  it_file ff;                   
  BERC berc;                     
  Real_Timer tt;                 
  Rayleigh_RNG rayleigh_rng;
  QAM qam;
qam.set_M(M);
tt.tic();
   
  //set the size of matrix
  bit_error_rate.set_size(EbN0dB.length(),false);
  N0.set_size(EbN0dB.length(),false);
  Transmit_bit_bvec.set_size(Number_of_bits,false);
  Transmit_symbol_cvec.set_size(Number_of_bits/Rm,false);
  Transmit_symbol_colum_cmat.set_size(Kc,1,false);
  Receive_symbol_colum_cmat_1.set_size(Kc,1,false);
  Receive_symbol_cvec.set_size(Kc,false);
  Receive_bit_colum_bvec.set_size(Kc*Rm,false);
  Receive_bit_bvec.set_size(Number_of_bits,false);

  //constellation 
  cmat c_cmat;
  bmat b_bmat;
  c_cmat.set_size(M,1,false);
  b_bmat.set_size(M,Rm,false);

  for(i=0;i<M;i++){
      bvec temp=dec2bin(Rm,i);//decimal index to binary
	  b_bmat.set_row(i,temp);
	  c_cmat.set_row(i,qam.modulate_bits(b_bmat.get_row(i)));
  }
  cvec x=c_cmat.get_col(0);
  //get the transmitting power
  spow=0;
  for(i=0;i<M;i++)
  {
      spow+=pow(abs(x(i)),2.0);
  }
  spow=spow/M;
RNG_randomize();

/******************************************************************
 Get the transmitted data symbols
******************************************************************/
  //data bits to be transmitted
  Transmit_bit_bvec=randb(Number_of_bits);
  //data symbols after QAM
  Transmit_symbol_cvec=qam.modulate_bits(Transmit_bit_bvec);
  
  // start the simulation loops

  for (m=0; m<EbN0dB.length(); m++) 
  {
	for(i=0;i<Number_of_bits/(Kc*Rm);i++)
	{
		for(j=0;j<Kc;j++)
		{
			Transmit_symbol_colum_cmat(j,0)=Transmit_symbol_cvec(Kc*i+j);
	     }
		/******************************************************************
		go through the channel
		******************************************************************/
	     N0(m)=spow*Kc*pow(EbN0(m),-1.0)/Rm;
		//noise matrix generation
		Noise_mat_r=sqrt(N0(m)/2)*randn(Kc,1);
		Noise_mat_i=sqrt(N0(m)/2)*randn(Kc,1);
		Noise_cmat=to_cmat(Noise_mat_r,Noise_mat_i);
		//fading channel matrix generation
	  	H_mat_r=randn(Kc,Kc);
	     H_mat_i=randn(Kc,Kc);
	     H_cmat=to_cmat(H_mat_r,H_mat_i);
		//received signal vector
		Receive_symbol_colum_cmat=H_cmat*Transmit_symbol_colum_cmat+Noise_cmat;
		/******************************************************************
		obtain the equalization matrix
******************************************************************/

//ZF
cmat F; 
		F=inv(H_cmat);
		
/******************************************************************
		   Detection
		******************************************************************/
        //multiplid by the equalization matrix
         Receive_symbol_colum_cmat_1=F*Receive_symbol_colum_cmat;
		Receive_symbol_cvec=Receive_symbol_colum_cmat_1;
	    //detect x_0,x_1,...,x_Kc-1
         vec Distance_vec;
	     Distance_vec.set_size(M,false);
	    for(j=0;j<Kc;j++)
		{
		  Distance_vec.zeros();
		  for(int ii=0;ii<M;ii++)
		  {
			  Distance_vec(ii)+=sqr(Receive_symbol_cvec(j)-x(ii));
		  }
          int ii_min=min_index(Distance_vec);
		  Receive_symbol_cvec(j)=x(ii_min);
		}

		//QAM demoluation
		Receive_bit_colum_bvec=qam.demodulate_bits(Receive_symbol_cvec);
		//restore the data bit stream
		for(j=0;j<Kc*Rm;j++)
		{
	        Receive_bit_bvec(Rm*Kc*i+j)=Receive_bit_colum_bvec(j); 
	    }
	}
	berc.clear();
    berc.count(Transmit_bit_bvec,Receive_bit_bvec); 
    bit_error_rate(m) = berc.get_errorrate();

  }
  tt.toc_print();
  cout<<"BER="<<bit_error_rate<<endl;
  system("pause");
  return 0;
}


