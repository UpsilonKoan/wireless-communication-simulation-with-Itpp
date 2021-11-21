//----- VectorChannel.cpp -----
// Programmed by LTE Simulation Group of IMC, SWJTU
#include "stdafx.h"
#include "itpp/itcomm.h"
#include "VectorChannel.h"
#include<iostream>
#include<fstream>
#include<math.h>
#include<complex>
typedef std::complex<double> complex;

using std::cin;
using std::cout;
using std::endl;
using namespace itpp;
using namespace std;

#define   N   4//N个点 
#define   T   3  //T次拟合 
#define   W   1 //权函数
#define   pi  3.141592653
#define   PRECISION   0.0000001

// *************************************** Class CVectorChannel Definition *****************************************

//--- Definition of CVectorChannel constructor
CVectorChannel::CVectorChannel()
{
	mNumOfMSAntenna=2;                              // MS antenna Number
	mNumOfBSAntenna=2;                              // BS antenna Number
	mAntennaSpacingBS=4;                            // d
	mAntennaSpacingMS=0.5;                          // d
	AOAPerPathMS=22.5;                              // 22.5 or 67.5or-67.5
	AODPerPathBS=50;                                // 50  20
	mAngleSpreadPerPathMS="35,35,35,35,35,35";      // ChannelCoefficient is calculated by path
	mAngleSpreadPerPathBS=" 2, 2, 2, 2, 2, 2";
	NumOfSample=1024;
	mModelCase=1;
	mNumOfPath=6;
	mVelocityOfMS="30";
	mFrequencyOfCarrier=1.8e9;
	mDoppler=50;
}

//--- Definition of CVector destructor
CVectorChannel::~CVectorChannel()
{}

//--- Definition of SetCVectorChannelParameter()
void CVectorChannel::SetCVectorChannelParameter()
{
 // Set up the parameters of Function Member of Vector Channel
}

//--- Definition of ChannelCorrelation()
mat CVectorChannel::ChannelCorrelation(double AngleSpread,double Angle,double RatioDistance,int NumOfAntenna)
  {
    double PowerTotalPath=0;
    double anglespread1=720;
    int MacroPathNum=1000;
    mat AngleEachPath(1,MacroPathNum);
	mat AngleEachPath1(1,MacroPathNum);
    mat FAI(1,MacroPathNum);
	mat PowerEachPath(1,MacroPathNum);
	mat correlation(NumOfAntenna,NumOfAntenna);
	cmat temp_M1(NumOfAntenna,1);
	cmat temp_M2(1,NumOfAntenna);
    cmat correlation1(NumOfAntenna,NumOfAntenna);
	cmat correlation2(NumOfAntenna,NumOfAntenna);
	correlation2.zeros();//very important to init....

	for(int i=0;i<1000;i++)
	{
		  AngleEachPath1(0,i)=Angle-anglespread1+2*anglespread1*(i+1)/1000;
		  AngleEachPath(0,i)=2*pi*(Angle-anglespread1+2*anglespread1*(i+1)/1000)/360;  	  
		  FAI(0,i)=RatioDistance*sin(AngleEachPath(i));	  
    }

	for(int i=0;i<1000;i++)
	{
 		PowerEachPath(0,i)=1/(AngleSpread*sqrt(2.0))*
			exp(-sqrt(2.0)*abs(AngleEachPath1(i)-Angle)/AngleSpread)*2*720/1000;  
	}

	for(int i=0;i<1000;i++)
	{
	    PowerTotalPath=PowerTotalPath+PowerEachPath(0,i);
	}

	/*ofstream f1("PowerEachPath.txt");
    for(int ff=0;ff<1000;ff=ff+1)
	{
    f1<<(PowerEachPath(0,ff))<<endl;
	}
	cout<<PowerTotalPath<<endl;*/

	cmat temp(NumOfAntenna,NumOfAntenna);

     for(int m=0;m<MacroPathNum;m++)
	 {
	  for(int n=0;n<NumOfAntenna;n++)
	  {
		  temp_M1(n,0).real(cos((FAI(0,m)*2*pi*(n))));
		  temp_M1(n,0).imag(sin((FAI(0,m)*2*pi*(n))));
	  }
       temp_M2=conj(transpose(temp_M1));
	   temp=temp_M1*temp_M2;
	   correlation1=PowerEachPath(0,m)*temp;
	   correlation2=correlation1+correlation2;
	 }
    correlation=(abs(correlation2))/PowerTotalPath;
	//cout<<correlation<<endl;
	return correlation;
  }

// *************************************** Functions used by scalarchannel() START *****************************************

mat init(vec x,vec y,int n) 
{ 
	int   i; 
	mat x_y(4,2);
	for(i=0;i <n;i++) 
	{
		x_y(i,0)=x(i);
		x_y(i,1)=y(i); 
	} 
	return(x_y);
} 

mat get_A(mat  x_y,int   n) 
{
	mat matrixA(N,T+1);
	int   i,j; 
	for(i=0;i <N;i++)
	{
		for(j=0;j <T+1;j++)
		{
			matrixA(i,j)=W*pow(x_y(i,0),j);
		}
	}
	return matrixA;
} 

vec  get_y(mat trans_A,mat x_y,int   n) 
{
	int   i,j;
	double  temp; 
	vec y(4);
	for(i=0;i <n;i++)
	{
		temp=0;
		for(j=0;j <N;j++)
		{
			temp+=trans_A(i,j)*x_y(j,1);
		}
		y(i)=temp;
	}
	return y;
} 

mat cons_formula(mat coef_A,vec  y) 
{
	mat coef_form(T+1,T+2);
	for(int i=0;i <T+1;i++)
	{
		for(int j=0;j <T+2;j++)
		{
			if(j==T+1)
				coef_form(i,j)=y(i);
			else
				coef_form(i,j)=coef_A(i,j);
		}
	}
	return coef_form;
} 

void  convert(mat argu,int n) 
{ 
	int   i,j,k,p,t; 
	double  rate,temp; 
	for(i=1;i <n;i++)
	{
		for(j=i;j <n;j++)
		{
			if(argu(i-1,i-1)==0)
			{
				for(p=i;p <n;p++)
				{
					if(argu(p,i-1)!=0)  break;
				}
				if(p==n)
				{ 
					cout<<"error"<<endl;				
				}
				for(t=0;t <n+1;t++)
				{
					temp=argu(i-1,t); 
					argu(i-1,t)=argu(p,t); 
				    argu(p,t)=temp;
				}
			}
			rate=argu(j,i-1)/argu(i-1,i-1);
			for(k=i-1;k <n+1;k++)
			{
				argu(j,k)-=argu(i-1,k)*rate; 
				if(fabs(argu(j,k)) <=PRECISION) 
					argu(j,k)=0;
			}
		} 
	}
}

vec compute(mat argu,int   n) 
{
	int   i,j;
	double  temp;
	vec root(4);
	for(i=n-1;i>=0;i--)
	{
		temp=argu(i,n);
		for(j=n-1;j>i;j--)
		{
			temp-=argu(i,j)*root(j);
		}
		root(i)=temp/argu(i,i);
	}
	return root;
} 

vec process(vec x_point,vec   y_point) 
{ 
	vec  result;
	mat x_y(4,2),matrix_A(N,T+1),trans_A(T+1,N),coef_A(T+1,T+1),coef_formu(T+1,T+2);
	vec y(T+1);
	x_y=init(x_point,y_point,N);
	matrix_A=get_A(x_y,N);
	trans_A=transpose(matrix_A);
	coef_A=trans_A*matrix_A;
	y=get_y(trans_A , x_y , T+1);
	coef_formu=cons_formula(coef_A , y);
	convert(coef_formu , T+1);
	result=compute(coef_formu , T+1);
	return result;
} 
// *************************************** Functions used by scalarchannel() END *****************************************

//--- Definition of ScalarChannel()
cvec CVectorChannel::ScalarChannel(cvec velocity, double fre, double doppler)
{
	double frequency=1.8e9;
	cvec v(1);
	v(0).imag(0.0);
	v(0).real(30.0);
	vec ve=10e2*abs(mVelocityOfMS)/3600;
	int NumofRealPoint=8;
	double mDoppler=frequency*ve(0)/3e8;
	int chaNumofRealPointelp=NumOfSample;
	int fs=1500;
	while(NumofRealPoint)
	{
		if(NumofRealPoint<2*mDoppler*chaNumofRealPointelp/fs)
			NumofRealPoint=2*NumofRealPoint;
		else break;
	}
  int ifftnum=ceil_i(NumofRealPoint*fs/(2*mDoppler));//points of ifft
  double FreInterval=2*mDoppler;
  FreInterval=FreInterval/NumofRealPoint; //frequency interval after fft.
  double TimeInterval=1.0/fs;
  //Randomize the random number generator
  RNG_randomize();
  //Generates random Gaussian vector:
   cvec I_input_time=randn_c(NumofRealPoint);
   cvec Q_input_time=randn_c(NumofRealPoint);
   cvec I_input_freq=fft(I_input_time);//I q do fft transform
   cvec Q_input_freq=fft(Q_input_time);
   vec sez(NumofRealPoint+1);// fiter
   sez(1)=1.5/(pi*mDoppler);
   vec f(NumofRealPoint+1);
   f(1)=0;
   for (int j=2;j<=NumofRealPoint/2;j++)
   {
	   f(j)=(j-1)*FreInterval;
	   sez(j)=1.5/(pi*mDoppler*sqrt(1-pow(f(j)/mDoppler,2)));
	   sez(NumofRealPoint-j+2)=sez(j);
   }
   int NN=NumofRealPoint;
   int kk=3;
   vec result(T+1);//coeffienct of the result
   vec  x(4),y(4);
   int ii=NN/2-kk;
   int t=0;
   while (( ii>=(NN/2-kk))&&(ii<=(NN/2)))
   {
	   x(t)=ii;
	   y(t)=sez(ii);
	   t=t+1;
	   ii=ii+1;
   }
   result=process(x,y);
   sez(NN/2+1)=result(0)+result(1)*(NN/2+1)+result(2)*(NN/2+1)*(NN/2+1)+result(3)*(NN/2+1)*(NN/2+1)*(NN/2+1);
   cvec I_output_freq(NumofRealPoint);

   I_output_freq=ones_c(NumofRealPoint);
   vec sez2=concat(sez.right(NumofRealPoint-1),0.0);
   for(int k1=0;k1<NumofRealPoint;k1++)
   {
	   I_output_freq(k1)=sez2(k1)*I_input_freq(k1);
   }
   
   cvec Q_output_freq;
   Q_output_freq=ones_c(NumofRealPoint);
   vec sez3=concat(sez.right(NumofRealPoint-1),0.0);
   for(int k2=0;k2<NumofRealPoint;k2++)
   {
	   Q_output_freq(k2)=sez3(k2)*Q_input_freq(k2);
   }

   cvec I_temp(ifftnum);
   I_temp=concat(I_output_freq(0,(NumofRealPoint/2)-1),zeros_c(ifftnum-NumofRealPoint),I_output_freq(NumofRealPoint/2,NumofRealPoint-1));
   cvec I_output_time(ifftnum);
   I_output_time=ifft(I_temp,ifftnum);
   
   cvec Q_temp=concat(Q_output_freq(0,NumofRealPoint/2-1),zeros_c(ifftnum-NumofRealPoint),Q_output_freq(NumofRealPoint/2,NumofRealPoint-1));
   cvec Q_output_time=ifft(Q_temp);
   
   cvec vi(ifftnum);
   vi=to_cvec(0.0,1.0);
   cvec Tchannel(ifftnum);
   Tchannel=Q_output_time*vi(0);
   cvec channel(chaNumofRealPointelp);
   channel=Tchannel.left(chaNumofRealPointelp);
   return channel;
}

//--- Definition of RunScalarChannel()
cmat CVectorChannel::RunScalarChannel()
 {
     cmat scalarch(mNumOfMSAntenna*mNumOfBSAntenna,NumOfSample);
	 cvec ctemp;

	 for(int s1=0;s1<mNumOfMSAntenna*mNumOfBSAntenna;s1++)
	 {
		 ctemp=ScalarChannel(mVelocityOfMS,mFrequencyOfCarrier,mDoppler);
		 scalarch.set_row(s1,ctemp);
	 }
	 return scalarch;
 }

//--- Definition of TransmissionCoefficientPerPath()
cmat CVectorChannel::TransmissionCoefficientPerPath(int path)
{
	mat corrMS(mNumOfMSAntenna,mNumOfMSAntenna);
    mat corrBS(mNumOfBSAntenna,mNumOfBSAntenna);
	mat corrMB(mNumOfMSAntenna*mNumOfBSAntenna,mNumOfMSAntenna*mNumOfBSAntenna);
	mat cholMB(mNumOfMSAntenna*mNumOfBSAntenna,mNumOfMSAntenna*mNumOfBSAntenna);
	mat transMB(mNumOfMSAntenna*mNumOfBSAntenna,mNumOfMSAntenna*mNumOfBSAntenna);
	cmat coefficientp(mNumOfMSAntenna*mNumOfBSAntenna,NumOfSample);
	cmat ScalarC(mNumOfMSAntenna*mNumOfBSAntenna,NumOfSample);

    corrMS=CVectorChannel::ChannelCorrelation(mAngleSpreadPerPathMS(path),AOAPerPathMS,mAntennaSpacingMS,mNumOfMSAntenna);
    corrBS=CVectorChannel::ChannelCorrelation(mAngleSpreadPerPathBS(path),AODPerPathBS,mAntennaSpacingBS,mNumOfBSAntenna);
    corrMB=kron(corrBS,corrMS);
	cholMB=chol(corrMB);
	transMB=transpose(cholMB);
	ScalarC=RunScalarChannel();
	coefficientp=transMB*ScalarC;
 	return coefficientp;
}


cmat CVectorChannel::RunCVectorChannel()// To MIMO
{
  cmat *coefficient,*Tempcoefficient;
  coefficient=new cmat[mNumOfPath];
  Tempcoefficient=new cmat[mNumOfPath];
  for(int t=0;t<mNumOfPath;t++)
  {
    	Tempcoefficient[t].set_size(mNumOfBSAntenna*mNumOfMSAntenna,NumOfSample);
	    Tempcoefficient[t]=TransmissionCoefficientPerPath(t);
        coefficient[t].set_size(mNumOfBSAntenna,mNumOfMSAntenna*NumOfSample);
        coefficient[t]=reshape(Tempcoefficient[t],mNumOfBSAntenna,mNumOfMSAntenna*NumOfSample);
  }
  return *Tempcoefficient;
}