//----- PhyChannel.cpp -----
// Programmed by LTE Simulation Group of IMC, SWJTU
#include "stdafx.h"
#include <complex>
#include <iomanip>
#include <stdio.h>
#include "math.h"
#include "itpp/itcomm.h"
#include <itpp/itbase.h>
#include <itpp/base/mat.h> 
#include "PhyChannel.h"

using std::cin;
using std::cout;
using std::endl;
using namespace std;
using namespace itpp;

#define PI 3.14159265


//The functions class CPhyChannel

CPhyChannel::CPhyChannel()
{  
	
}

CPhyChannel::~CPhyChannel()
{
	
}




void CPhyChannel::GenScramblingSeq()     // Generate scrambling sequence
{

}

void CPhyChannel::Scrambling()           //Scrambling
{   
	int i,n_RNIT,q,ns,Nc,N_cell_id,Scram_length;
	double c_init;
	bvec Scram_input,Scram_output,x1,x2,Scram_sequence,temp;

    n_RNIT=2;             //Radio network temporary identifier
	q=0;
	ns=20;                //Slot number within a radio frame
	N_cell_id=0;
	Nc=60;
    Scram_input=randb(30);
    Scram_length=Scram_input.length();
	Scram_output.set_size(Scram_length,false);
    Scram_sequence.set_size(Scram_length,false);
	x1.set_size(Scram_length+Nc,false);
	x2.set_size(Scram_length+Nc,false);
	temp.set_size(1,false);
	
	//Pseudo-random sequence generation
	c_init=n_RNIT*pow(2.0,14)+q*pow(2.0,13)+(int)ceil(double(ns/2))*pow(2.0,q)+N_cell_id;
	cout<<"c_init="<<c_init<<endl;

    // The initialization of the first m-sequence 
	for(i=0;i<=30;i++)
	{
		if(i==0)
			x1(i)=1;
		else
			x1(i)=0;
	}
//	cout<<"x1="<<x1<<" "<<endl;


    // The initialization of the second m-sequence
	for(i=0;i<=30;i++)
	{
		x2(i)=(int)c_init%2;
		c_init=c_init/2;
	}
//	cout<<"x2="<<x2<<" "<<endl;   


	for(i=0;i<x1.length()-31;i++)
	{ 
	  x1(i+31)=mod((x1(i+3)+x1(i)),2);
      x2(i+31)=mod((x2(i+3)+x2(i+2)+x2(i+1)+x2(i)),2);
	}	
    
	cout<<"x1="<<x1<<endl;

	for(i=0;i<Scram_length;i++)
	{
	   Scram_sequence(i)=mod((x1(i+Nc)+x2(i+Nc)),2);
       Scram_output(i)=mod((Scram_input(i)+Scram_sequence(i)),2);
	}
	
	cout<<"Scram_input   ="<<Scram_input<<endl;
	cout<<"Scram_sequence="<<Scram_sequence<<endl;
	cout<<"Scram_output  ="<<Scram_output<<endl;
//	return Scram_output;
}

void CPhyChannel::Modulation()          //Modulation Mapper
{
    bvec inModuData;
	cvec outModuData;
	//PP str;
	
	inModuData=mMoPhyChInputData;
	//cout<<"调制输入信息=  "<<inModuData<<"调制输入**********"<<endl;

	if (mModulType==0)
	{		
		//ModuOfQPSK
		QPSK qpsk;                //The MDPSK modulator class
		outModuData = qpsk.modulate_bits(inModuData);
        //cout<<"调制输出信息= "<<outModuData<<"output*************"<<endl;

	}	    
	else if (mModulType==1)
	{
		//ModuOf16QAM
		QAM qam(16);
        outModuData= qam.modulate_bits(inModuData);
        //cout<<"Modulate symbols"<<outModuData<<endl;
	}
	else if (mModulType==2)
	{
		//ModuOf64QAM
		QAM qam(64);
        outModuData= qam.modulate_bits(inModuData);
        //cout<<"Modulate symbols"<<outModuData<<endl;
	}
    mTempDataSpace.set_size(1,outModuData.size(),false);
	//cout<<"rows="<<mTempDataSpace.rows()<<endl;
	for(int i=0;i<outModuData.size();i++)
	   mTempDataSpace(0,i)=outModuData(i);
	//***********************test
	mTempDataSpace2.set_size(1,outModuData.size(),false);
	for(int i=0;i<outModuData.size();i++)
	   mTempDataSpace2(0,i)=mTempDataSpace(0,i);
    //cout<<"outModuData"<<outModuData<<endl;
    //cout<<"rows="<<mTempDataSpace.rows()<<endl;
}

void CPhyChannel::MIMODec()           // MIMO detection
{
	typedef std::complex<double> complex;

    cmat inDecData="1.2846-0.3201i 1.2846-0.3201i;0.5046-0.1448i 0.5046-0.1448i;1.4653+2.1131i 1.4653+2.1131i;-2.0554-2.7233i -2.0554-2.7233i";                     // data from  receiver
	cmat outDecData;
	//cmat H=zeros(n,m);                           // matrix of channel 需要信道端给出 m transmitter,n receiver
	cmat H="-2.1707+0.4282i   0.5077+0.0403i   0.3803-0.3775i   0.0000+0.1184i;-0.0592+0.8956i   1.6924+0.6771i  -1.0091-0.2959i  -0.3179+0.3148i;-1.0106+0.7310i   0.5913+0.5689i  -0.0195-1.4751i   1.0950+1.4435i;0.6145+0.5779i  -0.6436-0.2556i  -0.0482-0.2340i  -1.8740-0.3510i";
    
	double sigma= 0.0040;                     // noise factor
	//PP str;
	//------- the transpose of H --------------------
	int p0;
	double a;
	int i ,j,m=4,n=4;
	cmat H0;        // m defind as the number of row,n defind as the number of column
	
	H0=transpose(H);
	
   //----------G1=inv(H_1'*H_1+((sigma/sqrt(10)).^2)*eye(2*Tx_n))*H_1'; ------------------------------
     cmat G;
	 G=inv(H0*H+(sigma*sigma)*eye(n))*H0;
	 mat G0=zeros(m,n);

	 //-------------------[gk1 m0]=min(sum(abs(G1).^2,2));-------------------------
	 mat  sum=zeros(m,1);
     for (i=0;i<m;i++)
		for (j=0;j<n;j++)
		{
			G0(i,j)=(abs(G(i,j))*abs(G(i,j))); 
	     
			/*cout<<G0(i,j)<<'\n';*/
	     sum(i)=sum(i)+G0(i,j);
	     }

	/* cout<<sum;*/
	a=sum(0);
	p0=0;
	for (i=1;i<m;i++)
	{if (sum(i)<a)
	    {a=sum(i);
	     p0=i;}
	}	
    //-----------------------------------------------------

	
      //int n;          // the column of H
	  mat m1(m,1);
	  cvec det_mmse,det_mmse2;
	  //double det,det1;
	  //double det2,det3;
	  int loop;
	  cmat result=zeros_c(m,inDecData.cols());         //根据具体情况改变 ‘1’
	  for (loop=0;loop<m;loop++)
	  {
		  m1(loop)=p0;
		  /*cout<<'\n'<<p0<<'\n';*/
		  //----------------------det_mmseosic=G1(m1(loop),:)*r_mmse;----
          det_mmse=G(m1(loop),m1(loop),0,n-1)*inDecData;
          /*cout<<G(m1(loop),m1(loop),0,n-1);*/
		  //if //16QAM modulation
		  QAM qam(16);
		  det_mmse2=qam.modulate_bits(qam.demodulate_bits(det_mmse));


		  //------------------r_mmse = r_mmse - det_mmseosic*H_2(:, m1(loop));----
		  cmat det_mmse3=zeros_c(n,inDecData.cols());            //根据具体情况改变  ‘1’
		  det_mmse3.set_row(0,det_mmse2);
		  //cout<<det_mmse3<<'\n'<<det_mmse3(0);
		  inDecData=inDecData-H(0,n-1,m1(loop),m1(loop))*det_mmse3(0,0,0,inDecData.cols()-1);
		  
		  //-------------result2(m1(loop),:)=det_mmseosic;-----------
		  result.set_row(m1(loop),det_mmse2);
		  //cout<<result<<'\n'<<"indecdata"<<inDecData<<det_mmse<<'\n';
		  //-------------H_2(:, m1(loop))=zeros(2*Rx_n,1);--------




		  H.set_submatrix(0,n-1,m1(loop),m1(loop),complex(0,0));
		  //-------------- G1=inv(H_2'*H_2+((sigma/sqrt(10)).^2)*eye(2*Tx_n))*H_2';
          H0=transpose(H);
          G=inv(H0*H+(sigma*sigma)*eye(n))*H0;
		  // -----------temp2=sum(abs(G1).^2,2)---------------------
		  mat  sum=zeros(m,1);
		  mat G0=zeros(m,n);
          for (i=0;i<m;i++)
		      for (j=0;j<n;j++)
		         {G0(i,j)=abs(G(i,j))*abs(G(i,j)); 
	              sum(i)=sum(i)+G0(i,j);
	             }
          //cout<<m1(loop)<<'\n'<<sum<<'\n';
		  //--------------temp2(m1(1:loop)) = 1e10;---------------------
		  int k=loop;
		  while (k>=0)
		      {sum(m1(k))=1e3;
		       k--;}
		  //------------------[gk1 m0]=min(temp2);------------------------------
		  //cout<<'\n'<<"sum"<<sum;
		  a=sum(0);
		  p0=0;
		  for (i=1;i<m;i++)
			  if(sum(i)<a)
			  {a=sum(i);
		       p0=i;}
	  }
	  outDecData=result;
	  //cout<<'\n'<<"outDecData="<<outDecData;
//	  return outDecData;
}


void CPhyChannel::Demodulation(double noisevar)      //Demodulation Mapper
{
	int i;
	cvec inDemoduData;
	vec inDemoduData1;
	vec outDemoduData;
	//inDemoduData=outDecData;
    //PP str;
    inDemoduData.set_size(mTempDataSpace.cols(),false);
	inDemoduData1.set_size(mTempDataSpace.cols(),false);
	for(i=0;i<inDemoduData.length();i++)
	    inDemoduData(i)=mTempDataSpace(0,i);
	//cout<<"ffffffff"<<endl;

	if (mModulType==0)
	{
		//DemoduOfQPSK
		QPSK qpsk;                //The MDPSK modulator class
		outDemoduData = qpsk.demodulate_soft_bits(inDemoduData,noisevar);
		 //demodulate_bits(inDemoduData);
		//cout<<"demodulate bits"<<outDemoduData<<"&&&&&&&&"<<endl;
	}
	else if (mModulType==1)
	{
		//DemoduOf16QAM
		QAM qam(16);
		outDemoduData=qam.demodulate_soft_bits(inDemoduData,noisevar);
			//(inDemoduData,noisevar);
        //cout<<"demodulate bits"<<outDemoduData<<endl;
	}
	else if (mModulType==2)
	{
		//DemoduOf64QAM
        QAM qam(64);
		outDemoduData=qam.demodulate_soft_bits(inDemoduData,noisevar);
			//demodulate_bits(inDemoduData,noisevar);
        //cout<<"demodulate bits"<<outDemoduData<<endl;
	}
  
	mTempDataSpace1.set_size(1,outDemoduData.length(),false);
	for(i=0;i<outDemoduData.length();i++)
	    mTempDataSpace1(0,i)=outDemoduData(i);

}

void CPhyChannel::DeScrambling()      //Inverse of Scrambling
{
//	return 0;
}


//The functions class CPUSch


CPUSCh::CPUSCh()
{
  mAmpScaFactor=0; 
}
CPUSCh::~CPUSCh()
{

}


void CPUSCh::InniCPhyChannel(LCP para,GP para1)
 {
    mUpDownlink=para.mUpDownlink;
	mPRBConfig=para.mDownCPType;
	mModulType=para.mDownModType;
	//mReduVersiNum=para.mReduVersiNum;
	mSoftbuffSize=para.mSoftBuffSize;
	mRBNum=para1.mNumofUpRB;
	mCarrierType=para1.mSubcarrierType;
	//cout<< "CPUSCh:mUpDownlink ="<<mUpDownlink<<endl;

	mULSCh.InniCTranChannel(mSoftbuffSize);

 }



void CPUSCh::Run()
{
  bvec ScrambleData;
  cvec ModulData;
  cmat LayerMapData,PrecodeData,REmapperData;
  mULSCh.Run();
  mMoPhyChInputData=mULSCh.mpTrChOutputData; //从传输信道得到数据
  Scrambling();
  //ModulData = Modulation();
//  PrecodeData = Precoding();
  REmapper(); 

  mMoPhyChOutputData=REmapperData;

}
void CPUSCh::DRun()
{  
	bvec ErrorNum;
  bvec DeScrambleData;
  bvec DeModulData;
  cmat DeLayerMapData,DePrecodeData,DeREmapperData;
  //mPhyChInputData=GenScramblingSeq();  //从接收机接受数据
  DeREmapperData=mDePhyChInputData; 
  DePrecoding();
  //DeModulData = Demodulation ();
   DeScrambling();
  mDePhyChOutputData=DeScrambleData;

}

void CPUSCh::MultiDataControl()         //Data and Control multiplexing
{
//	return 0;
}

void CPUSCh::ChannelInterl()            //Channel Interleaver
{
//	return 0;
}

void CPUSCh::Precoding()      // Precoding
{
//	mPrecodedSymbolUp=fft(mMoPhyChOutputData);
//	return 0;
}	

void CPUSCh::DeMultiDataControl()      //Demultiplexing 
{
//	return 0;
}
	 
void CPUSCh::DeChannelInterl()         //Channel Deinterleaver
{
//	return 0;
}
	
    
void  CPUSCh::DePrecoding()    // Inverse of precoding for PUSCH: FFT
{
//	return 0;
}
   

void CPUSCh::REmapper()
{
    double scale_factor=0.34;
	cvec Z;
	cmat out;
   
	int i,j,k,l,Nsc,Nslot=1;
	int num=mRBNum;
	int symb=7*Nslot; 
    int m=length(Z);

	if(mCarrierType==0) Nsc=12;
	else Nsc=24;

	for(l=0;l<3;l++)
	{
	    for(k=0;k<num*Nsc;k++)
	    {
		   out(k,l)=scale_factor*Z(k+l*num*Nsc);
	    }
	}

    for(l=4;l<symb;l++)
	{
	    for(k=0;k<num*Nsc;k++)
	    {
		   out(k,l)=scale_factor*Z(k+(l-1)*num*Nsc);
	    }
	}
	
   	mTempDataSpace.set_size(out.rows(),out.cols(),false);
	for(i=0;i<mTempDataSpace.rows();i++)
       for(j=0;j<mTempDataSpace.cols();j++)
	     mTempDataSpace(i,j)=out(i,j);

//	return out;
}

void CPUSCh::DeREmapper()
{
	cmat intput;
	cvec out;
   
	int j,k,l,Nsc,Nslot=8;//8?
	int num=mRBNum;
	int symb=7*Nslot; 

	if(mCarrierType==0) Nsc=12;
	else Nsc=24;

	for(l=0;l<3;l++)
	{
	    for(k=0;k<num*Nsc;k++)
	    {
		   out(k+l*num*Nsc)=intput(k,l);
	    }
	}

    for(l=4;l<symb;l++)
	{
	    for(k=0;k<num*Nsc;k++)
	    {
		   out(k+(l-1)*num*Nsc)=intput(k,l);;
	    }
	}
		
	mTempDataSpace.set_size(1,out.length(),false);	
       for(j=0;j<mTempDataSpace.cols();j++)
	     mTempDataSpace(0,j)=out(j);
//return out;
}



//The functions class CPDSch


CPDSCh::CPDSCh()
{

}
CPDSCh::~CPDSCh()
{

}


void CPDSCh::InniCPhyChannel(LCP para,GP para1)
 {
    mUpDownlink=para.mUpDownlink;
	mPRBConfig=para.mDownCPType;
	mModulType=para.mDownModType;
	//mReduVersiNum=para.mReduVersiNum;
	mSoftbuffSize=para.mSoftBuffSize;
	mStreamNum=para.mStreamNum;
	mLayerNum=para.mLayerNum;
	mRBNum=para1.mNumofDwRB;
	RxAntennaNum=para1.mNumofRxAntenna;
	TxAntennaNum=para1.mNumofTxAntenna;
	mCarrierType=para1.mSubcarrierType;
	//cout<< "CPDSCh:mUpDownlink ="<<mUpDownlink<<endl;

	mDLSCh.InniCTranChannel(mSoftbuffSize);

 }



void CPDSCh::Precoding()
{
	int i,j,k;
	int PathNum = 6;//径数*****
	double Segma = 0.04;//噪声功率 *****
	TxAntennaNum = 4;//发射天线数*****
	RxAntennaNum =2;//接收天线数*****
	mLayerNum = 2;//层数***
	int CodeBookSize;
	int mSubBandWidth = 50;//由带宽决定，mmSubBandWidth*SubBandNum就是用的子载波的数量*****
	int mSubBandNum = 1;
	RBWidth = 12;
	int ModulatedDataLen=mMoPhyChOutputData.cols();
	int LayerSymbolNum=ModulatedDataLen/mLayerNum;//每层符号数
	cmat LayerMappedSymbol(mLayerNum,LayerSymbolNum);
	cmat PrecodedSymbol(TxAntennaNum,mmSubBandWidth*RBWidth);
	

	cmat ChannelTimeDomain[1024];//只考虑一个服务扇区，信道****
	cmat ChannelFreqDomain[1024];//频域信道***

	cmat UsedChannelFreqDomain[1024];//使用1-600的子载波
	cmat RBChannelFreqDomain[10];//10个subband，假设一个用户使用所有subband
	cmat CodeBook[16];//预编码矩阵
	ivec CodeBookIndex(10);
	cmat UsedCodeBook[10];//10个subband最后选择的码本
	cmat MSEMatrix[11][16];//每个subband每个码本的判决矩阵，16是可能的最大codebook size
	mat Det(11,16);//每个判决矩阵的行列式,不是复数,16是可能的最大codebook size
	

	/**************************************************************************/

	if (TxAntennaNum==2)//2天线码本
	{
		if(mLayerNum ==1)
		{
			CodeBookSize = 4;
			vec RealValue00 = "0.707107 0.707107 0.707107 0.707107";
			vec RealValue10 = "0.707107 -0.707107 0.000000 0.000000";
			vec ImagValue00 = "0.000000 0.000000 0.000000 0.000000";
			vec ImagValue10 = "0.000000 0.000000 0.707107 -0.707107";
			for(i=0;i<CodeBookSize;i++)
			{
				CodeBook[i].set_size(2,1,0);
			}
			for(i=0;i<CodeBookSize;i++)
			{
				//column 1
				CodeBook[i](0,0).real(RealValue00(i));
				CodeBook[i](1,0).real(RealValue10(i));	
				CodeBook[i](0,0).imag(ImagValue00(i));
				CodeBook[i](1,0).imag(ImagValue10(i));
			}
		}
		if(mLayerNum == 2)
		{
			CodeBookSize = 2;
			vec RealValue00 = "0.500000 0.500000 ";
			vec RealValue10 = "0.500000 0.000000 ";
			vec ImagValue00 = "0.000000 0.000000 ";
			vec ImagValue10 = "0.000000 0.500000 ";
			vec RealValue01 = "0.500000 0.500000 ";
			vec RealValue11 = "-0.500000 0.000000 ";
			vec ImagValue01 = "0.000000 0.000000 ";
			vec ImagValue11 = "0.000000 -0.500000 ";
			for(i=0;i<CodeBookSize;i++)
			{
				CodeBook[i].set_size(2,2,0);
			}
			for(i=0;i<CodeBookSize;i++)
			{
				//column 1
				CodeBook[i](0,0).real(RealValue00(i));
				CodeBook[i](1,0).real(RealValue10(i));	
				CodeBook[i](0,0).imag(ImagValue00(i));
				CodeBook[i](1,0).imag(ImagValue10(i));
				//column 2
				CodeBook[i](0,1).real(RealValue01(i));
				CodeBook[i](1,1).real(RealValue11(i));	
				CodeBook[i](0,1).imag(ImagValue01(i));
				CodeBook[i](1,1).imag(ImagValue11(i));
			}
		}
	}

	if(TxAntennaNum==4)//4天线码本
	{
		CodeBookSize = 16;
		if(mLayerNum==1)//rank1,4行一列
		{
			vec RealValue00 = "0.500000 0.500000 0.500000 0.500000 0.500000 0.500000 0.500000 0.500000 0.500000 0.500000 0.500000 0.500000 0.500000 0.500000 0.500000 0.500000 ";//实部的值
			vec RealValue10 = "0.500000 0.000000 -0.500000 0.000000 0.353553 -0.353553 -0.353553 0.353553 0.500000 0.000000 -0.500000 0.000000 0.500000 0.500000 -0.500000 -0.500000 ";
			vec RealValue20 = "0.500000 0.500000 0.500000 -0.500000 0.000000 0.000000 0.000000 0.000000 -0.500000 0.500000 -0.500000 0.500000 0.500000 -0.500000 0.500000 -0.500000 ";
			vec RealValue30 = "0.500000 0.000000 -0.500000 0.000000 -0.353553 0.353553 0.353553 -0.353553 -0.500000 0.000000 0.500000 0.000000 -0.500000 0.500000 0.500000 -0.500000 ";
			vec ImagValue00 = "0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 ";//虚部的值
			vec ImagValue10 = "0.000000 0.500000 0.000000 -0.500000 0.353553 0.353553 -0.353553 -0.353553 0.000000 0.500000 0.000000 -0.500000 0.000000 0.000000 0.000000 0.000000 ";
			vec ImagValue20 = "0.000000 0.000000 0.000000 0.000000 0.500000 -0.500000 0.500000 -0.500000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 ";
			vec ImagValue30 = "0.000000 -0.500000 0.000000 0.500000 0.353553 0.353553 -0.353553 -0.353553 0.000000 0.500000 0.000000 -0.500000 0.000000 0.000000 0.000000 0.000000 ";	
		
			for(i=0;i<CodeBookSize;i++)
			{
				CodeBook[i].set_size(4,1,0);
			}
			for(i=0;i<CodeBookSize;i++)
			{
				//column 1
				CodeBook[i](0,0).real(RealValue00(i));
				CodeBook[i](1,0).real(RealValue10(i));
				CodeBook[i](2,0).real(RealValue20(i));
				CodeBook[i](3,0).real(RealValue30(i));
				CodeBook[i](0,0).imag(ImagValue00(i));
				CodeBook[i](1,0).imag(ImagValue10(i));
				CodeBook[i](2,0).imag(ImagValue20(i));
				CodeBook[i](3,0).imag(ImagValue30(i));
			}
		}
		if(mLayerNum == 2)
		{
			//column 1
			vec RealValue00 ="0.353553 0.353553 0.353553 0.353553 0.353553 0.353553 0.353553 0.353553 0.353553 0.353553 0.353553 0.353553 0.353553 0.353553 0.353553 0.353553 ";
			vec RealValue10 ="0.353553 0.000000 -0.353553 0.000000 0.250000 -0.250000 -0.250000 0.250000 0.353553 0.000000 -0.353553 0.000000 0.353553 0.353553 -0.353553 -0.353553 ";
			vec RealValue20 ="0.353553 0.353553 0.353553 -0.353553 0.000000 0.000000 0.000000 0.000000 -0.353553 0.353553 -0.353553 0.353553 0.353553 -0.353553 0.353553 -0.353553 ";
			vec RealValue30 ="0.353553 0.000000 -0.353553 0.000000 -0.250000 0.250000 0.250000 -0.250000 -0.353553 0.000000 0.353553 0.000000 -0.353553 0.353553 0.353553 -0.353553 ";
			vec ImagValue00 ="0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 ";
			vec ImagValue10 ="0.000000 0.353553 0.000000 -0.353553 0.250000 0.250000 -0.250000 -0.250000 0.000000 0.353553 0.000000 -0.353553 0.000000 0.000000 0.000000 0.000000 ";
			vec ImagValue20 ="0.000000 0.000000 0.000000 0.000000 0.353553 -0.353553 0.353553 -0.353553 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 ";
			vec ImagValue30 ="0.000000 -0.353553 0.000000 0.353553 0.250000 0.250000 -0.250000 -0.250000 0.000000 0.353553 0.000000 -0.353553 0.000000 0.000000 0.000000 0.000000 ";
			//column 2
			vec RealValue01 ="0.353553 0.000000 -0.353553 0.000000 -0.250000 0.250000 0.000000 0.000000 0.353553 0.000000 -0.353553 0.353553 0.353553 -0.353553 0.353553 -0.353553 ";
			vec RealValue11 ="-0.353553 0.353553 0.353553 0.353553 0.000000 0.000000 0.250000 -0.250000 0.353553 -0.353553 -0.353553 0.000000 0.353553 0.353553 0.353553 0.353553 ";
			vec RealValue21 ="-0.353553 0.000000 0.353553 0.000000 -0.250000 0.250000 0.353553 0.353553 0.353553 0.000000 0.353553 0.353553 -0.353553 0.353553 0.353553 -0.353553 ";
			vec RealValue31 ="0.353553 0.353553 -0.353553 0.353553 0.353553 0.353553 0.250000 -0.250000 0.353553 0.353553 0.353553 0.000000 0.353553 0.353553 -0.353553 -0.353553 ";
			vec ImagValue01 ="0.000000 -0.353553 0.000000 0.353553 -0.250000 -0.250000 -0.353553 0.353553 0.000000 -0.353553 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 ";
			vec ImagValue11 ="0.000000 0.000000 0.000000 0.000000 0.353553 -0.353553 -0.250000 -0.250000 0.000000 0.000000 0.000000 0.353553 0.000000 0.000000 0.000000 0.000000 ";
			vec ImagValue21 ="0.000000 0.353553 0.000000 0.353553 0.250000 0.250000 0.000000 0.000000 0.000000 0.353553 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 ";
			vec ImagValue31 ="0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.250000 0.250000 0.000000 0.000000 0.000000 0.353553 0.000000 0.000000 0.000000 0.000000 ";

			for(i=0;i<CodeBookSize;i++)
			{
				CodeBook[i].set_size(4,2,0);
			}
			for(i=0;i<CodeBookSize;i++)
			{
				//column 1
				CodeBook[i](0,0).real(RealValue00(i));
				CodeBook[i](1,0).real(RealValue10(i));
				CodeBook[i](2,0).real(RealValue20(i));
				CodeBook[i](3,0).real(RealValue30(i));
				CodeBook[i](0,0).imag(ImagValue00(i));
				CodeBook[i](1,0).imag(ImagValue10(i));
				CodeBook[i](2,0).imag(ImagValue20(i));
				CodeBook[i](3,0).imag(ImagValue30(i));
				//column 2
				CodeBook[i](0,1).real(RealValue01(i));
				CodeBook[i](1,1).real(RealValue11(i));
				CodeBook[i](2,1).real(RealValue21(i));
				CodeBook[i](3,1).real(RealValue31(i));
				CodeBook[i](0,1).imag(ImagValue01(i));
				CodeBook[i](1,1).imag(ImagValue11(i));
				CodeBook[i](2,1).imag(ImagValue21(i));
				CodeBook[i](3,1).imag(ImagValue31(i));
			}
			
		}
		if(mLayerNum == 3)
		{
			//column 1
			vec RealValue00 ="0.353553 0.353553 0.353553 0.353553 0.353553 0.353553 0.353553 0.353553 0.353553 0.353553 0.353553 0.353553 0.353553 0.353553 0.353553 0.353553 ";
			vec RealValue10 ="0.353553 0.000000 -0.353553 0.000000 0.250000 -0.250000 -0.250000 0.250000 0.353553 0.000000 -0.353553 0.000000 0.353553 0.353553 -0.353553 -0.353553 ";
			vec RealValue20 ="0.353553 0.353553 0.353553 -0.353553 0.000000 0.000000 0.000000 0.000000 -0.353553 0.353553 -0.353553 0.353553 0.353553 -0.353553 0.353553 -0.353553 ";
			vec RealValue30 ="0.353553 0.000000 -0.353553 0.000000 -0.250000 0.250000 0.250000 -0.250000 -0.353553 0.000000 0.353553 0.000000 -0.353553 0.353553 0.353553 -0.353553 ";
			vec ImagValue00 ="0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 ";
			vec ImagValue10 ="0.000000 0.353553 0.000000 -0.353553 0.250000 0.250000 -0.250000 -0.250000 0.000000 0.353553 0.000000 -0.353553 0.000000 0.000000 0.000000 0.000000 ";
			vec ImagValue20 ="0.000000 0.000000 0.000000 0.000000 0.353553 -0.353553 0.353553 -0.353553 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 ";
			vec ImagValue30 ="0.000000 -0.353553 0.000000 0.353553 0.250000 0.250000 -0.250000 -0.250000 0.000000 0.353553 0.000000 -0.353553 0.000000 0.000000 0.000000 0.000000 ";
			//column 2
			vec RealValue01 ="0.353553 0.000000 -0.353553 0.000000 0.250000 -0.250000 0.000000 0.000000 0.353553 0.353553 -0.353553 0.353553 0.353553 0.353553 -0.353553 -0.353553 ";
			vec RealValue11 ="0.353553 0.353553 0.353553 0.353553 0.353553 0.353553 0.250000 -0.250000 0.353553 0.000000 0.353553 0.000000 0.353553 0.353553 0.353553 0.353553 ";
			vec RealValue21 ="-0.353553 0.000000 0.353553 0.000000 -0.250000 0.250000 0.353553 0.353553 0.353553 0.353553 -0.353553 0.353553 -0.353553 0.353553 0.353553 -0.353553 ";
			vec RealValue31 ="-0.353553 0.353553 -0.353553 0.353553 0.000000 0.000000 0.250000 -0.250000 0.353553 0.000000 0.353553 0.000000 0.353553 -0.353553 0.353553 -0.353553 ";
			vec ImagValue01 ="0.000000 -0.353553 0.000000 0.353553 -0.250000 -0.250000 -0.353553 0.353553 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 ";
			vec ImagValue11 ="0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 -0.250000 -0.250000 0.000000 -0.353553 0.000000 0.353553 0.000000 0.000000 0.000000 0.000000 ";
			vec ImagValue21 ="0.000000 0.353553 0.000000 0.353553 -0.250000 -0.250000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 ";
			vec ImagValue31 ="0.000000 0.000000 0.000000 0.000000 -0.353553 0.353553 0.250000 0.250000 0.000000 -0.353553 0.000000 0.353553 0.000000 0.000000 0.000000 0.000000 ";
			//column 3
			vec RealValue02 ="0.353553 0.353553 0.353553 -0.353553 -0.250000 0.250000 0.250000 -0.250000 -0.353553 0.000000 -0.353553 0.000000 0.353553 -0.353553 0.353553 -0.353553 ";
			vec RealValue12 ="-0.353553 0.000000 0.353553 0.000000 0.000000 0.000000 0.000000 0.000000 0.353553 -0.353553 -0.353553 -0.353553 -0.353553 0.353553 0.353553 -0.353553 ";
			vec RealValue22 ="-0.353553 0.353553 0.353553 0.353553 -0.250000 0.250000 0.250000 -0.250000 -0.353553 0.000000 0.353553 0.000000 0.353553 0.353553 0.353553 0.353553 ";
			vec RealValue32 ="0.353553 0.000000 0.353553 0.000000 0.353553 0.353553 0.353553 0.353553 0.353553 0.353553 0.353553 0.353553 0.353553 0.353553 -0.353553 -0.353553 ";
			vec ImagValue02 ="0.000000 0.000000 0.000000 0.000000 -0.250000 -0.250000 0.250000 0.250000 0.000000 -0.353553 0.000000 0.353553 0.000000 0.000000 0.000000 0.000000 ";
			vec ImagValue12 ="0.000000 -0.353553 0.000000 -0.353553 0.353553 -0.353553 0.353553 -0.353553 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 ";
			vec ImagValue22 ="0.000000 0.000000 0.000000 0.000000 0.250000 0.250000 -0.250000 -0.250000 0.000000 0.353553 0.000000 -0.353553 0.000000 0.000000 0.000000 0.000000 ";
			vec ImagValue32 ="0.000000 0.353553 0.000000 0.353553 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 ";

			for(i=0;i<CodeBookSize;i++)
			{
				CodeBook[i].set_size(4,3,0);
			}
			for(i=0;i<CodeBookSize;i++)
			{
				//column 1
				CodeBook[i](0,0).real(RealValue00(i));
				CodeBook[i](1,0).real(RealValue10(i));
				CodeBook[i](2,0).real(RealValue20(i));
				CodeBook[i](3,0).real(RealValue30(i));
				CodeBook[i](0,0).imag(ImagValue00(i));
				CodeBook[i](1,0).imag(ImagValue10(i));
				CodeBook[i](2,0).imag(ImagValue20(i));
				CodeBook[i](3,0).imag(ImagValue30(i));
				//column 2
				CodeBook[i](0,1).real(RealValue01(i));
				CodeBook[i](1,1).real(RealValue11(i));
				CodeBook[i](2,1).real(RealValue21(i));
				CodeBook[i](3,1).real(RealValue31(i));
				CodeBook[i](0,1).imag(ImagValue01(i));
				CodeBook[i](1,1).imag(ImagValue11(i));
				CodeBook[i](2,1).imag(ImagValue21(i));
				CodeBook[i](3,1).imag(ImagValue31(i));
				//column 2
				CodeBook[i](0,2).real(RealValue02(i));
				CodeBook[i](1,2).real(RealValue12(i));
				CodeBook[i](2,2).real(RealValue22(i));
				CodeBook[i](3,2).real(RealValue32(i));
				CodeBook[i](0,2).imag(ImagValue02(i));
				CodeBook[i](1,2).imag(ImagValue12(i));
				CodeBook[i](2,2).imag(ImagValue22(i));
				CodeBook[i](3,2).imag(ImagValue32(i));
			}
		}
		if(mLayerNum == 4)
		{
			//column 1
			vec RealValue00 ="0.353553 0.353553 0.353553 -0.353553 0.353553 0.353553 0.353553 0.353553 0.353553 0.353553 0.353553 0.353553 0.353553 0.353553 0.353553 0.353553 ";
			vec RealValue10 ="0.353553 0.000000 0.353553 0.000000 0.250000 -0.250000 -0.250000 0.250000 0.353553 0.000000 -0.353553 0.000000 0.353553 0.353553 0.353553 -0.353553 ";
			vec RealValue20 ="0.353553 0.353553 0.353553 0.353553 0.000000 0.000000 0.000000 0.000000 -0.353553 0.353553 -0.353553 0.353553 0.353553 -0.353553 0.353553 -0.353553 ";
			vec RealValue30 ="0.353553 0.000000 0.353553 0.000000 -0.250000 0.250000 0.250000 -0.250000 -0.353553 0.000000 0.353553 0.000000 -0.353553 0.353553 -0.353553 -0.353553 ";
			vec ImagValue00 ="0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 ";
			vec ImagValue10 ="0.000000 0.353553 0.000000 -0.353553 0.250000 0.250000 -0.250000 -0.250000 0.000000 0.353553 0.000000 -0.353553 0.000000 0.000000 0.000000 0.000000 ";
			vec ImagValue20 ="0.000000 0.000000 0.000000 0.000000 0.353553 -0.353553 0.353553 -0.353553 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 ";
			vec ImagValue30 ="0.000000 -0.353553 0.000000 0.353553 0.250000 0.250000 -0.250000 -0.250000 0.000000 0.353553 0.000000 -0.353553 0.000000 0.000000 0.000000 0.000000 ";
			//column 2
			vec RealValue01 ="0.353553 0.000000 -0.353553 0.000000 0.250000 -0.250000 0.000000 0.000000 0.353553 0.000000 -0.353553 0.353553 0.353553 -0.353553 -0.353553 -0.353553 ";
			vec RealValue11 ="0.353553 0.353553 0.353553 0.353553 0.353553 0.353553 0.250000 -0.250000 0.353553 0.353553 -0.353553 0.000000 0.353553 0.353553 0.353553 0.353553 ";
			vec RealValue21 ="-0.353553 0.000000 0.353553 0.000000 -0.250000 0.250000 0.353553 0.353553 0.353553 0.000000 0.353553 0.353553 -0.353553 0.353553 0.353553 -0.353553 ";
			vec RealValue31 ="-0.353553 0.353553 -0.353553 0.353553 0.000000 0.000000 0.250000 -0.250000 0.353553 -0.353553 0.353553 0.000000 0.353553 0.353553 0.353553 -0.353553 ";
			vec ImagValue01 ="0.000000 -0.353553 0.000000 0.353553 -0.250000 -0.250000 -0.353553 0.353553 0.000000 -0.353553 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 ";
			vec ImagValue11 ="0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 -0.250000 -0.250000 0.000000 0.000000 0.000000 0.353553 0.000000 0.000000 0.000000 0.000000 ";
			vec ImagValue21 ="0.000000 0.353553 0.000000 0.353553 -0.250000 -0.250000 0.000000 0.000000 0.000000 0.353553 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 ";
			vec ImagValue31 ="0.000000 0.000000 0.000000 0.000000 -0.353553 0.353553 0.250000 0.250000 0.000000 0.000000 0.000000 0.353553 0.000000 0.000000 0.000000 0.000000 ";
			//column 
			vec RealValue02 ="0.353553 0.353553 0.353553 0.353553 0.000000 0.000000 -0.250000 0.250000 -0.353553 0.353553 -0.353553 0.000000 0.353553 0.353553 0.353553 -0.353553 ";
			vec RealValue12 ="-0.353553 0.000000 -0.353553 0.000000 -0.250000 0.250000 0.353553 0.353553 0.353553 0.000000 0.353553 0.353553 -0.353553 0.353553 -0.353553 -0.353553 ";
			vec RealValue22 ="0.353553 0.353553 0.353553 -0.353553 0.353553 0.353553 0.250000 -0.250000 0.353553 0.353553 -0.353553 0.000000 0.353553 0.353553 0.353553 0.353553 ";
			vec RealValue32 ="-0.353553 0.000000 -0.353553 0.000000 -0.250000 0.250000 0.000000 0.000000 -0.353553 0.000000 0.353553 -0.353553 0.353553 -0.353553 0.353553 -0.353553 ";
			vec ImagValue02 ="0.000000 0.000000 0.000000 0.000000 -0.353553 0.353553 0.250000 0.250000 0.000000 0.000000 0.000000 0.353553 0.000000 0.000000 0.000000 0.000000 ";
			vec ImagValue12 ="0.000000 -0.353553 0.000000 -0.353553 0.250000 0.250000 0.000000 0.000000 0.000000 -0.353553 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 ";
			vec ImagValue22 ="0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.250000 0.250000 0.000000 0.000000 0.000000 -0.353553 0.000000 0.000000 0.000000 0.000000 ";
			vec ImagValue32 ="0.000000 0.353553 0.000000 0.353553 -0.250000 -0.250000 -0.353553 0.353553 0.000000 -0.353553 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 ";
			//column 3
			vec RealValue03 ="0.353553 0.000000 -0.353553 0.000000 -0.250000 0.250000 0.250000 -0.250000 -0.353553 0.000000 0.353553 0.000000 -0.353553 0.353553 0.353553 -0.353553 ";
			vec RealValue13 ="-0.353553 0.353553 -0.353553 0.353553 0.000000 0.000000 0.000000 0.000000 0.353553 -0.353553 0.353553 -0.353553 0.353553 -0.353553 0.353553 -0.353553 ";
			vec RealValue23 ="-0.353553 0.000000 0.353553 0.000000 -0.250000 0.250000 0.250000 -0.250000 -0.353553 0.000000 0.353553 0.000000 0.353553 0.353553 -0.353553 -0.353553 ";
			vec RealValue33 ="0.353553 0.353553 0.353553 0.353553 0.353553 0.353553 0.353553 0.353553 0.353553 0.353553 0.353553 0.353553 0.353553 0.353553 0.353553 0.353553 ";
			vec ImagValue03 ="0.000000 0.353553 0.000000 -0.353553 -0.250000 -0.250000 0.250000 0.250000 0.000000 -0.353553 0.000000 0.353553 0.000000 0.000000 0.000000 0.000000 ";
			vec ImagValue13 ="0.000000 0.000000 0.000000 0.000000 0.353553 -0.353553 0.353553 -0.353553 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 ";
			vec ImagValue23 ="0.000000 -0.353553 0.000000 -0.353553 0.250000 0.250000 -0.250000 -0.250000 0.000000 0.353553 0.000000 -0.353553 0.000000 0.000000 0.000000 0.000000 ";
			vec ImagValue33 ="0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 ";
			
			for(i=0;i<CodeBookSize;i++)
			{
				CodeBook[i].set_size(4,4,0);
			}
			for(i=0;i<CodeBookSize;i++)
			{
				//column 1
				CodeBook[i](0,0).real(RealValue00(i));
				CodeBook[i](1,0).real(RealValue10(i));
				CodeBook[i](2,0).real(RealValue20(i));
				CodeBook[i](3,0).real(RealValue30(i));
				CodeBook[i](0,0).imag(ImagValue00(i));
				CodeBook[i](1,0).imag(ImagValue10(i));
				CodeBook[i](2,0).imag(ImagValue20(i));
				CodeBook[i](3,0).imag(ImagValue30(i));
				//column 2
				CodeBook[i](0,1).real(RealValue01(i));
				CodeBook[i](1,1).real(RealValue11(i));
				CodeBook[i](2,1).real(RealValue21(i));
				CodeBook[i](3,1).real(RealValue31(i));
				CodeBook[i](0,1).imag(ImagValue01(i));
				CodeBook[i](1,1).imag(ImagValue11(i));
				CodeBook[i](2,1).imag(ImagValue21(i));
				CodeBook[i](3,1).imag(ImagValue31(i));
				//column 3
				CodeBook[i](0,2).real(RealValue02(i));
				CodeBook[i](1,2).real(RealValue12(i));
				CodeBook[i](2,2).real(RealValue22(i));
				CodeBook[i](3,2).real(RealValue32(i));
				CodeBook[i](0,2).imag(ImagValue02(i));
				CodeBook[i](1,2).imag(ImagValue12(i));
				CodeBook[i](2,2).imag(ImagValue22(i));
				CodeBook[i](3,2).imag(ImagValue32(i));
				//column 4
				CodeBook[i](0,3).real(RealValue03(i));
				CodeBook[i](1,3).real(RealValue13(i));
				CodeBook[i](2,3).real(RealValue23(i));
				CodeBook[i](3,3).real(RealValue33(i));
				CodeBook[i](0,3).imag(ImagValue03(i));
				CodeBook[i](1,3).imag(ImagValue13(i));
				CodeBook[i](2,3).imag(ImagValue23(i));
				CodeBook[i](3,3).imag(ImagValue33(i));
			}			
		}
	}
	//以上为码本产生
	//初始化信道向量
	for(i=0;i<1024;i++)
	{
		ChannelTimeDomain[i].set_size(RxAntennaNum,TxAntennaNum);
		ChannelTimeDomain[i].zeros();//置零
	}
	for(i=0;i<1024;i++)//后面的0都没用
	{
		ChannelFreqDomain[i].set_size(RxAntennaNum,TxAntennaNum);
		ChannelFreqDomain[i].zeros();//置零
	}
	//产生6径信道
	for(i=0;i<PathNum;i++)
	{
		ChannelTimeDomain[i] = randn_c(RxAntennaNum,TxAntennaNum);//矩阵赋给矩阵
	}	
	//信道fft
	for(i=0;i<RxAntennaNum;i++)	
	{
		for(j=0;j<TxAntennaNum;j++)
		{
			cvec Temp1(1024);
			cvec Temp2(1024);
			for(k=0;k<1024;k++)
				Temp1(k) = ChannelTimeDomain[k](i,j);
			Temp2 = fft(Temp1,1024);
			for(k=0;k<1024;k++)
				ChannelFreqDomain[k](i,j) = Temp2(k);
		}
	}
	//取出使用的子载波,使用1-600子载波，共600个，10个subband,每个子带5个RB
	for(j=0,i=1;i<mSubBandNum*mmSubBandWidth*RBWidth+1;i++,j++)
	{
		UsedChannelFreqDomain[j] = ChannelFreqDomain[i];
	}
	//信道平均到每个子带上？？？
	for(i=0;i<mSubBandNum;i++)
	{
		cmat sum;
		sum.set_size(RxAntennaNum,TxAntennaNum,0);
		sum.zeros();
		for(j=0;j<RBWidth*mmSubBandWidth;j++)
		{
			sum = sum + ChannelFreqDomain[i*RBWidth*mmSubBandWidth + j];
		}
		RBChannelFreqDomain[i] = sum/(RBWidth*mmSubBandWidth);
	}
	//预编码矩阵决定
	for(i=0;i<mSubBandNum;i++)//10个子带
	{
		cmat ChannelTemp;//等效信道矩阵
		for(j=0;j<CodeBookSize;j++)
		{		
			ChannelTemp = RBChannelFreqDomain[i] * CodeBook[j];
			MSEMatrix[i][j] = itpp::inv(itpp::hermitian_transpose(ChannelTemp)*ChannelTemp + itpp::eye_c(mLayerNum));//公式不完整，还有发射功率和噪声功率
			Det(i,j) = itpp::det(MSEMatrix[i][j]).real();//虚部=0
		}
		//当前RB的各个DET求最小值，注意Det的列数大于码本size，后边的值很小，itpp提供的函数会选出后边没定义的最小值
		//CodeBookIndex(i) = itpp::min_index(Det.get_row(i));自己写一个
		CodeBookIndex(i) = 0;
		for(j=1;j<CodeBookSize;j++)
		{				
			if(Det(i,CodeBookIndex(i)) > Det(i,j))
				CodeBookIndex(i) = j;
		}
		UsedCodeBook[i] = CodeBook[CodeBookIndex(i)];//最佳码本
							 						 
	}
	
	//写入每个子带的预编码信息
	for(i=0;i<mSubBandNum;i++)
	{
		mPrecodingMatrixIndex[i] =  CodeBookIndex(i);
		mPrecodingMatrix[i] = UsedCodeBook[i];
	}

	/*输出
	for(i=0;i<mSubBandNum;i++)
	{
		cout << "RB Number             "<< i << endl;
		cout << "所选码本序号          "<< mPrecodingMatrixIndex << endl;//输出10个RB选择码本的序号
		cout << "precode matrix chosed " << endl << mPrecodingMatrix[i] << endl;//输出10个RB选择的预编码矩阵
		cout << "**********************"<<endl;
	}*/
	
	for(i=0;i<LayerSymbolNum;i++)
	{
		LayerMappedSymbol(1,i) =  mMoPhyChOutputData(2*i+1);
		LayerMappedSymbol(2,i) =  mMoPhyChOutputData(2*i+2);
	}
	for(i=0;i<LayerSymbolNum;i++)
	{
//		PrecodedSymbol.set_col(i,LayerMappedSymbol.get_col(i)*mPrecodingMatrix[i/mmSubBandWidth]);
	}
//	mPrecodedSymbol = PrecodedSymbol;

}


void CPDSCh::REmapper()
{
    //double scale_factor=0.34;
	int temp;
	cmat out;
	ivec Nsubframe;
	cvec input;
    
	int i,j,k,l,Nsc,len,Ndl;
	int num=mRBNum;
	
	if(mCarrierType==0) Nsc=12;
	else Nsc=24;
    //cout<<"mCarrierType="<<mCarrierType<<endl;
	
	out.set_size(num*Nsc,140,false);

    input.set_size(mTempDataSpace.cols(),false);
    for(i=0;i<mTempDataSpace.cols();i++)
		input(i)=mTempDataSpace(0,i);
	/*input.set_size(50,false);
	randn_c(50,input);*/
    //cout<<"REinput= "<<input<<endl;

	RElen=input.length();

	temp=(int)ceil((double)input.length()/(num*Nsc));
	//Ndl=(int)ceil((double)temp/14);
	//cout<<"temp"<<temp<<endl;
	
	out.zeros();
    len=mod(input.length(),num*Nsc);
	if(mUpDownlink==0)
	{
//		Ndl=1;
		//Nsubframe.set_length(10);
		//Nsubframe="0 1 2 3 4 5 6 7 8 9";
		
		//for(i=0;i<Ndl;i++)
		//{
			for(l=0;l<temp;l++)
			{   if(l<temp-1)
			    {
				 for(k=0;k<num*Nsc;k++)
				  {   
					/*cout<<"Nsubframe"<<Nsubframe<<endl;*///ahdsfjasfakdkjasdkjfdksdfkjsfdkjkjdfadssadsdfsfas
				   //if(l<14)
					//out(k,(mod(l,14)+14*Nsubframe(0)))=input(k+l*num*Nsc);
				   //else if
					//out(k,(mod(l,14)+14*Nsubframe(1)))=input(k+l*num*Nsc);
					  out(k,l)=input(k+l*num*Nsc);
					  //cout<<"xia biao "<<k+l*num*Nsc<<endl;
				  }
			    }
			   else
			  {
				  if(len==0)
				  {
					  for(k=0;k<num*Nsc;k++)
					  {
						  out(k,l)=input(k+l*num*Nsc);
					  }
				  }
				  else
				  {
					  for(k=0;k<len;k++)
				     {  
						 out(k,l)=input(k+l*num*Nsc);
				   /*if(l<14)
					out(k,(mod(l,14)+14*Nsubframe(0)))=input(k+l*num*Nsc);
				   else
					out(k,(mod(l,14)+14*Nsubframe(1)))=input(k+l*num*Nsc);*/
				    	
				     }
				  }
			 }
	   	  }
       
		//}
	}
    //cout<<"out= "<<out<<endl;
	//cout<<"RElen= "<<RElen<<endl;
	//cout<<"out="<<out.rows()<<endl;
    if(mUpDownlink==1)
	{
		Ndl=4;
		Nsubframe.set_length(Ndl);
		Nsubframe="0 4 5 9";
		for(i=0;i<Ndl;i++)
		{
			for(l=0;l<14*Ndl;l++)
			{
				for(k=0;l<num*Nsc;k++)
				{
					out(k,mod(l,14)+14*Nsubframe(i))=input(k+l*num*Nsc);
				}
			}
		}	
	}
    if(mUpDownlink==2)
	{
		Ndl=6;
		Nsubframe.set_length(Ndl);
		Nsubframe="0 3 4 5 8 9";
		for(i=0;i<Ndl;i++)
		{
			for(l=0;l<14*Ndl;l++)
			{
				for(k=0;l<num*Nsc;k++)
				{
					out(k,mod(l,14)+14*Nsubframe(i))=input(k+l*num*Nsc);
				}
			}
		}	
	}
    if(mUpDownlink==3)
	{
		Ndl=6;
		Nsubframe.set_length(Ndl);
		Nsubframe="0 5 6 7 8 9";
		for(i=0;i<Ndl;i++)
		{
			for(l=0;l<14*Ndl;l++)
			{
				for(k=0;l<num*Nsc;k++)
				{
					out(k,mod(l,14)+14*Nsubframe(i))=input(k+l*num*Nsc);
				}
			}
		}
	}
	if(mUpDownlink==4)
	{
		Ndl=7;
		Nsubframe.set_length(Ndl);
		Nsubframe="0 4 5 6 7 8 9";
		for(i=0;i<Ndl;i++)
		{
			for(l=0;l<14*Ndl;l++)
			{
				for(k=0;l<num*Nsc;k++)
				{
					out(k,mod(l,14)+14*Nsubframe(i))=input(k+l*num*Nsc);
				}
			}
		}
	}
	if(mUpDownlink==5)
	{
		Ndl=8;
		Nsubframe.set_length(Ndl);
		Nsubframe="0 3 4 5 6 7 8 9";
		for(i=0;i<Ndl;i++)
		{
			for(l=0;l<14*Ndl;l++)
			{
				for(k=0;l<num*Nsc;k++)
				{
					out(k,mod(l,14)+14*Nsubframe(i))=input(k+l*num*Nsc);
				}
			}
		}	
	}
	if(mUpDownlink==6)
	{
		Ndl=3;
		Nsubframe.set_length(Ndl);
		Nsubframe="0 5 9";
		for(i=0;i<Ndl;i++)
		{
			for(l=0;l<14*Ndl;l++)
			{
				for(k=0;l<num*Nsc;k++)
				{
					out(k,mod(l,14)+14*Nsubframe(i))=input(k+l*num*Nsc);
				}
			}
		}
		
	}
	/*for(l=0;l<3;l++)
	{
	    for(k=0;k<num*Nsc;k++)
	    {
		   out(k,l)=Z(k+l*num*Nsc);
	    }
	}

    for(l=4;l<symb;l++)
	{
	    for(k=0;k<num*Nsc;k++)
	    {
		   out(k,l)=Z(k+(l-1)*num*Nsc);
	    }
	}*/

	mMoPhyChOutputData.set_size(out.rows(),out.cols(),false);
	for(i=0;i<out.rows();i++)
       for(j=0;j<out.cols();j++)
	     mMoPhyChOutputData(i,j)=out(i,j);
	//cout<<"REmapping"<<mMoPhyChOutputData<<endl;
	
//	return out;
}

void CPDSCh::DeREmapper()
{
	cmat input;
	cvec out;
   
	int i,j,k,l,temp,Nsc,len;
	int num=mRBNum;
	
	if(mCarrierType==0) Nsc=12;
	else Nsc=24; 

    //k=0;
	input.set_size(mDePhyChInputData.rows(),mDePhyChInputData.cols(),false);
    for(i=0;i<input.rows();i++)
	{
       for(j=0;j<input.cols();j++)
	   {
		  // cout<<"mDePhyChInputData="<<mDePhyChInputData(i,j)<<endl;
		   input(i,j)=mDePhyChInputData(i,j);
	       /*if(input(i,j).real()!=0||input(i,j).imag()!=0)
			   k++;*/
	   }
	}
	//cout<< "input RE 解映射"<< input<<endl;

	out.set_size(RElen,false);
	out.zeros();

	temp=(int)ceil((double)out.length()/(num*Nsc));
	len=mod(out.length(),num*Nsc);
    
	for(l=0,i=0;l<temp;l++)
	{   if(l<temp-1)
	    {
			for(k=0;k<num*Nsc;k++)
			  {
				out(i)=input(k,l);
				i++;
			  }
		}
	    
		else
		{   if(len==0)
	        {
				for(k=0;k<num*Nsc;k++)
			   {  //cout<<"i=  "<<i<<endl;
				  out(i)=input(k,l);	
				  i++;
			   }
	        }
		    else
			{
				for(k=0;k<len;k++)
			   { // cout<<"i=  "<<i<<endl;
				  out(i)=input(k,l);	
				  i++;
			   }
	        }
		    
		}
    }
	//for(j=0,k=0;j<input.cols();j++)
	//  for(i=0;i<input.rows();i++)		  
	//   //if(input(i,j).real()!=0||input(i,j).imag()!=0)
	//	 {
	//		out(k)=input(i,j);
	//		k++;
	//	 }


    //cout<<"out="<<out<<endl;
	//cout<<"length = "<<out.length()<<endl;
    //cout<<"i=  "<<i<<endl;
	//cout<<"******************"<<endl;

	mTempDataSpace.set_size(1,out.length(),false);
       for(j=0;j<mTempDataSpace.cols();j++)
	     mTempDataSpace(0,j)=out(j);

//	return out;
}


void CPDSCh::Run()
{
    mDLSCh.Run();
	mMoPhyChInputData=mDLSCh.mpTrChOutputData; //从传输信道得到数据
	
    //Scrambling();
    Modulation();

    REmapper();


}
void CPDSCh::DRun(double noisevar)
{  
    int i;
    //mPhyChInputData=GenScramblingSeq();  //从接收机接受数据
	//double error_sum;
	DeREmapper();
	//cout<<"调制后"<<mTempDataSpace2.get_row(0)<<endl<<"******************"<<endl;
	//cout<<"%%%%%%%%%%%%%"<<mTempDataSpace2.cols()<<endl;
	//cout<<"解RE 映射后 "<<mTempDataSpace.get_row(0)<<endl;
	//cout<<"&&&&&&&&&&&&&"<<mTempDataSpace.cols()<<endl;
    Demodulation(noisevar);
    //cout<<"调制前比特"<<endl<< mMoPhyChInputData<<endl;
	//cout<<"解调后比特"<<endl<<mTempDataSpace1.get_row(0)<<endl;
	//cout<<"zhou compare"<<endl<<mMoPhyChInputData-to_bvec(mTempDataSpace1.get_row(0))<<endl;
	//error_sum=0;
	//for(i=0;i<mMoPhyChInputData.length();i++)
	//{error_sum=error_sum+abs((int)mMoPhyChInputData(i)-(int)mTempDataSpace1(0,i));
	//}
	//cout<<"total error "<<error_sum<<endl;
	//cout<<"error ratio "<<error_sum/mTempDataSpace1.cols()<<endl;
	//cout<<mMoPhyChInputData.length()<<endl<<mTempDataSpace1.cols();
	mDLSCh.mpTrChInputData.set_size(mTempDataSpace1.cols(),false);
	 for(i=0;i<mDLSCh.mpTrChInputData.length();i++)
       mDLSCh.mpTrChInputData(i)=mTempDataSpace1(0,i);

	mDLSCh.DRun();

}


   
void  CPDSCh::LayerMapping()   // Layer mapping
{
//	return 0;
}
    
void  CPDSCh::DePrecoding()    // Inverse of precoding for PUSCH: FFT
{
//	return 0;
}
   
void  CPDSCh::DeLayerMapping()   // Inverse of Layer mapping
{
//	return 0;
}




//88*****************************************************************88
//上行控制信道类成员函数
CPUCCh::CPUCCh() //构造函数
{ //mFormatofPUCCh=1; // PUCCH format: 0—1, 1a, 1b; 1—2, 2a, 2b
 }

CPUCCh::~CPUCCh(){}


void CPUCCh::InniCPhyChannel(LCP para,GP para1)
 {
    mUpDownlink=para.mUpDownlink;
	mUpCPType=para.mUpCPType;             //Uplink 0:Normal;1:Extended 
	mFormatofPUCCh=para.mFormatofPUCCh;
	mRBNum=para1.mNumofDwRB;
	//cout<< "mFormatofPUCCh="<<mFormatofPUCCh<<endl;
}

void CPUCCh::Run()
{
	/*int L=24,Z=50; 
	 bvec bits = randb(60),v,debits;
      Array<bvec> SegCrcbits;
     CTranChannel a;
     v=a.AttachCrc(bits,L);
     cout<<"p="<< v<<endl;
     SegCrcbits=a.SegCodBlockCRC(v,Z); 
  cout<<"segbits="<<SegCrcbits<<endl;
  debits=a.DeSegAttachCrc(SegCrcbits,v);
  cout<<"debits="<<debits.size()<<endl;
  debits=a.DeAttachCrc( L,debits);
  cout<<"debits2="<<debits<<endl;
  v=bits-debits;
  cout<<"b2="<<v<<endl;*/


UCI.Run();
UCI.DRun();
int mFormatofPUCCh=4,j;//mFormatofPUCCh的取值要在用户层修改
bvec c=randb(20),b;//
//int numsubframe;
vec d01;
cvec d02,space,d03;
cpl d10,d0;
 

if(mFormatofPUCCh==1||mFormatofPUCCh==2||mFormatofPUCCh==3)
{
 if(mFormatofPUCCh==1||mFormatofPUCCh==2)
 {cvec dd(c.size());
 BPSK bpsk;
 d01=bpsk.modulate_bits(c);
for(j=0;j<d01.size();j++)
{dd(j)=cpl(d01(j),0);}
d02=dd;
//cout<<"p="<< d01<<endl;
}
else
  {QPSK qpsk;
 d02=qpsk.modulate_bits(c);
//cout<<"c="<< c<<endl;
 }
//
for(j=0;j<d02.size();j++)
{ d0=d02(j);
//cout<<"111="<<d0<<endl;
	
   ComplexSymShift1(d0);//
   //cout<<"space="<<mCycShiftSeq<<endl;
   mCycShiftSeq=BlockWiseSpread(mCycShiftSeq);
   //cout<<"mCycShiftSeq="<<mCycShiftSeq.size()<<endl;
 }
}//PUCCH format 1,1a
else if(mFormatofPUCCh==4||mFormatofPUCCh==5||mFormatofPUCCh==6)
{    
	 b=ScrambLing(c(0,19));//
     //cout<<"b="<< b<<endl;

	if(mFormatofPUCCh==4)   //PUCCH format 2/2b , QPSK modulation
	{QPSK qpsk;
    d02=qpsk.modulate_bits(b);
    //cout<<"p"<< d02<<endl;
	}
	else if(mFormatofPUCCh==5)
	{
    BPSK bpsk;
	QPSK qpsk;
    d02=qpsk.modulate_bits(b);
	d01=bpsk.modulate_bits(c(20,20));
    d10=cpl(d01(0),0);
	//cout<<"space="<<d10<<endl;
	}
	else if(mFormatofPUCCh==6)
	{QPSK qpsk;
    d02=qpsk.modulate_bits(b);
	d03=qpsk.modulate_bits(c(20,21));
    //cout<<"p"<< d02<<endl;
	}
  
    mCycShiftSeq=ComplexSymShift2(d02);}//
	////cout<<"zzz="<<space<<endl;
}
void CPUCCh:: DRun()
{
//DeREmapper();
int d01,d0;//mFormatofPUCCh的取值要在用户层修改
bvec b,d02;
//int numsubframe;
cvec space,d03;
cpl d10;

if(mFormatofPUCCh==1||mFormatofPUCCh==2||mFormatofPUCCh==3)
{ 
	DeBlockWiseSpread(mCycShiftSeq);
	//cout<<"mCycShiftSeq="<<mCycShiftSeq<<endl;
    d03=DeComplexSymShift1(mCycShiftSeq);
	//cout<<"mCycShiftSeq="<<d03<<endl;
    
	if(mFormatofPUCCh==1||mFormatofPUCCh==2)
{ 
	d01=d03(0).real();
	 if(d01==1)
	 {d0=0;}
	 else if(d01==-1)
	 {d0=1;}
	 //cout<<"bits="<<d0<<endl;
	}
else
  {QPSK qpsk;
  d02=qpsk.demodulate_bits(d03);
 //cout<<"bits="<< d02<<endl;
	 }
}

else if(mFormatofPUCCh==4||mFormatofPUCCh==5||mFormatofPUCCh==6)
{ 
    DeComplexSymShift2( mCycShiftSeq);

	 //b=ScrambLing(c(0,19));
  //   cout<<"b="<< b<<endl;

	if(mFormatofPUCCh==4)   //PUCCH format 2/2b , QPSK modulation
	{QPSK qpsk;
    d02=qpsk.demodulate_bits(mCycShiftSeq);
    //cout<<"p"<< d02<<endl;
	}
	else if(mFormatofPUCCh==5)
	{
    BPSK bpsk;
	QPSK qpsk;
    d02=qpsk.demodulate_bits(mCycShiftSeq);
	//d01=bpsk.demodelate_bits(c(20,20));
   /* d10=cpl(d01(0),0);
	cout<<"space="<<d10<<endl;*/
	}
	else if(mFormatofPUCCh==6)
	{QPSK qpsk;
    d02=qpsk.demodulate_bits(mCycShiftSeq);
	/*d03=qpsk.modulate_bits(c(20,21));
    cout<<"p"<< d02<<endl;*/
	}}
 cout<<"decode_bits"<< d02<<endl;
 //    ComplexSymShift2(d02);
	////cout<<"zzz="<<space<<endl;
     
}

void  CPUCCh::ComplexSymShift1(cpl d0)//Complex-valued modulation symbols cyclically shift

{
//const double PI=3.14159265;
// //cpl d0=cpl(0,1);
// double a=1;   //求循环移位序列，需要确定a）
// int n,Nseq=12;// cyclically shifted length
// cvec y(Nseq);
// cmat ruv(a,Nseq);
//
// vec phase="-1 1 3 -3 3 3 1 1 3 1 -3 3"; 
// for (int n=0;n<Nseq;n++)
// {ruv(0,n)=cpl(cos(phase(n)*PI/4+a*n),sin(phase(n)*PI/4+a*n));
// //cout<<"1="<<ruv(0,n)<<endl;
//  y(n)=d0*ruv(0,n);}
// mCycShiftSeq=y;
// //cout<<"2="<< y<<endl; 
};

cvec  CPUCCh::DeComplexSymShift1(cvec y)//Complex-valued modulation symbols cyclically shift
//void main()
{//const double PI=3.14159265;
 cvec d0(1);
 int a=1;   //求循环移位序列，需要确定a）
 int Nseq=12;// cyclically shifted length
 cmat ruv(a,Nseq);

 vec phase="-1 1 3 -3 3 3 1 1 3 1 -3 3"; 
 ruv(0,0)=cpl(cos(phase(0)*PI/4+a*0),sin(phase(0)*PI/4+a*0));
  d0(0)=y(0)/ruv(0,0);
return d0;
};
cvec CPUCCh::ComplexSymShift2(cvec d)
//void main()
{//const double PI=3.14159265;
 //cpl d0=cpl(0,1);
 int a=1;   //求循环移位序列，需要确定a）
 int n,Nseq=12,NscRB=12;// cyclically shifted length,NSCRB 常规和扩展都是12
 cvec y(Nseq);
 cmat ruv(a,Nseq);
 cvec z(10*NscRB);

 vec phase="-1 1 3 -3 3 3 1 1 3 1 -3 3";
 for(n=0;n<10;n++)
 for (int j=0;j<NscRB;j++)
 {ruv(0,j)=cpl(cos(phase(j)*PI/4+a*j),sin(phase(j)*PI/4+a*j));
 //cout<<"1="<<ruv(0,n)<<endl;
  z(Nseq*n+j)=d(n)*ruv(0,j);}
 //cout<<"2="<< y<<endl; 
 return z;
//}//cyclic shift variable a
};

void CPUCCh::DeComplexSymShift2(cvec z)
//void main()
{//const double PI=3.14159265;
 //cpl d0=cpl(0,1);
 int a=1;   //求循环移位序列，需要确定a） a double???
 int n=0,Nseq=12,NscRB=12;// cyclically shifted length,NSCRB 常规和扩展都是12
 cvec y(Nseq);
 cmat ruv(a,Nseq);
 cvec d(11);

  vec phase="-1 1 3 -3 3 3 1 1 3 1 -3 3";
 cout<<"z.size()="<<z.size()<<endl;
 ruv(0,0)=cpl(cos(phase(0)*PI/4+a*0),sin(phase(0)*PI/4+a*0));
 for (int n=0;n<10;n++)
 {d(n)=z(12*n)/ruv(0,0);
 //cout<<"d(n)="<<d(n)<<endl;
 } 
  mCycShiftSeq=d;
  //cout<<" mCycShiftSeq="<< mCycShiftSeq<<endl;
//}//cyclic shift variable a
};
cvec CPUCCh::BlockWiseSpread(cvec y)
//void main()
{
int	Nsf=4,Nseq=12;
int mp,m,n;
imat w="+1 +1 +1 +1;+1 -1 +1 -1;+1 -1 -1 +1";
//cout<<"imatw="<<w.get_row(0)<<endl;
ivec ns="1:5";//时隙数目
int mOrthSeqIndexnoc=0;//需要从参数获得
cvec yy(Nseq);
cvec z(96);
//cvec y=" 2 3 4 1 2 1 3 4 2 3 1 2";

//确定是常规还是扩展
//int mPRBConfig ;  //在基类已经设置 Resource block Configuration(0- Normal cp/1- Extended cp)

//cout<<"imatw="<<ComplexSymShift( )<<endl;

//

for(mp=0;mp<2;mp++)
  {for(m=0;m<4;m++)
    {for(n=0;n<12;n++)
       {z(mp*Nsf*Nseq+m*Nseq+n)=w(mOrthSeqIndexnoc,m)*y(n);
         //cout<<"imatw="<<z(mp*Nsf*Nseq+m*Nseq+n)<<endl;
      }
    }
   }
//cout<<"imatw="<<z<<endl;
return z;
};

bvec  CPUCCh::ScrambLing(bvec b)
{bvec bp(20);
 bvec c=randb(20); //需要得知加扰序列c(i)

//cout<<"imatw="<<b<<"\n"<<"c="<<c<<endl;
 for(int j=0;j<20;j++)
 {
bp(j)=mod((b(j)+c(j)),2);
//cout<<"bp="<<bp(j)<<endl;
 }
return bp;
}
void CPUCCh::DeBlockWiseSpread(cvec z)   // Inverse of block of complex-valued symbols spread
{int i;
cvec y(12);
for(i=0;i<12;i++)
{y(i)=z(i);}
 mCycShiftSeq=y;
}
cvec DeREmapper()         // Inverse of mapping
{return 0;}
//88*****************************************************************88
//下行控制信道类成员函数
CPDCCh::CPDCCh()
{mIMbit=0;}

CPDCCh::~CPDCCh(){}


void CPDCCh:: Run()
{ cout<<"imatw="<<1<<endl;
CDCI a;
a.Run();
a.DRun();
/*BlockMultiplexed();

Scrambling();

Modulation ();

LayerMapping();

Precoding();

REmapper();    */     //RE mapping
}


void CPDCCh:: DRun()
{ 
//DeREmapper();
//
//DePrecoding();
//
//DeLayerMapping();
//
//DeScrambling();
//
//DeBlockMultiplexed();
}


cvec CPDCCh::BlockMultiplexed()
{
return 0;}
cvec CPDCCh::LayerMapping()
{
return 0;}
cvec CPDCCh::Precoding()
{
return 0;}
cvec CPDCCh:: DePrecoding()
{
return 0;}
cvec CPDCCh:: DeLayerMapping()
{
return 0;}
cmat DeModulation ()
{
return 0;}
 cmat DeScrambling ()
{
return 0;}
cvec CPDCCh:: DeBlockMultiplexed()
{ return 0;
}

