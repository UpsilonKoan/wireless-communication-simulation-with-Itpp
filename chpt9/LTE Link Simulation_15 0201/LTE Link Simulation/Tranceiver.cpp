//----- Tranceiver.cpp -----
// Programmed by LTE Simulation Group of IMC, SWJTU
#include "stdafx.h"
#include "itpp/itcomm.h"
#include "Tranceiver.h"
#include "ParaInit.h"
 
using std::cin;
using std::cout;
using std::endl;
using namespace itpp;

// *************************************** Class CTranceiver Definition *****************************************

//--- Definition of CTranceiver constructor
CTranceiver::CTranceiver()
{}

//--- Definition of CTranceiver destructor
CTranceiver::~ CTranceiver()
{}

//--- Definition of Run()
int CTranceiver::Run()
{
	double N0;

	/*VectorChannel.RunCVectorChannel();
    if(mRSTypeState[0]||mRSTypeState[1]||mRSTypeState[2]||mRSTypeState[3]) ReferenceSignal.Run();
    if(mPUSChState) PUSCh.Run();
	if(mPUCChState) PUCCh.Run();*/
	if(mPDSChState) PDSCh.Run();
	//if(mPDCChState) PDCCh.Run();
	//if(mPUSChState||str.mPUCChState) GenOFDMSig();
	if(mPDSChState||mPDCChState) GenOFDMSig(PDSCh.mMoPhyChOutputData);
    //VectorChannel.Convolution ();
    AWGNChannel(N0);
	/*if(mRSTypeState[0]||mRSTypeState[1]||mRSTypeState[2]||mRSTypeState[3]) ReferenceSignal.DRun();
	if(mPUSChState) PUSCh.DRun();
	if(mPUCChState) PUCCh.DRun();*/
	if(mPDSChState||mPDCChState) DemOFDMSig(TempSpace);
	PDSCh.mDePhyChInputData=TempSpace;
	if(mPDSChState) PDSCh.DRun(N0);
	//if(mPDCChState) PDCCh.DRun();
    return Judgement();
}

//--- Definition of initRun()
void CTranceiver::initRun(GP &a,LCP &b,RSP &c,VCP &d)
{
    initGP("in1.txt",a);
    initLCP("in2.txt",b);
    initRSP("in3.txt",c);
    initVCP("in4.txt",d);
	initTranceiver(a);
    PDSCh.InniCPhyChannel(b,a);
}

//--- Definition of initTranceiver()
void CTranceiver::initTranceiver(GP a)
{
    mRSTypeState[0]=a.mRSTypeState[0];
	mRSTypeState[1]=a.mRSTypeState[1];
	mRSTypeState[2]=a.mRSTypeState[2];
	mRSTypeState[3]=a.mRSTypeState[3];
	mPDCChState=a.mPDCChState;
	mPDSChState=a.mPDSChState;
	mNumofDwRB=a.mNumofDwRB; 
	/*mPUCChState=a.mPUCChState;
	mPUSChState=a.mPUSChState;*/
}
//--- Definition of GenOFDMSig()
void CTranceiver::GenOFDMSig(cmat ResourceGrid)   //cmat ResourceGrid(mNumofDwRB*12,20*7)
{
	int i,j,k;
	int CP=2;
	int Nslot=20;
    cmat ofdmsym1(mNumofDwRB*12,Nslot*7),ofdmsym2(CP,Nslot*7),ofdmsym(mNumofDwRB*12+CP,Nslot*7);
	for(int ns=0;ns<Nslot;ns++)
		for(int l=0;l<7;l++)
			ofdmsym1.set_col(ns*7+l,ifft(ResourceGrid.get_col(ns*7+l)));
    for(j=0;j<Nslot*7;j++)
		for(i=mNumofDwRB*12-CP,k=0;i<mNumofDwRB*12;i++)
		{
			ofdmsym2(k,j)=ofdmsym1(i,j);
			k++;
		}
	//ofdmsym2=ofdmsym1(mNumofDwRB*12-144,mNumofDwRB*12-1,0,20*7-1);
    ofdmsym=concat_vertical(ofdmsym2,ofdmsym1); // Add CP;

	TempSpace.set_size(ofdmsym.rows(),ofdmsym.cols(),false);
	for(i=0;i<TempSpace.rows();i++)
       for(j=0;j<TempSpace.cols();j++)
	     TempSpace(i,j)=ofdmsym(i,j);
	
}
//--- Definition of GenSCFDMASig()
void CTranceiver::GenSCFDMASig()
{}//----- Hold In Abeyance -----

//--- Definition of AWGNChannel()
void CTranceiver::AWGNChannel(double &N0)
{
	int i,j;
	int CP=2;
	int Nslot=20;
	int ML=2;    //QPSKµ÷ÖÆ½×Êý
	double SumPower=0,SignalPower,EbN0,sigma,CodeRate;
	cmat AWGNnoise(mNumofDwRB*12+CP,Nslot*7),AWGNsymbol(mNumofDwRB*12+CP,Nslot*7);

	EbN0 = (double)pow(10.0, EbN0dB/10.0); //Calculate Eb/N0 in a linear scale instead of dB.
	CodeRate = PDSCh.mDLSCh.mCodeRate;
    randn_c(mNumofDwRB*12+CP,Nslot*7,AWGNnoise);

	for(i=0;i<TempSpace.rows();i++)
	{
		SumPower = SumPower + TempSpace(i,1).real() * TempSpace(i,1).real() +TempSpace(i,1).imag() * TempSpace(i,1).imag();
	}
	
	SignalPower = SumPower/TempSpace.rows();
    sigma = sqrt(SignalPower/(2*EbN0*CodeRate*ML));
	N0 = 2*sigma*sigma;
    AWGNsymbol = TempSpace + sigma * AWGNnoise;
	
/**************** The method to add AWGN by IT++ *********************************
  AWGN_Channel awgn_channel;     //The AWGN channel class
    RNG_randomize();                //Randomize the random number generators in it++:
    N0 = 2*sigma*sigma;
	awgn_channel.set_noise(N0);
	AWGNsymbol.clear();
	for(i=0;i<TempSpace.rows();i++)
	{
		received_symbols= awgn_channel(TempSpace.get_row(i));
		for(j=0;j<AWGNsymbol.cols();j++)
			AWGNsymbol(i,j) = received_symbols(j);
	}
*********************************************************************************/

	TempSpace.set_size(AWGNsymbol.rows(),AWGNsymbol.cols(),false);
	for(i=0;i<TempSpace.rows();i++)
       for(j=0;j<TempSpace.cols();j++)
	   {
		   TempSpace(i,j)=AWGNsymbol(i,j);
	   }
}
//--- Definition of DemOFDMSig()
void CTranceiver::DemOFDMSig(cmat ofdmsym)  //ofdmsym(mNumofDwRB*12+144,20*7), ofdmsym is TempSpace actually
{
	int CP=2,Nslot=20;
	int i,j,k;
	cmat demofdmsym(mNumofDwRB*12,Nslot*7),ResourceGrid(mNumofDwRB*12,Nslot*7);
    //demofdmsym=ofdmsym(0,143,0,20*7-1);
	for(j=0;j<Nslot*7;j++)
		for(i=CP,k=0;i<mNumofDwRB*12+CP;i++)
		{
			demofdmsym(k,j)=ofdmsym(i,j);
			k++;
		}

	for(int ns=0;ns<Nslot;ns++)
		for(int l=0;l<7;l++)
			ResourceGrid.set_col(ns*7+l,fft(demofdmsym.get_col(ns*7+l)));

	TempSpace.set_size(ResourceGrid.rows(),ResourceGrid.cols(),false);
	for(i=0;i<TempSpace.rows();i++)
       for(j=0;j<TempSpace.cols();j++)
	     TempSpace(i,j)=ResourceGrid(i,j);
	//cout<<"ResourceGrid= "<< ResourceGrid <<endl;
}
//--- Definition of Judgement()
int CTranceiver::Judgement()
{
	mNumberOfRetr = PDSCh.mDLSCh.mRedunVersionNum;    // Update the number of retransmission
	mInputData  = PDSCh.mDLSCh.mUserData;             // The bits transmitted
	mOutputData = PDSCh.mDLSCh.mpTrChOutputData;      // The bits received	
	if ( PDSCh.mDLSCh.mCorrectness == 1)              // CRC result(mCorrectness) is right
	{
		PDSCh.mDLSCh.mRedunVersionNum = 0;
		if ( mOutputData != mInputData )
		{
			mNumberOfError++;
		}
	}
	else      // mCorrectness=0
		if ( PDSCh.mDLSCh.mRedunVersionNum < 3 )      // RV is less than 3, RV++
			PDSCh.mDLSCh.mRedunVersionNum++;
		else  // When mRedVersionNum = 3 
		{
			PDSCh.mDLSCh.mRedunVersionNum = 0;        // Transmission fails with 3 times of retransmission, block error happens
			mNumberOfError++;
		}
	return PDSCh.mDLSCh.mRedunVersionNum;
}