//----- PParaInit.h -----
// Programmed by LTE Simulation Group of IMC, SWJTU
#ifndef ParaInit_HPP
#define ParaInit_HPP

#include "stdafx.h"
#include "itpp/itcomm.h"
#include <string>
#include <fstream>

using std::cin;
using std::cout;
using std::endl;
using namespace itpp;
using namespace std;

typedef struct GlobalParameter
{
	 int mNumofTxAntenna;       // The number of Transmitter antenna,1/2/4.
     int mNumofRxAntenna;       // The number of Receiver antenna,1/2/4.
     int mPUSChState;           // The state of PUSCH:on=1; off=0
     int mPUCChState;           // The state of PUCCH:on=1; off=0
     int mPDSChState;           // The state of PDSCH:on=1; off=0
     int mPDCChState;           // The state of PDCCH:on=1; off=0
     int mRSTypeState[4];       // The state of Reference signals:0:no reference signals transmitted,1:Demodulation reference signal,2:Sounding reference signal
                                // 3:common reference signal,4:dedicated reference signal
     int mNumofUpRB;            // Uplink transmission bandwidth configured,[6,110]
     int mNumofDwRB;            // Downlink transmission bandwidth configured,[6,110]
     int mSubcarrierType;       // The choice of Subcarrier spacing (0: 15kHz 1: 7.5kHz) 
     int mNumofSimulationLoop;  // The number of simulation loop
}GP;

typedef struct LogicChannelParameter
{
     int mUpDownlink;           // from 0 to 6;
     int mUpCPType;             // Uplink 0:Normal;1:Extended 
	 int mDownCPType;           // Downlink 0:Normal;1:Extended 
	 int mUpModType;            // Uplink Modulation 0: QPSK 1: 16QAM 2: 64QAM
	 int mDownModType;          // Downlink Modulation 0: QPSK 1: 16QAM 2: 64QAM
	 int mFormatofPUCCh;        // PUCCH format
	 int mFormatofPDCCh;        // PDCCH format
	 int mTranBlock;            // The size of transmission block
	 int mSoftBuffSize;         // The size of Soft Buffer        
	 int mReduVersiNum;         // The redundancy version number
	 int mStreamNum;            // Number of stream
	 int mLayerNum;             // Number of Layer

}LCP;

typedef struct  ReferenceSignalParameter
{
     int mGroupnumber;          // The group number u of reference sequence (Defined by users)
	 int mLengthsequence;       // Length of reference sequence to determine m,v=0 when 1<=m<=5,v=0»ò1£¬when 6<=m<=110(Define the RS scheme)
     //int mRSState[4];         // Reference signal state(four)

}RSP;

typedef struct VectorChannelParameter
{  
     int mModelCase;            // Model
	 cvec mVelocityOfMS;        // Velocity Of MS
	 double mFrequencyOfCarrier;// Frequency of carrier wave
	 int mNumOfPath;            // Number of path
	 vec mPathDelays;           // Delay each path
	 double mDoppler;           // Doppler
	 double mAntennaSpacingBS;  // Antenna Spacing in BS
	 double mAntennaSpacingMS;  // Antenna Spacing in MS
	 double AODPerPathBS;       // AOD Per Path in BS
	 double AOAPerPathMS;       // AOA Per Path in MS
	 vec mAngleSpreadPerPathMS; // Angle Spread Per Path in MS
	 vec mAngleSpreadPerPathBS; // Angle Spread Per Path in BS
	 int NumOfSample;           // Number of sample

}VCP;

//If these parameters are needed, functions below will be called, then datas are obtained from the file

//--- Definition of initGP()
void initGP(string filename,GP &a)
{
    
    ifstream fin(filename.c_str());
    if(!fin) return;

	int a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14;
 
    fin>>a1>>a2>>a3>>a4>>a5>>a6>>a7>>a8>>a9>>a10>>a11>>a12>>a13>>a14;
    {
		a.mNumofTxAntenna=a1;
		a.mNumofRxAntenna=a2;
		a.mPUSChState=a3;
		a.mPUCChState=a4;
		a.mPDSChState=a5;
		a.mPDCChState=a6;
		a.mRSTypeState[0]=a7;
		a.mRSTypeState[1]=a8;
		a.mRSTypeState[2]=a9;
		a.mRSTypeState[3]=a10;
		a.mNumofUpRB=a11;
		a.mNumofDwRB=a12;
		a.mSubcarrierType=a13;
		a.mNumofSimulationLoop=a14;
    }
    fin.close();

}

//--- Definition of initLCP()
void initLCP(string filename,LCP &b)
{
    ifstream fin(filename.c_str());
    if(!fin) return;

	int a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12;
 
    fin>>a1>>a2>>a3>>a4>>a5>>a6>>a7>>a8>>a9>>a10>>a11>>a12;
    {
		b.mUpDownlink=a1;
		b.mUpCPType=a2;
		b.mDownCPType=a3;
		b.mUpModType=a4;
		b.mDownModType=a5;
		b.mFormatofPUCCh=a6;
		b.mFormatofPDCCh=a7;
		b.mTranBlock=a8;
		b.mSoftBuffSize=a9;
		b.mReduVersiNum=a10;
		b.mStreamNum=a11;
		b.mLayerNum=a12;
    }
    fin.close();
}

//--- Definition of initRSP()
void initRSP(string filename,RSP &c)
{
    ifstream fin(filename.c_str());
    if(!fin) return;

	int a1,a2;
 
    fin>>a1>>a2;
    {
		c.mGroupnumber=a1;
		c.mLengthsequence=a2;
	/*	c.mRSState[0]=a3;
		c.mRSState[1]=a4;
		c.mRSState[2]=a5;
		c.mRSState[3]=a6;*/
    }
    fin.close();
}

//--- Definition of initVCP()
void initVCP(string filename,VCP &d)
{
	ifstream fin(filename.c_str());
    if(!fin) return;

	int a1,a4,a13;
	double a3,a6,a7,a8,a9,a10;
	cvec a2;
	vec a5,a11,a12;
 
    fin>>a1>>a2>>a3>>a4>>a5>>a6>>a7>>a8>>a9>>a10>>a11>>a12>>a13;
    {
		d.mModelCase=a1;
		d.mVelocityOfMS=a2;
		d.mFrequencyOfCarrier=a3;
		d.mNumOfPath=a4;
		d.mPathDelays=a5;
		d.mDoppler=a6;
		d.mAntennaSpacingBS=a7;
		d.mAntennaSpacingMS=a8;
		d.AODPerPathBS=a9;
		d.AOAPerPathMS=a10;
		d.mAngleSpreadPerPathMS=a11;
		d.mAngleSpreadPerPathBS=a12;
		d.NumOfSample=a13;
    }
    fin.close();
}

#endif