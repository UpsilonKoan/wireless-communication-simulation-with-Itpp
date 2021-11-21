//----- PhyChannel.h -----
// Programmed by LTE Simulation Group of IMC, SWJTU
#ifndef PhyChannel_HPP
#define PhyChannel_HPP

#include "stdafx.h"
#include <complex>
#include <iostream>
#include "itpp/itcomm.h"
#include <complex>
#include "ParaInit.h"
#include "TranChannel.h"

typedef std::complex<double> complex;
typedef std::complex<double> cpl;

using std::cin;
using std::cout;
using std::endl;
using namespace itpp;

//************************************* Class CPhyChannel Declaration *************************************************************

class CPhyChannel
{
	/***** Data Members *****/
protected:   
   cmat mTempDataSpace;                    // Temporary data space
   mat mTempDataSpace1;                    // Temporary data space
   cmat mTempDataSpace2;                   // test 
   vec mScrambSeq;                         // Scrambling sequence
   int mUpDownlink;
   int mPRBConfig;                         // Resource block Configuration(0- Normal cp/1- Extended cp)
   int mModulType;                         // Modulation type
   int mRBNum;                             // Number of RB
   int mSoftbuffSize;
   int mCarrierType;                       // Subcarrier type (0:12,1:24)

public:
	int RElen;                             // 用于保存RE的长度
	bvec mMoPhyChInputData;                // Modulation Input data from CTranChannel
    cmat mDePhyChInputData;                // Demodulation input data from ofdm
    cmat mMoPhyChOutputData;               // Physical channel output data to ofdm
	bvec mDePhyChOutputData;               // Physical channel output data to CTranChannel
	cmat mPrecodingMatrix[10];             // 每个子带选择的预编码矩阵
	ivec mPrecodingMatrixIndex[10];        // 每个子带选择的预编码矩阵索引号
 
	
	/***** Function Members *****/
 public:
   CPhyChannel();
   virtual ~ CPhyChannel(); 
   void GenScramblingSeq();                // Generate scrambling sequence
   void Scrambling();                      // Scrambling
   void Modulation ();                     // Modulation Mapper
   void MIMODec();                         // MIMO detection  
   void Demodulation (double noisevar);    // Demodulation Mapper, include MIMO detection
   void DeScrambling();                    // Inverse of Scrambling

};


class CPDSCh:public CPhyChannel
{
		/***** Data Members *****/
  private:
   int mStreamNum;                         // Number of stream
   int mLayerNum;                          // Number of Layer
   int RxAntennaNum;                       // Number of antenna in Rx
   int TxAntennaNum;                       // Number of antenna in Tx  
       
   //************************zhou!!************************
   int mSubBandNum;                        // 子带个数
   int mmSubBandWidth;                     // 每个子带RB数
   int RBWidth;                            // 每个RB子载波数
    cmat mPrecodingMatrix[10];             // 每个子带选择的预编码矩阵
	ivec mPrecodingMatrixIndex[10];        // 每个子带选择的预编码矩阵索引号
	cmat *mPrecodedSymbol;                 // 预编码输出数据
	//**********************************

  public:
    CDLSCh mDLSCh;                         // Object of DL-SCH class

	/***** Function Members *****/
  public:
    CPDSCh();
    virtual~CPDSCh();
    void Run();                            // Run all Modulation functions of CPDSch 
    void DRun(double noisevar);            // Run all Demodulation functions of CPDSch
	void InniCPhyChannel(LCP para,GP para1);
    void Precoding();
    void LayerMapping();                   // Layer mapping
    void DePrecoding();                    // Inverse of precoding for PUSCH: FFT
    void DeLayerMapping();                 // Inverse of Layer mapping
	void REmapper();                       // RE mapper
    void DeREmapper();                     // Inverse of RE mapper
};


class CPUSCh:public CPhyChannel
{
	/***** Data Members *****/
   private:
     int mAmpScaFactor;                   // Amplitude scaling factor
	 cmat *mPrecodedSymbolUp;             // 预编码输出数据

   public:
	 CULSCh mULSCh;                       // Object of UL-SCH class

   public:
	 CPUSCh();
	 virtual~CPUSCh();
     void Run();         //Run all modulation functions of CPUSch
	 void DRun();         //Run all demodulation functions of CPUSch
     void InniCPhyChannel(LCP para,GP para1);

   public:
	 void MultiDataControl();         //Data and Control multiplexing
	 void ChannelInterl();            //Channel Interleaver
	 void Precoding();               //Precoding for PUSCH: DFT
	 void DeMultiDataControl();      //Demultiplexing 
	 void DeChannelInterl();         //Channel Deinterleaver
	 void DePrecoding();            //Deprecoding
	 void REmapper();             // RE mapper
     void DeREmapper();           // Inverse of RE mapper
};


class CPUCCh: public CPhyChannel           //上行控制信道类设计
{
 private:
      int mFormatofPUCCh;       // PUCCH format: 0—1, 1a, 1b; 1—2, 2a, 2b
	  int mUpCPType; //normal cp or extended
	  
      cvec mCycShiftSeq;         //Cyclically shifted sequence
      int mResourceIndex;       // Resource index
      int mOrthSeqIndexnoc;        // Orthogonal sequence index
      cmat mSymbolofMod;      // Complex-valued symbol
 protected:

	  CUCI UCI;         // Object of UCI class

 public:
      CPUCCh();
      virtual~ CPUCCh ();
	  void Run();            //Run all Modulation functions of CPUCch
	  void DRun();            //Run all Demodulation functions of CPUCch
 public:
	  cvec BmodofPUCCh();      // Block of bits QPSK/BPSK modulation
	  void ComplexSymShift1(cpl d0);  //Complex-valued modulation symbols cyclically shift
	  cvec ComplexSymShift2(cvec d); 
	  cvec BlockWiseSpread(cvec y);   //Block of complex-valued symbols spread
	  bvec ScrambLing(bvec b);
	  //cvec ReMapping();         //Modulation symbols are mapped
	  cvec DeBmodofPUCCh();      // Block of bits QPSK/BPSK demodulation
	  cvec DeComplexSymShift();  // Inverse of Complex-valued modulation symbols cyclically shift
	  void DeBlockWiseSpread(cvec);   // Inverse of block of complex-valued symbols spread
	  void DeComplexSymShift2(cvec z);
	  cvec DeComplexSymShift1(cvec y);
	  void InniCPhyChannel(LCP para,GP para1);
	  //cvec DeReMapping();         // Inverse of mapping
};
class CPDCCh: public CPhyChannel       //下行控制信道类设计
{
 private:
     int mIMbit;              // Bit number in a subframe transmitted on PDCCH 
     int mNPDCCh;             //Number of PDCCHs transmitted in the subframe 
     cmat mSymbolofMod;     //Complex-valued modulation symbols
     int mLayerNum;//number of layers
//protected:
//	  CDCI DCI;         // Object of DCI class

 public:
     CPDCCh();
     virtual ~ CPDCCh ();
     void Run();            //Run all Modulation functions of CPDCch
     void DRun();            //Run all Demodulation functions of CPDCch
 public:
     cvec BlockMultiplexed();   //Block of bits shall be multiplexed
     cvec LayerMapping();         //Modulation symbols are mapped
     cvec Precoding();         //Modulation symbols are precoded
     //cvec REmapper();         //RE mapping
     cvec ComplexSymShift();  //Complex-valued modulation symbols cyclically shift
     cvec DeBlockMultiplexed();   // Inverse of Block of bits multiplexed
     cvec DeLayerMapping();         // Inverse of mapping
     cvec DePrecoding();         // Inverse of precoding
     cvec DeComplexSymShift();  // Inverse of Complex-valued modulation symbols cyclically shift
	 void InniCPhyChannel(LCP para,GP para1);
	 //cvec DeREmapper();         //Inverse of RE mapping

};
#endif