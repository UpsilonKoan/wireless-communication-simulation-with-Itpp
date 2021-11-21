//----- Tranceiver.h -----
// Programmed by LTE Simulation Group of IMC, SWJTU
#include "stdafx.h"
#include "itpp/itcomm.h"
#include "ReferenceSignal.h"
#include "PhyChannel.h"
#include "VectorChannel.h"
#include "ParaInit.h"

using std::cin;
using std::cout;
using std::endl;
using namespace itpp;
using namespace std;

//************************************* Class CTranceiver Declaration *************************************************************
class CTranceiver
{
		/***** Data Members *****/
protected:
		 bvec mInputData;               // Input Data
		 bvec mOutputData;              // Output Data
		 cmat TempSpace;                // Buffer
		 int mPUSChState;               // The state of PUSCH:on=1; off=0
		 int mPUCChState;               // The state of PUCCH:on=1; off=0
		 int mPDSChState;               // The state of PDSCH:on=1; off=0
		 int mPDCChState;               // The state of PDCCH:on=1; off=0
		 int mRSTypeState[4];           // The state of Reference signals
public:
	//int mRedVersionNum;               // Number of Redundancy Version ***
	int mNumberOfError;                 // Number of error   **************
	int mNumberOfRetr;                  // Number of retransmission ******* 
	int mNumofDwRB;                     // Number of RB
	double EbN0dB;                      // Eb/N0    
public:
	CReferenceSignal ReferenceSignal;   // Object of reference signal class
	CVectorChannel  VectorChannel;      // Object of vector channel class
	CPDSCh PDSCh;                       // Object of pdsch class
	CPUSCh PUSCh;                       // Object of pusch class
	//----- Hold In Abeyance -----
	//CPDCCh PDCCh;                     // Object of pdcch class
	//CPUCCh PUCCh;                     // Object of pucch class

		/***** Function Members *****/
	    /***** Constructors *****/
public:
	CTranceiver();
	/*------------------------------------------------------------------------
      Construct a Tranceiver object.
      Precondition:  None.
      Postcondition: An empty Tranceiver object has been constructed
	------------------------------------------------------------------------*/
	virtual ~ CTranceiver();
	/*-----------------------------------------------------------------------
      Class destructor 
	  Precondition:  None
	  Postcondition: The Tranceiver object has been destructed.
    ------------------------------------------------------------------------*/
	int Run();
	/*-----------------------------------------------------------------------
      Run all functions of CTranceiver 
	  Precondition:  Every function member works properly
	  Postcondition: Return the result of Judgement which is the number of RV
    ------------------------------------------------------------------------*/
	void initRun(GP &a,LCP &b,RSP &c,VCP &d);
	/*-----------------------------------------------------------------------
      Initialize Run function with essential parameters, i.e GP, LCP, RSP and
	  VCP, which are read from in1.txt and so on
	  Precondition:  Data file is available and ParaInit.h works fine
	  Postcondition: Essential parameters for CTranceiver are obtained
    ------------------------------------------------------------------------*/
protected:
	void initTranceiver(GP a);
	/*-----------------------------------------------------------------------
      Initialize Tranceiver with parameters GP to update mRSTypeState ...
	  Precondition:  Data file is available and ParaInit.h works fine
	  Postcondition: Essential parameters such as mRSTypeState,mPDSChState
	  and so on are obtained
    ------------------------------------------------------------------------*/
	void GenOFDMSig(cmat ResourceGrid);   
	/*-----------------------------------------------------------------------
      Generate OFDM Signals
	  Precondition:  ResourceGrid is inputed
	  Postcondition: OFDM Signals are generated and stored to TempSpace 
    ------------------------------------------------------------------------*/
	void GenSCFDMASig(); //generate SC-FDMA signal
	/*-----------------------------------------------------------------------
      //----- Hold In Abeyance -----
    ------------------------------------------------------------------------*/
	void  AWGNChannel(double &N0) ;        // additive white Gaussian noise
	/*-----------------------------------------------------------------------
      Transmit the signals through AWGN channel
	  Precondition:  Data in TempSpace which actually are OFDM Signals
	  Postcondition: Data in TempSpace are changed which is signals received now
    ------------------------------------------------------------------------*/
	void DemOFDMSig(cmat ofdmsym);
	/*-----------------------------------------------------------------------
      Demodulate OFDM Signals
	  Precondition:  ResourceGrid is inputed
	  Postcondition: OFDM Signals are generated and stored to TempSpace 
    ------------------------------------------------------------------------*/
	int  Judgement();           // Judgment
	 /*---------------------------------------------------------------------------------
	   Judge CRC result and whether Input is equal to Output,then whether a new block 
	   should be transmitted
	   Precondition:  Current RV for this transmission, CRC result, and the buffer of 
	                  Input and Output is available to access
	   Postcondition: RV, Number of Block Error and Number of Retransmission is updated
	 ---------------------------------------------------------------------------------*/
};