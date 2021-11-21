//----- LTE.cpp -----
// Programmed by LTE Simulation Group of IMC, SWJTU

#include <string>
#include <fstream>
#include "stdafx.h"
#include "itpp/itcomm.h"
#include "Tranceiver.h"
#include "ParaInit.h"
#include "TranChannel.h"
#include "VectorChannel.h"
#include "PhyChannel.h"
#include "ParaInit.h"

using namespace itpp;
using namespace std;

int main()
{
	//----- Scalar Definition -----

	int mCounter;                              // Loop Counter
	int mNumberOfError;                        // Number Of Error
	const int mNumOfSimulationLoop = 100000;   // Number Of Simulation Loop For One dB Vaule
	int mNumberOfRetr[mNumOfSimulationLoop];   // Number Of Retransmission For Each Simulation Loop

	//----- Vector Definition -----
	vec EbN0dB,block_error_rate;               // Eb/N0 and BLER

	//----- Hold In Abeyance -----

	//int mNumofTxAntenna;        // Number Of Tx Antenna   
	//int mNumofRxAntenna;        // Number Of  Rx Antenna
	//int mPhyChannelState[4] ;   // Physical Channel State
	//int mRSState[4] ;           // Reference Signal State
	//int mNumofUpRB;             // Uplink Transmission Bandwidth Configured
	//int mNumofDwRB;             // Downlink Transmission Bandwidth Configured
	//int mSubcarrierType;        // Subcarrier Spacing
	//float blockerrorrate;       // Block error Rate
	//cmat mDataSpace;            // The data space

	//----- Class Definition -----
	
	CTranceiver Tranceiver;       // Object Of CTranceiver
	GP gp_str;                    // Global Parameters
	LCP lcp_str;                  // Logic Channel Parameters
	RSP rsp_str;                  // Reference Signal Parameters
	VCP vcp_str;                  // Vector Channel Parameters

	//----- Program Begins -----
	
	EbN0dB = linspace(-4,-0.5,8);    // Simulate for 8 Eb/N0 values from -4 to -0.5 dB
	block_error_rate.set_size(EbN0dB.length(),false); //Allocate storage space for the result vector.
	//The "false" argument means "Do not copy the old content of the vector to the new storage area."
	for (int i=0; i<EbN0dB.length(); i++) 
	{
		//Show how the simulation progresses:
		cout << "Now simulating Eb/N0 value number " << i+1 << " of " << EbN0dB.length() << endl;
		//----- Initialization -----
		Tranceiver.initRun(gp_str,lcp_str,rsp_str,vcp_str);
		Tranceiver.PDSCh.mDLSCh.mRedunVersionNum=0;                       // Initialize the Reduncy Version Number for the first transmission
		mCounter=0;                                                       // Reset the Counter
		mNumberOfError=0;                                                 // Reset the Number of Block Error in main()
		Tranceiver.mNumberOfError=0;                                      // Reset the Number of Block Error in Class Tranceiver
		Tranceiver.EbN0dB=EbN0dB(i);                                      // Assignment Eb/N0 value
		ofstream fout("BLER.txt",ios::app);                                // Define the name of file
		cout << "Eb/N0 is " << EbN0dB(i) << endl;                         // Show the current Eb/N0 value
		while ( mCounter < mNumOfSimulationLoop )
		{
			if ( Tranceiver.Run() == 0 )                                  // The result of Tranceiver.Run() denotes the RV. When RV is 0, new block is transmitted
			{
				mNumberOfRetr[mCounter] = Tranceiver.mNumberOfRetr;       // Update the Number of Retransmission from the object tranceiver
				mNumberOfError = Tranceiver.mNumberOfError;               // Update the Number of Block Error from the object tranceiver
				mCounter++;
				if ( mNumberOfError == 100 ) break;                       // Break when block errors achieved 100
			}
		}
		block_error_rate(i) = (mNumberOfError+0.0)/mCounter;              // Caculate the BLER of current Eb/N0 value
		fout << Tranceiver.EbN0dB << "  " << block_error_rate(i) << endl; // Write to the file
		fout.close();
		cout << "\t Number of Block Error is " << mNumberOfError << endl; // Number of Block Error of current Eb/N0 value
		cout << "\t Number of Simulation Loop is " << mCounter << endl;   // Number of Simulation Loop
		cout << "\t BLER is " << block_error_rate(i) << endl;             // BLER of current Eb/N0 value
	}
	//Print the results:
  cout << endl;
  cout << "EbN0dB = " <<     EbN0dB       << " [dB]" << endl;
  cout << "BLER = "   << block_error_rate << endl;
  system("pause");
}