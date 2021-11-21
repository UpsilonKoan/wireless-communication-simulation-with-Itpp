//----- ReferenceSignal.h -----
// Programmed by LTE Simulation Group of IMC, SWJTU
#include "stdafx.h"
#include "itpp/itcomm.h"

using std::cin;
using std::cout;
using std::endl;
using namespace itpp;
using namespace std;

//************************************* Class CReferenceSignal Declaration *************************************************************

class CReferenceSignal
{
			/***** Data Members *****/
private:
	cvec  mRSignal;              // Reference signals  (Uplink: SRS, DMRS; DL: Common RS, Dedicated RS)
	//int mGroupnumber;          // The group number u of reference sequence (Defined by users)
	//int mLengthsequence;       // Length of reference sequence to determine m,v=0 when  ,v=0»ò1£¬when  (Define the RS scheme)
		  /***** Function Members *****/
public:
	CReferenceSignal();
	virtual ~ CReferenceSignal();
	int GenURSSequence();        // Generate uplink reference signal sequence
	int GenDRSSequence();        // generate downlink reference signal sequence
	int DemodRSmap();            // Demodulation reference signal map to resource element
    int SoundingRSmap();         // Sounding reference signal map to resource element
    int CommonRSmap();           // Common reference signal map to resource element
    int DedicatedRSmap();        // Dedicated reference signal map to resource element
    int GetDemodRS();            // Get demodulation reference signal from resource element
    int GetSoundingRS();         // Get Sounding reference signal from resource element
    int GetCommonRS();           // Get common reference signal from resource element
    int GetDedicatedRS();        // Get dedicated reference signal from resource element  
    int Run();                   // Run all functions of CReferenceSignal
    int  DRun();                 // Run all functions of CReferenceSignal
};
