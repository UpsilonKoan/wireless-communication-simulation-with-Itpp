//----- VectorChannel.h -----
// Programmed by LTE Simulation Group of IMC, SWJTU
#include "itpp/itcomm.h"
using namespace itpp;

//************************************* Class CTranceiver Declaration *************************************************************
class CVectorChannel             //-----------MIMO channel
{
				/***** Data Members *****/
public:
	int mModelCase;                // Modified Pedestrian A, Vehicular A, Pedestrian B, Single Path{00,01,10,11}
	int mNumOfPath;                // Number of Paths
	int NumOfSample;               // Number of Sampled
	int mNumOfMSAntenna;           // MS antenna Number
	int mNumOfBSAntenna;           // BS antenna Number
	double mFrequencyOfCarrier;    // Frequency of Carrier
	double mDoppler;               // Maximum of Doppler shift
	float mAntennaSpacingBS;       // 
	float mAntennaSpacingMS;       //
	double AOAPerPathMS;           // 22.5 or 67.5or-67.5
	double AODPerPathBS;           // 50  20
	vec mAngleSpreadPerPathMS;     // ChannelCoefficient is calculated by path
	vec mAngleSpreadPerPathBS;     // ChannelCoefficient is calculated by path
	cvec mVelocityOfMS;            // The velocity of MS
	//----- Hold In Abeyance -----
	//vec mPathDelays(6);
	//Complex *mChannelCoefficients;

        /***** Function Members *****/
	    /***** Constructors *****/
public:
	CVectorChannel();
    /*------------------------------------------------------------------------
      Construct a VectorChannel object.
      Precondition:  None.
      Postcondition: An empty VectorChannel object has been constructed
	------------------------------------------------------------------------*/
	virtual ~CVectorChannel();
	/*-----------------------------------------------------------------------
      Class destructor 
	  Precondition:  None
	  Postcondition: The VectorChannel object has been destructed.
    ------------------------------------------------------------------------*/
	void SetCVectorChannelParameter();
	/*-----------------------------------------------------------------------
      Set up the parameters of Function Member of Vector Channel
	  //----- Hold In Abeyance -----
	  Precondition:  None
	  Postcondition: None
    ------------------------------------------------------------------------*/
	mat ChannelCorrelation(double AngleSpread,double Angle,double RatioDistance,int NumOfAntenna);
	/*-----------------------------------------------------------------------
      Caculate the Correlation of the Channel
	  Precondition:  AngleSpread, Angle, RatioDistance and NumOfAntenna are
	                 given
	  Postcondition: Correlation of the Channel is returned
    ------------------------------------------------------------------------*/
	cvec ScalarChannel(cvec velocity, double fre, double doppler);
	/*-----------------------------------------------------------------------
      Caculate the scalar coeffiencts of the channel
	  Precondition:  velocity, frequency, maximum of doppler shift are given
	  Postcondition: Channel coefficiencts are returned
    ------------------------------------------------------------------------*/
	cmat RunScalarChannel();
	/*-----------------------------------------------------------------------
      All channel coefficiencts are caculated and stored to variable scalarch
	  Precondition:  ScalarChannel() works properly, total number of MS and
	                 BS antenna are given
	  Postcondition: Final channel coefficiencts are returned
    ------------------------------------------------------------------------*/
	cmat TransmissionCoefficientPerPath(int path);
	/*-----------------------------------------------------------------------
      Caculate the Transmission Coefficient of every path
	  Precondition:  None
	  Postcondition: coefficientp is returned which is for RunCVectorChannel()
    ------------------------------------------------------------------------*/
	cmat RunCVectorChannel();
	/*-----------------------------------------------------------------------
      Run all the functions of  CVectorChannel
	  Precondition:   Every function member works properly
	  Postcondition:  Return neccessary information of MIMO channel coefficient
    ------------------------------------------------------------------------*/

	//----- Hold In Abeyance -----
	//cmat TotalTransmissionCoefficient(); 
	//SetVectorParameter();
	//RunCVectorChannel();
};


