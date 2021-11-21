//----- TranChannel.h -----
// Programmed by LTE Simulation Group of IMC, SWJTU
#ifndef TranChannel_HPP
#define TranChannel_HPP

#include <iostream>
#include <string>
#include <fstream>
#include "itpp/itcomm.h"

using std::cin;
using std::cout;
using std::endl;
using namespace itpp;
using namespace std;

//************************************* Class CTranChannel Declaration *************************************************************
class CTranChannel
{
  /***** Data Members *****/
protected:
  imat mSpace;
  imat mSpace1;
  imat mTempDataSpace;       // Temporary Data Space
  mat mTempDataSpace1;       // Temporary Data Space
  imat mBufferSpace;         // Buffer RateMatch Data Space
  ivec mBufferSpace1;        // Buffer RateMatch Data Space
  ivec mTurboIntraInterlSeq; // Interleaver sequence in Turbo encoding
  imat SelectBitNum;         // The number of Selected bits
  mat mReceiverMemory;       // Buffer in Receiver, store the RVs
  long mTTI;                 // Transmission time interval
  long mTranBlockNum;        // Number of transport blocks in a TTI
  long mBitsPerTTI;          // Bits of every TTI
  long mCrcLen;              // Lengths of CRC attachment of transport blocks
  long mCrcCodBlock;         // Lengths of code blocks
  int mCodBlockNum;          // Number of code blocks
  int mLongCodBlockNum;      // Number of long code blocks
  int mShortCodBlockNum;     // Number of short code blocks
  int mLongCodBlockLen;      // Lengths of long code block
  int mShortCodBlockLen;     // Lengths of short code block
  int mCrcFillBitNum;        // Number of fill bits in CRC attachment
  int mLongCollectBitLen;    // Lengths of long code block after CollectBit
  int mShortCollectBitLen;   // Lengths of short code block after CollectBit
  int mLongRateMatchLen;     // Lengths of long code block after RateMatch
  int mShortRateMatchLen;    // Lengths of short code block after RateMatch
  int mTrChEncodeType;       // Type of transport encoder, 0—Turbo; 1—Convolutional
  int mpRateMatchCharacter;  // Character of rate match, 0—Turbo; 1—Convolutional
  int mSoftbuffSize;         // The size of soft buffer
  int mQm;                   // Modulation order
  int mNl;                   // Transmission layer order
  int mUpDownflag;           // Flag to differ uplink with downlink

public:
   vec  mpTrChInputData;        // Input data of transport channel
   bvec mpTrChOutputData;       // Output data of transport channel
   bvec mUserData;              // Store the User Data which is one of the provided precondition of judgement()
   bool mCorrectness;           // The final result of crc check, Notice that if and only if all the parity check is OK, mCorrectness=1 
   int  mRedunVersionNum;       // The redundancy version number
   int mTranBlockSize;          // The size of every transport block
   double mCodeRate ;           // Rate of Turbo coding, constant as 1/3

		/***** Function Members *****/
 public:
   CTranChannel();
   /*------------------------------------------------------------------------
     Construct a TranChannel object.
     Precondition:  None.
     Postcondition: An empty TranChannel object has been constructed
   ------------------------------------------------------------------------*/
   virtual ~CTranChannel();
   /*------------------------------------------------------------------------
     Class destructor.
     Precondition:  None.
     Postcondition: An TranChannel object has been destructed.
   ------------------------------------------------------------------------*/
   void init_turbo_list(string filename,imat &a);
   /*------------------------------------------------------------------------
     Initialize the Inner interleaver in turbo coding.
     Precondition:  File turbo_list.txt is provided.
     Postcondition: The parameters are readed from the file and then stored
	                to imat a.
   ------------------------------------------------------------------------*/
   void AttachCrc();
   /*------------------------------------------------------------------------
     CRC coding, attach the CRC parity bits.
	 Precondition:  Many.
     Postcondition: The data need to be transmitted is generated and then
	                CRC coded.
   ------------------------------------------------------------------------*/
   void EncodTrCh(int Counter);
   /*------------------------------------------------------------------------
     Encode the transport channel.
     Precondition:  The counter of segment block is provided, mTrChEncodeType
	                is given to choose the Transport channel coding type.
     Postcondition: The segment is coded by Turbo if (mTrChEncodeType==0)
	                else by Tail Biting Convolutional encoding.
   ------------------------------------------------------------------------*/
   void CodOfTailBiting(int Counter);
   /*------------------------------------------------------------------------
     Tail Biting Covolutional encoding.
     Precondition:  The counter of segment block is provided.
     Postcondition: The segment is encoded by Tail Biting Convolutional.
   ------------------------------------------------------------------------*/
   void CodOfTurbo(int Counter);
   /*------------------------------------------------------------------------
     Tail Biting Covolutional encoding.
     Precondition:  The counter of segment block is provided.
     Postcondition: The segment is encoded by Tail Biting Convolutional.
   ------------------------------------------------------------------------*/
   ivec TurboIntraInterl(int length);
   /*------------------------------------------------------------------------
     Set up the inner interleaver sequence in Turbo encoding.
     Precondition:  The lenth is provided and init_turbo_list() works properly.
     Postcondition: The interleaver sequence is obtained.
   ------------------------------------------------------------------------*/
   void MatchRate(int Counter);
   /*------------------------------------------------------------------------
     Rate Matching.
     Precondition:  SubBlockInterl(), CollectBit() and SelectBit() works
	                properly.
     Postcondition: The segment which is Counter number is rate matched.
   ------------------------------------------------------------------------*/
   void SubBlockInterl();
   /*------------------------------------------------------------------------
     Sub-block Interleaver.
     Precondition:  None.
     Postcondition: Each sub-block is interleaved.
   ------------------------------------------------------------------------*/
   void CollectBit(int Counter);
   /*------------------------------------------------------------------------
     Collecting bits.
     Precondition:  None.
     Postcondition: Systematic and Redundant Bits are collected.
   ------------------------------------------------------------------------*/
   void SelectBit(int Counter);
   /*------------------------------------------------------------------------
     Selecting and pruning bits.
     Precondition:  Systematic and Redundant Bits are collected.
     Postcondition: The bits are selected according to RV.
   ------------------------------------------------------------------------*/
   void ConcatCodeBlock(imat Input);
   /*------------------------------------------------------------------------
     Concatenated coded data block.
     Precondition:  All the CB blocks complete the encoding process.
     Postcondition: The CB blocks are concated.
   ------------------------------------------------------------------------*/
   void SeparCodeBlock();                    // Inverse of Concatenated coded data block
   void DeMatchRate(int Counter);            // Inverse of Rate matching
   void FillBit(int Counter);                // Inverse of selecting and pruning bits
   void SeparateBit(int Counter);            // Inverse of Collecting bits
   void DeSubBlockInterl(int Counter);       // Inverse of Sub-block Interleaver
   void DecodTrCh(int Counter);              // Decoding transmitting channel
   void DecodOfTailBiting(int Counter);      // Tail Biting Convolutional decoding
   void DecodOfTurbo(int Counter);           // Turbo decoding
   void CrcCheck();                          // CRC check
   void DeAttachCrc();                       // Inverse of CRC attchment
   void DeSegCodBlock();                     // Inverse of Segment coding blocks  
   void DeSegAttachCrc(int Counter);         // Inverse of CRC check of Segment coding blocks
};

//************************************* Class CDLSCh Declaration *************************************************************
class CDLSCh: public CTranChannel
{
private:
	int mCodBlockLoopCounter;                // Counter flag

public:
	CDLSCh();
	~CDLSCh();
    void InniCTranChannel(int SoftbuffSize);
    void Run();                             // Run all coding functions of CDLSch
    void DRun();                            // Run all decoding functions of CDLSch
    void SetCDLSChDataForPhysicalCh();      // Transmit data to Physical channel
    void OutCDLSChData();                   // Output the decoding data 
};

//************************************* Class CULSCh Declaration *************************************************************
class CULSCh: public CTranChannel
{
private:
	int mCodBlockLoopCounter;               // Counter flag

public:
	CULSCh();
    ~CULSCh();
    void InniCTranChannel(int SoftbuffSize);
 //  void InniCTranChannel(struct PhyChannel para);
    void Run();                             // Run all coding functions of CDLSch
    void DRun();                            // Run all decoding functions of CDLSch
    void SetCULSChDataForPhysicalCh();      // Transmit data to Physical channel
    void OutCULSChData();                   // Output the decoding data
};

//************************************* Class CUCI Declaration *************************************************************
class CUCI: public CTranChannel
{
private:
	cvec mUCQIOutputData;                  // Output data of CQI coding
	cvec mURIOutputData;                   // Output data of RI coding
	cvec mUACKOutputData;                  // Output data of ACK coding
	int  mConInforEncodeType;              // Type of control information encoder

public:
	CUCI ();
	virtual ~ CUCI ();
	void Run();                            // Run all encoding functions of CUCI
	void DRun();                           // Run all decoding functions of CUCI
	public:
	void SetUCIDataForPhysicalCh();        // Transmit data to Physical channel
	void OutUCIData();                     // Output the data 

public:
	void PNGen();                         // Generate PN sequence for UCI
	void CUCI::cqi_code_list(string filename,imat &a);
	void EncodforCQI();                   // Channel coding for CQI
	void EncodforRI();                    // Channel coding for RI
	void EncodforACK();                   // Channel coding for ACK
	void DecodforCQI();                   // Channel decoding for CQI
	void DecodforRI();                    // Channel decoding for RI
	void DecodforACK();                   // Channel decoding for ACK
};

//************************************* Class CDCI Declaration *************************************************************
class CDCI: public CTranChannel
{
public:
	CDCI();
	virtual ~ CDCI ();
	void Run();                           // Run all coding functions of CDCI
	void DRun();                          // Run all decoding functions of CDCI
	public:
	void SetDCIDataForPhysicalCh();       // Transmit data to Physical channel
	void OutDCIData();                    // Output the data 

protected:
	void PNGen ();                        // Generate PN sequence for DCI
};


#endif


