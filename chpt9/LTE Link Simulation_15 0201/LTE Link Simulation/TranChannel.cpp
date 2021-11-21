//----- TranChannel.cpp -----
// Programmed by LTE Simulation Group of IMC, SWJTU
#include "stdafx.h"
#include "itpp/itcomm.h"
#include <string>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "TranChannel.h"

using std::cin;
using std::cout;
using std::endl;
using namespace itpp;

# define Rate 0.5 // Coding Rate after Rate matching, 1/3 1/2 2/3 3/4 4/5, here is 1/3

//************************************* Class CTranChannel Defination *************************************************************

//--- Definition of CTranChannel constructor
CTranChannel::CTranChannel()
{
    mTranBlockSize=int (240/Rate);
	mTrChEncodeType=0;
	mQm=2;
	mNl=1;
	mUpDownflag=0;
	mCodeRate=(double)1/3;// Turbo Coding Rate,constant
}

//--- Definition of CTranChannel destructor
CTranChannel::~CTranChannel()
{
}

//--- Definition of init_turbo_list()
void CTranChannel::init_turbo_list(string filename,imat &a)
{
    
    ifstream fin(filename.c_str());
    if(!fin) 
	  {cout<<"can not open turbo_list.txt"<<endl;return ;}

	int a1,a2,a3,a4,i=0;
 
    for(int i=0;i<a.rows();i++)
	{
	  fin>>a1>>a2>>a3>>a4;
      {
		a(i,0)=a1;
		a(i,1)=a2;
		a(i,2)=a3;
		a(i,3)=a4;	
      }
	}
    fin.close();

}

//--- Definition of AttachCrc()
void CTranChannel::AttachCrc()       //CRC attchment
{
   //Local variables:
   int i, ii, k, s, r, KK, Bp, Kr;
   int UserDataLen;
   int CodBlockNum, SegCrcLen, LenAfterCrc, MixCodBlockLen=6144;
   int LongCodBlockLen, ShortCodBlockLen, LongCodBlockNum, ShortCodBlockNum, deltak, FillBitNum;
   bvec bits,coded_bits, segbits,coded_segbits,decoded_bits,CrcOutput;
   imat k_list;

   //Declarations of the CRC class:
   CRC_Code crc(std::string("CRC-24"));// Set CRC code to one of the standardpolynomials using the string value "CRC-24"

   // Initializations:
   UserDataLen = 240;
   mCrcLen = 24;
   k_list.set_size(188, 4, false);
   bits = randb(UserDataLen);     // Generate random bits which is transmitted
   mUserData = bits;              // Stored for comparison in judgment()

   // CRC attchment
   coded_bits = crc.encode(bits);   
   LenAfterCrc = coded_bits.length();

   // Compare with the maximum code block size to determine whether the code is divided into small pieces
   if (LenAfterCrc <= MixCodBlockLen)
   {
		SegCrcLen = 0;
        CodBlockNum = 1;
        Bp = LenAfterCrc;
    }	
  else
   {
	   SegCrcLen = 24;
       CodBlockNum = (int) ceil( (double) LenAfterCrc / (MixCodBlockLen - SegCrcLen) );
       Bp = LenAfterCrc + CodBlockNum * SegCrcLen;
   }
  
  
  // Read the parameters from turbo_list.txt to turbo_para: 
  init_turbo_list("turbo_list.txt",k_list);

  // Set the length of number of the long/short code block

   // First segmentation size
  for (i = 0; i < 188; i++)
  {
	  KK = k_list(i,1);
      if (CodBlockNum * KK >= Bp)
      {
		  LongCodBlockLen = KK;
          break;
	  }
  }
 
 
  // Determine whether there is one code block 
  if (CodBlockNum == 1)
  {
	  LongCodBlockNum = 1;
	  ShortCodBlockNum = 0;
	  ShortCodBlockLen = 0;
  }
  else 
	  // Second segmentation size
	  if(CodBlockNum > 1)  
	  { 
		  for (i = 0; i < 188; i++)
		  {   
			  KK = k_list(i,1);
			  if (KK == LongCodBlockLen)
			  {
				  ShortCodBlockLen = k_list(i-1,1);
				  break;
			  }
		  }
		  deltak = LongCodBlockLen - ShortCodBlockLen;
		  ShortCodBlockNum = (int) floor( (double) (CodBlockNum * LongCodBlockLen - Bp) / deltak);
		  LongCodBlockNum = CodBlockNum - ShortCodBlockNum;
	  }

   
   // Number of filler bits:
   FillBitNum = (int)LongCodBlockNum * LongCodBlockLen + ShortCodBlockNum * ShortCodBlockLen - Bp;
   

   //Insertion of filler bits
   CrcOutput.set_size(FillBitNum + LenAfterCrc, false);
   mSpace.set_size(CodBlockNum, LongCodBlockLen, false);
  
   // Set the filler bits to NULL:
   for (k = 0;k < FillBitNum; k++)
   {
	   CrcOutput(k) = NULL; 
   }
   
   k = FillBitNum;

   for (i = 0; i < LenAfterCrc; i++)
   {
	   CrcOutput(k) = coded_bits(i); 
	   k = k + 1;
   }
 
   s = 0;
   for (r = 0; r< CodBlockNum; r++)
   {
	   if(r < ShortCodBlockNum)
	   {
		   Kr = ShortCodBlockLen;
	   }
	   else
	   {
		   Kr = LongCodBlockLen;
	   }

	   if (CodBlockNum > 1)
	   { 
		   segbits = CrcOutput(s, Kr-SegCrcLen-1+s);
		   s = s + Kr - SegCrcLen;

		   // Attach an additional CRC sequence of L = 24 bits to each code block: 
		   coded_segbits= crc.encode(segbits);

		   // Store the bits in the public space named mSpace:
		   for(int ii=0; ii < Kr; ii++)
			   mSpace(r,ii) = coded_segbits(ii);
	   }
	   else
	   {
		   for(ii = 0; ii < Kr; ii++)
			   mSpace(r,ii) = CrcOutput(ii);
	   }
   }//end for (r = 0; r< CodBlockNum; r++)

   // Set the value of member variable:
   mCodBlockNum = CodBlockNum;
   mLongCodBlockNum = LongCodBlockNum;
   mShortCodBlockNum = ShortCodBlockNum;   
   mLongCodBlockLen = LongCodBlockLen;
   mShortCodBlockLen = ShortCodBlockLen;
   mCrcFillBitNum = FillBitNum;

}

//--- Definition of EncodTrCh()
void CTranChannel::EncodTrCh(int Counter)
{
  if(mTrChEncodeType==0)
	  CodOfTurbo(Counter);
  else
	  if(mTrChEncodeType==1)
 	    CodOfTailBiting(Counter);
}

//--- Definition of CodOfTailBiting()
void CTranChannel::CodOfTailBiting(int Counter)   //Tail Biting Convolutional encoding
{
   //Local variables:
   int i, j;
   ivec generator(3), get_bits, encoded_bits_ivec;
   bvec transmitted_bits, encoded_bits, decoded_bits, error_bits;
   vec tx_signal, rx_signal;

   //Declarations of the Convolutional Code class:
   BPSK bpsk;
   Convolutional_Code code;
   
   // Initializations:
   generator(0)=0133;
   generator(1)=0165;
   generator(2)=0171;
   code.set_generator_polynomials(generator,7);

   // Determine the number of code blocks and code block length:
   if(mCodBlockNum == 1)
   {
	   get_bits.set_size(mLongCodBlockLen,false);	   
   }
   else      
	   if(Counter<mShortCodBlockNum)
	    {
	      get_bits.set_size(mShortCodBlockLen,false);	      
	    }
       else
	   {
	      get_bits.set_size(mLongCodBlockLen,false);
	   }
     
   for(j = 0; j < get_bits.length(); j++)
   {
	   get_bits(j )= mSpace(Counter,j);
   }
   if(Counter == 0)
	   for(i=0;i<mCrcFillBitNum;i++)
	   {
		   get_bits(i) = 0;
	   }
   transmitted_bits.set_size(length(get_bits),false);
   transmitted_bits=to_bvec(get_bits);
   code.encode_tail(transmitted_bits, encoded_bits);
   encoded_bits_ivec=to_ivec(encoded_bits);
   mTempDataSpace.set_size(3,length(encoded_bits)/3,false);
   mTempDataSpace.clear();
   for(i=0,j=0;i<length(encoded_bits);i=i+3,j++)
   {
	   mTempDataSpace(0,j)=encoded_bits(i);
       mTempDataSpace(1,j)=encoded_bits(i+1);
       mTempDataSpace(2,j)=encoded_bits(i+2);
    } 
}

//--- Definition of CodOfTurbo()
void CTranChannel::CodOfTurbo(int Counter)    //Turbo encoding
{
   //Local variables:
   int i, j;
   int UncodedLen, constraint_length;
   ivec gen(2), encoded_bits1, get_bits, interleaver_sequence;
   bvec transmitted_bits, encoded_bits, in1, in2, d0, d1, d2;
   bmat parity1, parity2;  

   //Declarations of the turbo class:
   Turbo_Codec turbo; 

  // Determine the number of code blocks and code block length:
   if(mCodBlockNum == 1)
   {
	   get_bits.set_size(mLongCodBlockLen,false);
   }
   else
   {
	   if(Counter < mShortCodBlockNum)
	    
	      get_bits.set_size(mShortCodBlockLen,false);
	    
       else
	      get_bits.set_size(mLongCodBlockLen,false);
   }
   
   // Initializations:
   UncodedLen = get_bits.length();
   constraint_length = 4;
   gen(0) = 015;
   gen(1) = 013;
   transmitted_bits.set_size(UncodedLen,false);
   d0.set_size(UncodedLen,false);
   d1.set_size(UncodedLen,false);
   d2.set_size(UncodedLen,false);

   // Get the uncoded bits from the public space named mSpace:
   for(i = 0; i < UncodedLen; i++)
	   get_bits(i) = mSpace(Counter,i);   //Store the data into mSpace after segmentation
  
   // Set the filler bits to zero:
   if(Counter == 0)
	   for(i = 0; i < mCrcFillBitNum; i++)
		   get_bits(i) = 0;

   // Convert an ivec to a bvec:
   transmitted_bits = to_bvec(get_bits);
   
   // Call the function TurboIntraInterl() to generated interleaver sequence index
   mTurboIntraInterlSeq = TurboIntraInterl(UncodedLen);
   
   // Set the Turbo encoder/decoder parameters:
   turbo.set_parameters(gen, gen, constraint_length, mTurboIntraInterlSeq, 8, "LOGMAX", 0.75);

   // Encode the code block:
    turbo.encode(transmitted_bits,encoded_bits);
  
  // Store the encoded_bits in the public space named mTempDataSpace:
   mTempDataSpace.set_size(3,length(encoded_bits)/3,false);
   mTempDataSpace.clear();
   for(i=0,j=0;i<length(encoded_bits);i=i+3,j++)
	{   if(Counter==0 && j<mCrcFillBitNum)
       {
		mTempDataSpace(0,j)=2;
        mTempDataSpace(1,j)=2;
       }
	else
	   {
        mTempDataSpace(0,j)=encoded_bits(i);
        mTempDataSpace(1,j)=encoded_bits(i+1);
       }
        mTempDataSpace(2,j)=encoded_bits(i+2);
     }

}

//--- Definition of TurboIntraInterl()
ivec CTranChannel::TurboIntraInterl(int length)   //Interleaver in Turbo encoding
{   
	int i,j;
	int K,f1,f2; 
	ivec index;
	imat turbo_para;

	// Initializations:
	K = length;	
	index.set_size(K,false);
    turbo_para.set_size(188,4,false);
	turbo_para.clear();

	// Read the parameters from turbo_list.txt to turbo_para:
    init_turbo_list("turbo_list.txt",turbo_para);
	for(j=0;j<188;j++)
	{
		if(K == turbo_para(j,1))
		{
			f1 = turbo_para(j,2);
            f2 = turbo_para(j,3);
			break;
		}
	}

	// Calculate the interleaver sequence index:
	for(i = 0; i < K; i++)
	{      
		index(i)=(i * ((f1 + f2 * i) % K)) % K;//p95=> index(i)=(f1i+f2i^2) mod K
	}

    return index;
	
}

//--- Definition of MatchRate()
void CTranChannel::MatchRate(int Counter)
{
	
	SubBlockInterl();
	CollectBit(Counter);
    if(Counter<mShortCodBlockNum)
		mShortCollectBitLen=mTempDataSpace.cols();
	else
		mLongCollectBitLen=mTempDataSpace.cols();
	SelectBit(Counter);
}

//--- Definition of SubBlockInterl()
void CTranChannel::SubBlockInterl()  //Sub-block Interleaver
{
    int column, row, d, nd;
	ivec d0,d1,d2;
	ivec v0,v1,v2;
	imat y0,y1,y2;
	imat y0_con,y1_con,y2_con;
	int i,ii;
	int j=mTempDataSpace.cols();

    d0.set_size(mTempDataSpace.cols(),false);
    d1.set_size(mTempDataSpace.cols(),false);
    d2.set_size(mTempDataSpace.cols(),false);

    for(i=0;i<mTempDataSpace.cols();i++)
	{
		d0(i)=mTempDataSpace(0,i);
        d1(i)=mTempDataSpace(1,i);
        d2(i)=mTempDataSpace(2,i);
	}

	d=d0.length();
	column=32;
	row=(int)ceil(double(d)/double(column));
	nd=row*column-d;

	y0.set_size(row,column,false);y1.set_size(row,column,false);y2.set_size(row,column,false);
	y0_con.set_size(row,column,false);y1_con.set_size(row,column,false);y2_con.set_size(row,column,false);
	v0.set_size(row*column,false);v1.set_size(row*column,false);v2.set_size(row*column,false);
	mTempDataSpace.set_size(3,row*column,false);
	mTempDataSpace.clear();
    
	//mTempDataSpace.clear();
    y0.clear();y1.clear();y2.clear();
    y0_con.clear();y1_con.clear();y2_con.clear();
    v0.clear();v1.clear();v2.clear();

	for(i=0;i<row;i++)
	{
	    for(ii=0;ii<column;ii++)
		if((ii<nd)&&(i==0))
		{y0(i,ii)=2;y1(i,ii)=2;y2(i,ii)=2;}                                                     //pad dummy bits
		else 
		{y0(i,ii)=d0(i*column+ii-nd);y1(i,ii)=d1(i*column+ii-nd);y2(i,ii)=d2(i*column+ii-nd);}  //write the input bit sequence row by row
	}

	if(mTrChEncodeType==0)
	{
	   /********** Perform the inter-column permutaion for d0 and d1 **********/
	   y0.swap_cols(1,16);y0.swap_cols(2,8);y0.swap_cols(3,24);y0.swap_cols(5,20);y0.swap_cols(6,12);y0.swap_cols(7,28);
	   y0.swap_cols(9,18);y0.swap_cols(11,26);y0.swap_cols(13,22);y0.swap_cols(15,30);y0.swap_cols(19,25);y0.swap_cols(23,29);

	   y1.swap_cols(1,16);y1.swap_cols(2,8);y1.swap_cols(3,24);y1.swap_cols(5,20);y1.swap_cols(6,12);y1.swap_cols(7,28);
	   y1.swap_cols(9,18);y1.swap_cols(11,26);y1.swap_cols(13,22);y1.swap_cols(15,30);y1.swap_cols(19,25);y1.swap_cols(23,29);

	   /********** Output of the block interleaver for d0 and d1 column by column **********/
       for(i=0;i<column;i++)
	   {
		   v0.set_subvector(i*row,row+i*row-1,y0.get_col(i));
		   v1.set_subvector(i*row,row+i*row-1,y1.get_col(i));
	   }
	   
	   /********* Output of the sub-block interleaver for d2 **********/
	   for(i=0;i<column-1;i++)
		   v2.set_subvector(i*row,row+i*row-1,y2.get_col(i+1));
	   for(i=0;i<row;i++)
		   v2(i+row*(column-1))=y2(mod((i+1),row),0);
	}
	else 
	{
		/********** Perform the inter-column permutation patter for d0, d1 and d2 ***********/
		y0_con.set_col(0,y0.get_col(1));y0_con.set_col(1,y0.get_col(17));y0_con.set_col(2,y0.get_col(9));y0_con.set_col(3,y0.get_col(25));
		y0_con.set_col(4,y0.get_col(5));y0_con.set_col(5,y0.get_col(21));y0_con.set_col(6,y0.get_col(13));y0_con.set_col(7,y0.get_col(29));
		y0_con.set_col(8,y0.get_col(3));y0_con.set_col(9,y0.get_col(19));y0_con.set_col(10,y0.get_col(11));y0_con.set_col(11,y0.get_col(27));
		y0_con.set_col(12,y0.get_col(1));y0_con.set_col(13,y0.get_col(23));y0_con.set_col(14,y0.get_col(15));y0_con.set_col(15,y0.get_col(31));
		y0_con.set_col(16,y0.get_col(0));y0_con.set_col(17,y0.get_col(16));y0_con.set_col(18,y0.get_col(8));y0_con.set_col(19,y0.get_col(24));
		y0_con.set_col(20,y0.get_col(4));y0_con.set_col(21,y0.get_col(20));y0_con.set_col(22,y0.get_col(12));y0_con.set_col(23,y0.get_col(28));
		y0_con.set_col(24,y0.get_col(2));y0_con.set_col(25,y0.get_col(18));y0_con.set_col(26,y0.get_col(10));y0_con.set_col(27,y0.get_col(26));
        y0_con.set_col(28,y0.get_col(6));y0_con.set_col(29,y0.get_col(22));y0_con.set_col(30,y0.get_col(14));y0_con.set_col(31,y0.get_col(30));

		y1_con.set_col(0,y1.get_col(1));y1_con.set_col(1,y1.get_col(17));y1_con.set_col(2,y1.get_col(9));y1_con.set_col(3,y1.get_col(25));
		y1_con.set_col(4,y1.get_col(5));y1_con.set_col(5,y1.get_col(21));y1_con.set_col(6,y1.get_col(13));y1_con.set_col(7,y1.get_col(29));
		y1_con.set_col(8,y1.get_col(3));y1_con.set_col(9,y1.get_col(19));y1_con.set_col(10,y1.get_col(11));y1_con.set_col(11,y1.get_col(27));
		y1_con.set_col(12,y1.get_col(1));y1_con.set_col(13,y1.get_col(23));y1_con.set_col(14,y1.get_col(15));y1_con.set_col(15,y1.get_col(31));
		y1_con.set_col(16,y1.get_col(0));y1_con.set_col(17,y1.get_col(16));y1_con.set_col(18,y1.get_col(8));y1_con.set_col(19,y1.get_col(24));
		y1_con.set_col(20,y1.get_col(4));y1_con.set_col(21,y1.get_col(20));y1_con.set_col(22,y1.get_col(12));y1_con.set_col(23,y1.get_col(28));
		y1_con.set_col(24,y1.get_col(2));y1_con.set_col(25,y1.get_col(18));y1_con.set_col(26,y1.get_col(10));y1_con.set_col(27,y1.get_col(26));
        y1_con.set_col(28,y1.get_col(6));y1_con.set_col(29,y1.get_col(22));y1_con.set_col(30,y1.get_col(14));y1_con.set_col(31,y1.get_col(30));

		y2_con.set_col(0,y2.get_col(1));y2_con.set_col(1,y2.get_col(17));y2_con.set_col(2,y2.get_col(9));y2_con.set_col(3,y2.get_col(25));
		y2_con.set_col(4,y2.get_col(5));y2_con.set_col(5,y2.get_col(21));y2_con.set_col(6,y2.get_col(13));y2_con.set_col(7,y2.get_col(29));
		y2_con.set_col(8,y2.get_col(3));y2_con.set_col(9,y2.get_col(19));y2_con.set_col(10,y2.get_col(11));y2_con.set_col(11,y2.get_col(27));
		y2_con.set_col(12,y2.get_col(1));y2_con.set_col(13,y2.get_col(23));y2_con.set_col(14,y2.get_col(15));y2_con.set_col(15,y2.get_col(31));
		y2_con.set_col(16,y2.get_col(0));y2_con.set_col(17,y2.get_col(16));y2_con.set_col(18,y2.get_col(8));y2_con.set_col(19,y2.get_col(24));
		y2_con.set_col(20,y2.get_col(4));y2_con.set_col(21,y2.get_col(20));y2_con.set_col(22,y2.get_col(12));y2_con.set_col(23,y2.get_col(28));
		y2_con.set_col(24,y2.get_col(2));y2_con.set_col(25,y2.get_col(18));y2_con.set_col(26,y2.get_col(10));y2_con.set_col(27,y2.get_col(26));
        y2_con.set_col(28,y2.get_col(6));y2_con.set_col(29,y2.get_col(22));y2_con.set_col(30,y2.get_col(14));y2_con.set_col(31,y2.get_col(30));

		/***********  The output of the block interleaver ************/
		for(i=0;i<column;i++)
		{
			v0.set_subvector(i*row,row-1+i*row,y0_con.get_col(i));
			v1.set_subvector(i*row,row-1+i*row,y1_con.get_col(i));
			v2.set_subvector(i*row,row-1+i*row,y2_con.get_col(i));
		}
	}

    for(i=0;i<row*column;i++)
	{mTempDataSpace(0,i)=v0(i);mTempDataSpace(1,i)=v1(i);mTempDataSpace(2,i)=v2(i);}

}

//--- Definition of CollectBit()
void CTranChannel::CollectBit(int Counter)
{
    ivec v0,v1,v2;
	ivec w;
	int i;
	int stream_length = mTempDataSpace.cols();

	v0.set_size(stream_length,false);v1.set_size(stream_length,false);v2.set_size(stream_length,false);
	w.set_size(stream_length*3,false);

	for(i=0;i<stream_length;i++)
	{v0(i)=mTempDataSpace(0,i);v1(i)=mTempDataSpace(1,i);v2(i)=mTempDataSpace(2,i);}
	
	if(mTrChEncodeType==0)
	{
	    for(i=0;i<stream_length;i++)
		{
			w(i)=v0(i);
			w(2*i+stream_length)=v1(i);
			w(2*i+1+stream_length)=v2(i);
		}
	}
	else 
	{
	    for(i=0;i<stream_length;i++)
		{
		    w(i)=v0(i);
			w(i+stream_length)=v1(i);
			w(i+2*stream_length)=v2(i);
		}
	}
    /********* write the output of bits collection function into the shared dataspace *********/	
	mTempDataSpace.set_size(1,w.length(),false);
	//mReceiverMemory.set_length(w.length(),false);
	//mBufferSpace.set_size(w.length(),false);
	//mBufferSpace1.set_size(w.length(),false);
    mTempDataSpace.clear();
	for(i=0;i<w.length();i++)
	{
		mTempDataSpace(0,i)=w(i);
		mBufferSpace(Counter,i)=w(i);
	}

}

//--- Definition of SelectBit()
void CTranChannel::SelectBit(int Counter)   //Selecting and pruning bits
{
    
	int ncb,kw,length;
	ivec w,eout;
	int er;
	int g,r,k,row;
	int i,ii;
	if(mRedunVersionNum!=0)
	{
		if(Counter<mShortCodBlockNum)
		    length=mShortCollectBitLen;		    
		else
			length=mLongCollectBitLen;
        mTempDataSpace.set_size(1,length,false);
        for(i=0;i<mTempDataSpace.cols();i++)
		    mTempDataSpace(0,i)=mBufferSpace(Counter,i);

	}
    /************* Calculate the parametres used and initialize  *************/	
	kw = mTempDataSpace.cols();
	row=kw/3/32;
	g=mTranBlockSize/(mNl*mQm);
    r=mod(g,mCodBlockNum);
	w.set_size(kw,false); 
	w.clear();
	for(i=0;i<kw;i++)
		w(i)=mTempDataSpace(0,i);
	
	if(mUpDownflag==0)    // UpDownflag==0: Uplink transport channel; UpDownflag==1: Downlink transport channel;
		ncb=kw;
	else 
	{
	   if((int)floor((double)mSoftbuffSize/(double)mCodBlockNum)>kw) ncb=kw;
	   else ncb=(int)floor((double)mSoftbuffSize/(double)mCodBlockNum);
	}

	/******** Calculate the length of rate matching output bit sequence*********/
	k=row*(2*(int)ceil((double)ncb/(double)(8*row))*mRedunVersionNum+2);

	if(Counter<=mCodBlockNum-r-1)
	   {er=mNl*mQm*(int)floor((double)g/(double)mCodBlockNum);}
	else 
	   {er=mNl*mQm*(int)ceil((double)g/(double)mCodBlockNum);}  // er, the data length after pruning

	/********* Bit selection and pruning **********/
	eout.set_size(er,false);
	eout.clear();
    i=0;ii=0;
    if(mTrChEncodeType==0)
	{
		while(i<er)
		{
		   if(w(mod(k+ii,ncb))!=2)
		   {eout(i)=w(mod(k+ii,ncb));SelectBitNum(Counter,i)=mod(k+ii,ncb);i+=1;}
		   ii+=1;
		}

	}
	else 
	{
	    while(i<er)
	    {
		   if(w(mod(ii,kw))!=0)
		   {eout(i)=w(mod(ii,kw));SelectBitNum(Counter,i)=ii;i+=1;}
		   ii+=1;
		}
	}
   
	mTempDataSpace.set_size(1,er,true);
	mTempDataSpace.clear();
	for(i=0;i<er;i++)
		mTempDataSpace(0,i)=eout(i);

}

//--- Definition of ConcatCodeBlock()
void CTranChannel::ConcatCodeBlock(imat Input)  //Concatenated coded data block
{
    ivec f;
	int r,g,er;
	int i,ii,iii;
	
	f.set_size(mTranBlockSize,false);
	g=mTranBlockSize/(mNl*mQm);
	r=mod(mTranBlockSize/(mNl*mQm),mCodBlockNum);
	
	i=0;iii=0;
	while(i<mCodBlockNum)
	{
	    if(i<=mCodBlockNum-r-1)
		{er=mNl*mQm*(int)floor((double)g/(double)mCodBlockNum);mShortRateMatchLen=er;}
		else 
		{er=mNl*mQm*(int)ceil((double)g/(double)mCodBlockNum);mLongRateMatchLen=er;}

		ii=0;
		while(ii<er)
		{
		    f(iii)=Input(i,ii);
			ii+=1;
			iii+=1;
		}
		i+=1;/*iii=iii-1;*/
	}

	// The output of transport channel
	mpTrChOutputData.set_size(mTranBlockSize,false);
	mpTrChOutputData.clear();
	for(i=0;i<mTranBlockSize;i++)
		mpTrChOutputData(i)=f(i);

}

//--- Definition of SeparCodeBlock()
void CTranChannel::SeparCodeBlock()  // Inverse of Concatenated coded data block
{
    mat e;
	int length1,length2;
	int g,r;
	int i,ii,iii;
    
	g=mTranBlockSize/(mNl*mQm);
	r=mod(mTranBlockSize/(mNl*mQm),mCodBlockNum);
	length1=mNl*mQm*(int)floor((double)g/(double)mCodBlockNum);
	length2=mNl*mQm*(int)ceil((double)g/(double)mCodBlockNum);

    e.set_size(mCodBlockNum,length2,false);

	i=0;iii=0;
	while(i<mCodBlockNum)
	{
	    ii=0;
		if(i<=mCodBlockNum-r-1)
	    {
		    while(ii<length1)
		    {e(i,ii) = mpTrChInputData(iii);ii+=1;iii+=1;}
	    }
	    else 
	    {
		    while(ii<length2)
		    {e(i,ii)=mpTrChInputData(iii);ii+=1;iii+=1;}
	    }
		i+=1;
	}

	mTempDataSpace1.set_size(mCodBlockNum,length2+1,false);
	mTempDataSpace1.clear();
	for(i=0;i<mCodBlockNum;i++)
		for(ii=0;ii<length2;ii++)
			mTempDataSpace1(i,ii)=e(i,ii);

}

//--- Definition of DeMatchRate()
void CTranChannel::DeMatchRate(int Counter)    // Inverse of Rate matching
{
    FillBit(Counter);
	SeparateBit(Counter);
	DeSubBlockInterl(Counter);
}

//--- Definition of FillBit()
void CTranChannel::FillBit(int Counter)         //Inverse of selecting and pruning bits
{
    ivec w;
	//int er;
	int column,row,kw;
	//int g,r;
	int i,ii,iii;

	column=32;
 //   g=mTranBlockSize/(mNl*mQm);
	//r=mod(mTranBlockSize/(mNl*mQm),mCodBlockNum);

	//if(counter<=mCodBlockNum-r-1)
	//	  er=mNl*mQm*(int)floor((double)g/(double)mCodBlockNum);
	//else er=mNl*mQm*(int)ceil((double)g/(double)mCodBlockNum);

	//if(mUpDownflag==0)    // UpDownflag==0: Uplink transport channel; UpDownflag==1: Downlink transport channel;
	//ncb=kw;
	//else 
	//{
	//   if((int)floor((double)mSoftbuffSize/(double)mCodBlockNum)>kw) ncb=kw;
	//   else ncb=(int)floor((double)mSoftbuffSize/(double)mCodBlockNum);
	//}
	if(mTrChEncodeType==0)
	{
		if(Counter<mShortCodBlockNum)
		{
			row=(int)ceil(double(mShortCodBlockLen+4)/double(column));
			kw=3*row*column;
			w.set_size(kw,false);
			//nd=row*column-mLongCodBlockLen;
			//k=row*(2*(int)ceil((double)ncb/(double)(8*row))*mRedunVersionNum+2);
			ii=0;iii=0;

			for(i=0;i<kw;i++)
			{
				if(i==SelectBitNum(Counter,ii))
				{
					mReceiverMemory(Counter,i)=mTempDataSpace1(Counter,iii);
					/*w(i)=mSpace(counter,iii);*/
					ii+=1;iii+=1;
				}
				/*else w(i)=0;*/
			}			
		}
		else
		{   
			row=(int)ceil(double(mLongCodBlockLen+4)/double(column));
			kw=3*row*column;
			w.set_size(kw,false);


			ii=0;iii=0;
			for(i=0;i<kw;i++)
			{
				if(i==SelectBitNum(Counter,ii))
				{
					mReceiverMemory(Counter,i)=mTempDataSpace1(Counter,iii);
					/*w(i)=mSpace(counter,iii);*/
					ii+=1;iii+=1;
				}
				/*else{w(i)=0;}*/
			}			
		}
	}
	else
	{
		if(Counter<mShortCodBlockNum)
		{
			row=(int)ceil(double(mShortCodBlockLen)/double(column));
			kw=3*row*column;
			w.set_size(kw,false);
			//nd=row*column-mLongCodBlockLen;
			//k=row*(2*(int)ceil((double)ncb/(double)(8*row))*mRedunVersionNum+2);
			ii=0;iii=0;
			for(i=0;i<kw;i++)
			{
				if(i==SelectBitNum(Counter,ii))
				{
					mReceiverMemory(Counter,i)=mTempDataSpace1(Counter,iii);
					/*w(i)=mSpace(counter,iii);*/
					ii+=1;iii+=1;
				}
				/*else{w(i)=0;}*/
			}			
		}
		else
		{
			row=(int)ceil(double(mLongCodBlockLen)/double(column));
			kw=3*row*column;
			w.set_size(kw,false);
			ii=0;iii=0;
			for(i=0;i<kw;i++)
			{
				if(i==SelectBitNum(Counter,ii))
				{
					mReceiverMemory(Counter,i)=mTempDataSpace1(Counter,iii);
					/*w(i)=mSpace(counter,iii);*/
					ii+=1;iii+=1;
				}
				/*else{w(i)=0;}*/
			}			
		}
	}
	
}

//--- Definition of SeparateBit()
void CTranChannel::SeparateBit(int Counter)    //Inverse of Collecting bits
{
	vec v0,v1,v2;
	int kw,stream_length;
	int i;

    
	if(Counter<mShortCodBlockNum)
		    kw = mShortCollectBitLen;		    
		else
			kw = mLongCollectBitLen;
	//kw = mReceiverMemory.cols();
	stream_length=kw/3;
	v0.set_size(stream_length,false);v1.set_size(stream_length,false);v2.set_size(stream_length,false);

	if(mTrChEncodeType==0)
	{
		for(i=0;i<stream_length;i++)	
		{
			v0(i)=mReceiverMemory(Counter,i);
			v1(i)=mReceiverMemory(Counter,stream_length+2*i);
			v2(i)=mReceiverMemory(Counter,stream_length+2*i+1);
		}
	}
	else
	{
		for(i=0;i<stream_length;i++)
		{
			v0(i)=mReceiverMemory(Counter,i);
			v1(i)=mReceiverMemory(Counter,stream_length+i);
			v2(i)=mReceiverMemory(Counter,stream_length*2+i);
		}
	}
	
	mTempDataSpace1.set_size(3,stream_length,false);
	mTempDataSpace1.clear();
	for(i=0;i<stream_length;i++)
	{
		mTempDataSpace1(0,i)=v0(i);
		mTempDataSpace1(1,i)=v1(i);
		mTempDataSpace1(2,i)=v2(i);
	}
}

//--- Definition of DeSubBlockInterl()
void CTranChannel::DeSubBlockInterl(int Counter)   // Inverse of Sub-block Interleaver
{
    vec d0,d1,d2;
	vec v0,v1,v2;
	mat y0,y1,y2;
	mat y0_con,y1_con,y2_con;
	int stream_length,row,column;
	int nd,d;
	int i,ii;
    
	column=32;
	stream_length=mTempDataSpace1.cols();
	row=stream_length/column;
	v0.set_size(stream_length,false);v1.set_size(stream_length,false);v2.set_size(stream_length,false);
	y0.set_size(row,column,false);y1.set_size(row,column,false);y2.set_size(row,column,false);
	y0_con.set_size(row,column,false);y1_con.set_size(row,column,false);y2_con.set_size(row,column,false);

	for(i=0;i<stream_length;i++)
	{v0(i)=mTempDataSpace1(0,i);v1(i)=mTempDataSpace1(1,i);v2(i)=mTempDataSpace1(2,i);}

	if(mTrChEncodeType==0)
	{
		/************ Output of the DeInterleaver for v0, v1 ************/
		for(i=0;i<column;i++)
		{y0.set_col(i,v0(i*row,row-1+i*row));y1.set_col(i,v1(i*row,row-1+i*row));}

		/*********  Outpur of DeInterleaver for v2 **********/
	    for(i=0;i<column-1;i++)
			y2.set_col(i+1,v2(i*row,row+i*row-1));
	    for(i=0;i<row;i++)
			y2(mod(i+1,row),0)=v2(i+row*(column-1));

		/**************/
		y0.swap_cols(1,16);y0.swap_cols(2,8);y0.swap_cols(3,24);y0.swap_cols(5,20);y0.swap_cols(6,12);y0.swap_cols(7,28);
	    y0.swap_cols(9,18);y0.swap_cols(11,26);y0.swap_cols(13,22);y0.swap_cols(15,30);y0.swap_cols(19,25);y0.swap_cols(23,29);

	    y1.swap_cols(1,16);y1.swap_cols(2,8);y1.swap_cols(3,24);y1.swap_cols(5,20);y1.swap_cols(6,12);y1.swap_cols(7,28);
	    y1.swap_cols(9,18);y1.swap_cols(11,26);y1.swap_cols(13,22);y1.swap_cols(15,30);y1.swap_cols(19,25);y1.swap_cols(23,29);
	    
		/*****************/
		if(Counter<mShortCodBlockNum)
		{
			d=mShortCodBlockLen+4;
			nd=row*column-d;
			d0.set_size(d,false);d1.set_size(d,false);d2.set_size(d,false);
			for(i=0;i<row;i++)
				for(ii=0;ii<column;ii++)
				{
					if(i*column+ii-nd>=0)
					{d0(i*column+ii-nd)=y0(i,ii);d1(i*column+ii-nd)=y1(i,ii);d2(i*column+ii-nd)=y2(i,ii);}
				}
		}
		else
		{
			d=mLongCodBlockLen+4;
			nd=row*column-d;
			d0.set_size(d,false);d1.set_size(d,false);d2.set_size(d,false);
			for(i=0;i<row;i++)
				for(ii=0;ii<column;ii++)
				{
					if(i*column+ii-nd>=0)
					{d0(i*column+ii-nd)=y0(i,ii);d1(i*column+ii-nd)=y1(i,ii);d2(i*column+ii-nd)=y2(i,ii);}
				}
		}
	}
	else
	{
		for(i=0;i<column;i++)
		{y0.set_col(i,v0(i*row,row-1+i*row));y1.set_col(i,v1(i*row,row-1+i*row));y2.set_col(i,v2(i*row,row-1+i*row));}
		
		/********** Perform the inter-column permutation patter for d0, d1 and d2 ***********/
		y0_con.set_col(0,y0.get_col(1));y0_con.set_col(1,y0.get_col(17));y0_con.set_col(2,y0.get_col(9));y0_con.set_col(3,y0.get_col(25));
		y0_con.set_col(4,y0.get_col(5));y0_con.set_col(5,y0.get_col(21));y0_con.set_col(6,y0.get_col(13));y0_con.set_col(7,y0.get_col(29));
		y0_con.set_col(8,y0.get_col(3));y0_con.set_col(9,y0.get_col(19));y0_con.set_col(10,y0.get_col(11));y0_con.set_col(11,y0.get_col(27));
		y0_con.set_col(12,y0.get_col(1));y0_con.set_col(13,y0.get_col(23));y0_con.set_col(14,y0.get_col(15));y0_con.set_col(15,y0.get_col(31));
		y0_con.set_col(16,y0.get_col(0));y0_con.set_col(17,y0.get_col(16));y0_con.set_col(18,y0.get_col(8));y0_con.set_col(19,y0.get_col(24));
		y0_con.set_col(20,y0.get_col(4));y0_con.set_col(21,y0.get_col(20));y0_con.set_col(22,y0.get_col(12));y0_con.set_col(23,y0.get_col(28));
		y0_con.set_col(24,y0.get_col(2));y0_con.set_col(25,y0.get_col(18));y0_con.set_col(26,y0.get_col(10));y0_con.set_col(27,y0.get_col(26));
        y0_con.set_col(28,y0.get_col(6));y0_con.set_col(29,y0.get_col(22));y0_con.set_col(30,y0.get_col(14));y0_con.set_col(31,y0.get_col(30));

		y1_con.set_col(0,y1.get_col(1));y1_con.set_col(1,y1.get_col(17));y1_con.set_col(2,y1.get_col(9));y1_con.set_col(3,y1.get_col(25));
		y1_con.set_col(4,y1.get_col(5));y1_con.set_col(5,y1.get_col(21));y1_con.set_col(6,y1.get_col(13));y1_con.set_col(7,y1.get_col(29));
		y1_con.set_col(8,y1.get_col(3));y1_con.set_col(9,y1.get_col(19));y1_con.set_col(10,y1.get_col(11));y1_con.set_col(11,y1.get_col(27));
		y1_con.set_col(12,y1.get_col(1));y1_con.set_col(13,y1.get_col(23));y1_con.set_col(14,y1.get_col(15));y1_con.set_col(15,y1.get_col(31));
		y1_con.set_col(16,y1.get_col(0));y1_con.set_col(17,y1.get_col(16));y1_con.set_col(18,y1.get_col(8));y1_con.set_col(19,y1.get_col(24));
		y1_con.set_col(20,y1.get_col(4));y1_con.set_col(21,y1.get_col(20));y1_con.set_col(22,y1.get_col(12));y1_con.set_col(23,y1.get_col(28));
		y1_con.set_col(24,y1.get_col(2));y1_con.set_col(25,y1.get_col(18));y1_con.set_col(26,y1.get_col(10));y1_con.set_col(27,y1.get_col(26));
        y1_con.set_col(28,y1.get_col(6));y1_con.set_col(29,y1.get_col(22));y1_con.set_col(30,y1.get_col(14));y1_con.set_col(31,y1.get_col(30));

		y2_con.set_col(0,y2.get_col(1));y2_con.set_col(1,y2.get_col(17));y2_con.set_col(2,y2.get_col(9));y2_con.set_col(3,y2.get_col(25));
		y2_con.set_col(4,y2.get_col(5));y2_con.set_col(5,y2.get_col(21));y2_con.set_col(6,y2.get_col(13));y2_con.set_col(7,y2.get_col(29));
		y2_con.set_col(8,y2.get_col(3));y2_con.set_col(9,y2.get_col(19));y2_con.set_col(10,y2.get_col(11));y2_con.set_col(11,y2.get_col(27));
		y2_con.set_col(12,y2.get_col(1));y2_con.set_col(13,y2.get_col(23));y2_con.set_col(14,y2.get_col(15));y2_con.set_col(15,y2.get_col(31));
		y2_con.set_col(16,y2.get_col(0));y2_con.set_col(17,y2.get_col(16));y2_con.set_col(18,y2.get_col(8));y2_con.set_col(19,y2.get_col(24));
		y2_con.set_col(20,y2.get_col(4));y2_con.set_col(21,y2.get_col(20));y2_con.set_col(22,y2.get_col(12));y2_con.set_col(23,y2.get_col(28));
		y2_con.set_col(24,y2.get_col(2));y2_con.set_col(25,y2.get_col(18));y2_con.set_col(26,y2.get_col(10));y2_con.set_col(27,y2.get_col(26));
        y2_con.set_col(28,y2.get_col(6));y2_con.set_col(29,y2.get_col(22));y2_con.set_col(30,y2.get_col(14));y2_con.set_col(31,y2.get_col(30));

		/***********  The output of the block interleaver ************/
		if(Counter<mShortCodBlockNum)
		{
			d=mShortCodBlockLen;
			nd=row*column-d;
			for(i=0;i<row;i++)
				for(ii=0;ii<column;ii++)
				{
					if(i*column+ii-nd>=0)
					{d0(i*column+ii-nd)=y0(i,ii);d1(i*column+ii-nd)=y1(i,ii);d2(i*column+ii-nd)=y2(i,ii);}
				}
		}
		else
		{
			d=mLongCodBlockLen;
			nd=row*column-d;
			for(i=0;i<row;i++)
				for(ii=0;ii<column;ii++)
				{
					if(i*column+ii-nd>=0)
					{d0(i*column+ii-nd)=y0(i,ii);d1(i*column+ii-nd)=y1(i,ii);d2(i*column+ii-nd)=y2(i,ii);}
				}
		}
	}

    mTempDataSpace1.set_size(3,d0.length(),false);
	mTempDataSpace1.clear();
	for(i=0;i<mTempDataSpace1.cols();i++)
	{
		mTempDataSpace1(0,i)=d0(i);
		mTempDataSpace1(1,i)=d1(i);
		mTempDataSpace1(2,i)=d2(i);
	}

}

//--- Definition of DecodTrCh()
void CTranChannel::DecodTrCh(int Counter)             //Decoding transmitting channel
{
	if(mTrChEncodeType==0)
	 DecodOfTurbo(Counter);
  else
	  if(mTrChEncodeType==1)
	    DecodOfTailBiting(Counter);
}

//--- Definition of DecodOfTailBiting()
void CTranChannel::DecodOfTailBiting(int Counter)        //Tail Biting Convolutional decoding
{
	int i,j;
    Convolutional_Code code;
    ivec generator(3);
    generator(0)=0133;
    generator(1)=0165;
    generator(2)=0171;
    code.set_generator_polynomials(generator, 7);
    bvec decoded_output_bits;
    vec decoded_input_bits,d0,d1,d2;
  
   d0.set_size(mTempDataSpace1.cols(),false);
   d1.set_size(mTempDataSpace1.cols(),false);
   d2.set_size(mTempDataSpace1.cols(),false);

    for(i=0;i<mTempDataSpace1.cols();i++)
	{
		d0(i)=mTempDataSpace1(0,i);
        d1(i)=mTempDataSpace1(1,i);
        d2(i)=mTempDataSpace1(2,i);
    }

    decoded_input_bits.set_size(3*mTempDataSpace1.cols(),false);

    j=0;
    for(i=0;i<mTempDataSpace1.cols();i++)
	{
		decoded_input_bits(j)=d0(i);j++;
	    decoded_input_bits(j)=d1(i);j++;
	    decoded_input_bits(j)=d2(i);j++;
	}

   code.decode_tail(decoded_input_bits,decoded_output_bits);

}

//--- Definition of DecodOfTurbo()
void CTranChannel::DecodOfTurbo(int Counter)           //Turbo decoding
{
	//Local variables: 
	int i, j;
	int constraint_length, CodedLen;
    bvec decoded_output_bits;
    vec decoded_input_bits, d0, d1, d2;
    ivec gen(2), interleaver_sequence;

    //Declarations of the turbo class:
    Turbo_Codec turbo;

    // Determine the number of code blocks and code block length:
    if(mCodBlockNum == 1) 
	   CodedLen = mLongCodBlockLen;
    else      
	    if(Counter < mShortCodBlockNum)
          CodedLen = mShortCodBlockLen;
	    else
          CodedLen = mLongCodBlockLen;  

    // Initializations:
    gen(0) = 015; 
    gen(1) = 013;
	constraint_length = 4;
    d0.set_size(mTempDataSpace1.cols(),false);
    d1.set_size(mTempDataSpace1.cols(),false);
    d2.set_size(mTempDataSpace1.cols(),false);
	decoded_input_bits.set_size(3*mTempDataSpace1.cols(),false);

    // Set the Turbo encoder/decoder parameters:
    turbo.set_parameters(gen, gen, constraint_length, mTurboIntraInterlSeq, 8, "LOGMAX", 0.75);

    // Get the coded bits from the public space named mTempDataSpace1:
    for(i = 0;i < mTempDataSpace1.cols(); i++)
	{
		d0(i) = mTempDataSpace1(0,i);
        d1(i) = mTempDataSpace1(1,i);
        d2(i) = mTempDataSpace1(2,i);
    }

	// Set the bit counter to zero:
    j = 0;


    for(i = 0; i < mTempDataSpace1.cols(); i++)
	{
		decoded_input_bits(j) = d0(i); j++;
	    decoded_input_bits(j) = d1(i); j++;
	    decoded_input_bits(j) = d2(i); j++;
	}

   // Decode the code block:
   turbo.decode(decoded_input_bits,decoded_output_bits);
   
   // Store the encoded_bits in the public space named mTempDataSpace:
   mTempDataSpace.set_size(1,decoded_output_bits.length(),false);
   mTempDataSpace.clear();
   for(i = 0; i < decoded_output_bits.length(); i++)
      mTempDataSpace(0,i) = decoded_output_bits(i);
   
}

//--- Definition of CrcCheck()
void CTranChannel::CrcCheck()               //CRC check
{
   	int i;
	ivec get_bits;
	bvec crc_input_bits,crc_output_bits;
	CRC_Code crc(std::string("CRC-24"));

	get_bits.set_size(mTempDataSpace.cols(),false);
	crc_input_bits.set_size(mTempDataSpace.cols(),false);
    
    for(i=0;i<mTempDataSpace.cols();i++)
        get_bits(i)=mTempDataSpace(0,i);
	crc_input_bits=to_bvec(get_bits);

	// CRC check
	mCorrectness = crc.check_parity(crc_input_bits);
	//cout << "CRC check in function=" << mCorrectness << endl;
}

//--- Definition of DeSegAttachCrc()
void CTranChannel::DeSegAttachCrc(int Counter)         // Inverse of CRC attchment
{   
	ivec get_bits;
    bvec crc_input_bits,crc_output_bits;
	int i,j,length;
	CRC_Code crc(std::string("CRC-24"));

	get_bits.set_size(mTempDataSpace.cols(),false);
	crc_input_bits.set_size(mTempDataSpace.cols(),false);
    

    for(i=0;i<mTempDataSpace.cols();i++)
        get_bits(i)=mTempDataSpace(0,i);
	crc_input_bits=to_bvec(get_bits);

	crc.decode(crc_input_bits,crc_output_bits);
	
	for(i=0;i<crc_output_bits.length();i++)
	{
		mSpace1(Counter,i)=crc_output_bits(i);
		mpTrChOutputData(i)=crc_output_bits(i);
	}
      
	if(mCodBlockNum==1)
	{
		length=crc_output_bits.length()-mCrcFillBitNum;
		mpTrChOutputData.set_size(length,false);
	    for(i=mCrcFillBitNum,j=0;i<crc_output_bits.length();i++,j++)
          mpTrChOutputData(j)=crc_output_bits(i);
	}
}

//--- Definition of DeSegCodBlock()
void CTranChannel::DeSegCodBlock()      // Inverse of Segment coding blocks
{
    int i,j,k,length,l;
	ivec input;
	k=0;

    l = mShortCodBlockNum * mShortCodBlockLen + mLongCodBlockNum * mLongCodBlockLen - mCrcLen * mCodBlockNum;
	input.set_size(l,false);
	for(i=0;i<mCodBlockNum;i++)
	{
		
	  if(mCodBlockNum==1) 
	     length=mLongCodBlockLen-mCrcLen;
      else      
	    if(i<mShortCodBlockNum)
           length=mShortCodBlockLen-mCrcLen;
	    else
          length=mLongCodBlockLen-mCrcLen;  
      
      
	  for(j=0;j<length;j++)	
	  {input(k)=mSpace1(i,j);k++;}
	  
	}
   mTempDataSpace.set_size(1,input.length(),false);
   mTempDataSpace.clear();
   for(i=0;i<input.length();i++)
       mTempDataSpace(0,i) = input(i);
}

//--- Definition of DeAttachCrc()
void CTranChannel::DeAttachCrc()     // Inverse of CRC check of Segment coding blocks
{   
	int i,j,length;
	bvec crc_input_bits,crc_output_bits;
    CRC_Code crc(std::string("CRC-24"));

	crc_input_bits.set_size(mTempDataSpace.cols(),false);
	for(i=0;i<mTempDataSpace.cols();i++)
       crc_input_bits(i) = mTempDataSpace(0,i);

	crc.decode(crc_input_bits,crc_output_bits);

    length=crc_output_bits.length()-mCrcFillBitNum;
	mpTrChOutputData.set_size(length,false);
      for(i=mCrcFillBitNum,j=0;i<crc_output_bits.length();i++,j++)
        mpTrChOutputData(j)=crc_output_bits(i);
}


//The functions class CULSch
 

CULSCh::CULSCh()
{
	mCodBlockLoopCounter=0;
}


CULSCh::~CULSCh()
{

}


void CULSCh::Run()
{
  imat TempMatchRateOut;
  int TempMatrixCol;
  mTrChEncodeType=0;

  //GetData();
  //AttachCrc();
  //SegCodBlock();
  /*AttachCrc();*/
  TempMatrixCol=(int)ceil((double)mTranBlockSize/(double)mCodBlockNum);
  TempMatchRateOut.set_size(mCodBlockNum,TempMatrixCol,false);


  for(mCodBlockLoopCounter=0; mCodBlockLoopCounter<mCodBlockNum; mCodBlockLoopCounter++)
  {
	  cout<<"*******************mCodBlockLoopCounter**************************"<<mCodBlockLoopCounter<<endl;
	  EncodTrCh(mCodBlockLoopCounter);
	  //MatchRate(mCodBlockLoopCounter);      
	  //for(i=0;i<mTempDataSpace.cols();i++)
	  //TempMatchRateOut(mCodBlockLoopCounter,i)=mTempDataSpace(0,i);
   //   cout<<TempMatchRateOut<<"********"<<endl;
  }
  
  //ConcatCodeBlock(TempMatchRateOut);
  
}


void CULSCh::DRun()
{
  
 // /*GetData();*/
 // SeparCodeBlock();
 // for( mCodBlockLoopCounter=0; mCodBlockLoopCounter<mCodBlockNum; mCodBlockLoopCounter)
 // {
	//DeMatchRate();
	//DecodTrCh();
	//DeSegAttachCrc();

 // }
 // DeSegCodBlock();
 // DeAttachCrc();
 // mTempDataSpace;
}


void CULSCh::SetCULSChDataForPhysicalCh()
{
		
}


void CULSCh::InniCTranChannel(int SoftbuffSize)
{
	//mRedunVersionNum=ReduVersiNum;
	mSoftbuffSize=SoftbuffSize;
	//cout<<"CULSCh:mReduVersiNum="<<mRedunVersionNum<<endl;

}


//The functions class CDLSch
 

CDLSCh::CDLSCh()
{

}
CDLSCh::~CDLSCh()
{

}

void CDLSCh::InniCTranChannel(int SoftbuffSize)
{
	mSoftbuffSize=SoftbuffSize;
}


void CDLSCh::Run()
{
  //Local variables:
  int i,g,r;
  int TempMatrixCol;
  imat TempMatchRateOut;

  // Initializations:
  mTrChEncodeType = 0;
  if(mRedunVersionNum == 0)   // RV=0 denotes New Data should be transmitted
  {
    AttachCrc(); // Code block segmentation and code block CRC attachment:
	// Calculate the size of transmission block
    g = mTranBlockSize / (mNl * mQm);
    r = mod(g, mCodBlockNum);
    TempMatrixCol = mNl * mQm * (int) ceil ((double)g / (double)mCodBlockNum);
    TempMatchRateOut.set_size(mCodBlockNum, TempMatrixCol, false);
	SelectBitNum.set_size(mCodBlockNum, TempMatrixCol+1, false);
	mBufferSpace.set_size(mCodBlockNum, 3*mLongCodBlockLen+132, false);

	// Start the loop for all the code blocks:
    for(mCodBlockLoopCounter = 0; mCodBlockLoopCounter < mCodBlockNum; mCodBlockLoopCounter++)
	{
		// Encode a single block: 
		EncodTrCh(mCodBlockLoopCounter);
		
		// Rate matching:
		MatchRate(mCodBlockLoopCounter);
      
		for(i=0;i<mTempDataSpace.cols();i++)
		{ 
			TempMatchRateOut(mCodBlockLoopCounter, i) = mTempDataSpace(0, i);
		}
	}  
	// Code block concatenation:
	ConcatCodeBlock(TempMatchRateOut);
  }
  
  else  // Retransmission:
  {
	 // Calculate the size of transmission block
     g=mTranBlockSize/(mNl*mQm);
     r=mod(g,mCodBlockNum);
     TempMatrixCol=mNl*mQm*(int)ceil((double)g/(double)mCodBlockNum);
     TempMatchRateOut.set_size(mCodBlockNum,TempMatrixCol,false);
	 TempMatchRateOut.clear();
	 // Start the loop for all the code blocks:
	 for(mCodBlockLoopCounter = 0; mCodBlockLoopCounter < mCodBlockNum; mCodBlockLoopCounter++)
	 {
	     // Bit selection and transmission
		 SelectBit(mCodBlockLoopCounter);
		 for(i = 0; i < mTempDataSpace.cols(); i++)
		 {
			 TempMatchRateOut(mCodBlockLoopCounter, i) = mTempDataSpace(0, i);
		 }
	 }	 
	 // Code block concatenation:
	 ConcatCodeBlock(TempMatchRateOut);
  }
}


void CDLSCh::DRun()
{  
  
	mReceiverMemory.set_size(mCodBlockNum,mLongCollectBitLen,false);
	mSpace1.set_size(mCodBlockNum,mLongCodBlockLen-mCrcLen,false);
    
	// Clear the cache mReceiverMemory:
	if(mRedunVersionNum == 0)
	{
	    mReceiverMemory.zeros();
	}
	
	// Code block segmentation：
    SeparCodeBlock();

	// Start the loop for all the code blocks:
	for( mCodBlockLoopCounter=0; mCodBlockLoopCounter < mCodBlockNum; mCodBlockLoopCounter++)
	{
        // Rate de-matching:
		DeMatchRate(mCodBlockLoopCounter);
		
		// Decoding:
		DecodTrCh(mCodBlockLoopCounter); 	

        // CRC Check in the Sub-code block:
		CrcCheck();

        // Remove CRC sequence of Sub-code block:
		DeSegAttachCrc(mCodBlockLoopCounter); 
    }
	
	if(mCodBlockNum != 1)
	{ 
	    // Concatenate all the sub-code blocks into a transmission block:
		DeSegCodBlock();

	    // CRC Check in the transport block:
		CrcCheck();
	   
		// Remove CRC sequence of the transport block:
		DeAttachCrc();
	}
}
void CDLSCh::SetCDLSChDataForPhysicalCh()
{
		
}


//Class of UCI  传输信道上行控制信息类设计
CUCI:: CUCI ()
{}
CUCI::~ CUCI ()
{}
void CUCI::Run()
{ PNGen();
mConInforEncodeType=2;
if(mConInforEncodeType==1)
{EncodforACK();
}
if(mConInforEncodeType==2)
{ EncodforCQI();
}
if(mConInforEncodeType==3)
{ EncodforRI();
}
return;
}


void CUCI::DRun()
{ mConInforEncodeType=2;
if(mConInforEncodeType==1)
{DecodforACK();
}
if(mConInforEncodeType==2)
{ DecodforCQI();
}
if(mConInforEncodeType==3)
{ DecodforRI();
}
return;
}

 void CUCI::cqi_code_list(string filename,imat &a)
{
    
    ifstream fin(filename.c_str());
    if(!fin) {cout<<"cqi_code.txt"<<endl;return;}

int b[11],i=0,n=0;
 
    for(int i=0;i<a.rows();i++)

 
      {
		  for(n=0;n<11;n++) 
		  {fin>>b[n];
		  a(i,n)=b[n];}
      }

    fin.close();

}
void CUCI:: EncodforACK()
{
	int Qm=4,Q_ACK=20,i;   // Assumption
	bvec q_ACK=randb(Q_ACK);
	cout<<"q_ACK"<<q_ACK<<endl;
	Array<bvec> qp_ACK(Q_ACK/Qm);
	int k=0;
	for(i=0;i<Q_ACK;i=i+Qm)
	{
		qp_ACK(k)=q_ACK(i,i+Qm-1).transpose();
		k=k+1;
		cout<<"qp_ACK="<<qp_ACK<<endl;
	}
}

void CUCI:: EncodforCQI()
{
	int B=32,N_CQI=12;// CQ bits depends on the transmission format.less than 12
	int Q_CQI=10;
	bvec o=randb(N_CQI);
	bvec bb(B),qq(Q_CQI);
	Convolutional_Code code;
    ivec generator(3);
    generator(0)=0133;
    generator(1)=0165;
    generator(2)=0171;
    code.set_generator_polynomials(generator,7);

    bvec transmitted_bits, encoded_bits;
    transmitted_bits=o;
   
    code.encode_tailbite(transmitted_bits, encoded_bits);
    cout<<"transmitted_bits="<<transmitted_bits<<endl;
    cout<<"encoded_bits="<<encoded_bits<<endl;
    mpTrChOutputData=encoded_bits;
//
// imat M(32,11);
// cqi_code_list("cqi_code.txt",M);
// cout<<"aa="<<M<<endl;
//
//for(i=0;i<B;i++)
//{sum=0;
//for(n=0;n<N_CQI;n++)
//{kk=o(n);
//sum=sum+M(i,n)*kk;}
//cout<<"kk="<<sum<<endl;
// bb(i)=mod(sum,2);
// cout<<"bb(i)="<<bb(i)<<endl;
/*for(i=0;i<Q_CQI;i++)
{qq(i)=bb(mod(i,B));
cout<<"qq(i)="<<qq(i)<<endl;}*/
}
void CUCI:: EncodforRI()
{
	int Qm=4,Q_RI=20,i;   //Assumption
	bvec q_RI=randb(Q_RI);
	cout<<"q_RI"<<q_RI<<endl;
	Array<bvec> qp_RI(Q_RI/Qm);
	int k=0;
	for(i=0;i<Q_RI;i=i+Qm)
	{
		qp_RI(k)=q_RI(i,i+Qm-1).transpose();
		k=k+1;
		cout<<"qp_RI="<<qp_RI<<endl;
	}
}

void CUCI:: DecodforACK()
{
}

void CUCI:: DecodforCQI()
{ 
  Convolutional_Code code;
  ivec generator(3);
  generator(0)=0133;
  generator(1)=0165;
  generator(2)=0171;
  code.set_generator_polynomials(generator, 7);
  bvec decoded_input_bits,decoded_output_bits;
  vec decoded_input_bits1;
  BPSK bpsk;
  decoded_input_bits=mpTrChOutputData;
  decoded_input_bits1=bpsk.modulate_bits(decoded_input_bits);
  code.decode_tailbite(decoded_input_bits1,decoded_output_bits);
  cout<<"decoded_output_bits="<<decoded_output_bits<<endl;
  
}

void CUCI:: DecodforRI()
{
}

void CUCI:: PNGen()
{ rand();
}

///Class of DCI  传输信道下行控制信息类设计
CDCI::CDCI()
{}
CDCI:: ~ CDCI ()
{}
void CDCI:: Run()
{ /*PNGen();
AttachCrc();
EncodTrCh();
MatchRate();*/


	int A=10,L=16;///CRC-16
	vec coded_bits_vec;
	bvec a_bits = randb(A), coded_bits, segbits,coded_segbits,decoded_bits,b_bits,c_bits(A+L);
	bvec xrnti=zeros_b(16);//// need to be modified
    cout<<"xrnti="<<xrnti<<endl;
   
    CRC_Code crc(std::string("SDLC-16"));
    b_bits = crc.encode(a_bits);
	 
    Convolutional_Code code;
    ivec generator(3);
    generator(0)=0133;
    generator(1)=0165;
    generator(2)=0171;
    code.set_generator_polynomials(generator,7);

    bvec  encoded_bits;
    code.encode_tailbite(b_bits, encoded_bits);
    cout<<"a_bits="<<a_bits<<endl;
    cout<<"b_bits="<<b_bits<<endl;
    mpTrChOutputData=encoded_bits;
}


void CDCI:: DRun()
{
	bvec decoded_bits;
    //DeMatchRate();
    Convolutional_Code code;
    ivec generator(3);
    generator(0)=0133;
    generator(1)=0165;
    generator(2)=0171;
    code.set_generator_polynomials(generator, 7);
    bvec decoded_input_bits,decoded_output_bits;
    vec decoded_input_bits1;
    BPSK bpsk;
    decoded_input_bits=mpTrChOutputData;
    decoded_input_bits1=bpsk.modulate_bits(decoded_input_bits);
    code.decode_tailbite(decoded_input_bits1,decoded_output_bits);
    cout<<"decoded_output_bits="<<decoded_output_bits<<endl;
    CRC_Code crc(std::string("SDLC-16"));
    crc.decode(decoded_output_bits,decoded_bits );
	cout<<"decoded_bits="<<decoded_bits<<endl;
	return;
}