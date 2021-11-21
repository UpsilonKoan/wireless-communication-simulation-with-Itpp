//----- ReferenceSignal.cpp -----
// Programmed by LTE Simulation Group of IMC, SWJTU
#include "stdafx.h"
#include "itpp/itcomm.h"
#include "ReferenceSignal.h"
using std::cin;
using std::cout;
using std::endl;
using namespace itpp;

// *************************************** Class CTranceiver Definition *****************************************

CReferenceSignal::CReferenceSignal()
{}

CReferenceSignal:: ~ CReferenceSignal()
{}

int CReferenceSignal::GenURSSequence()
{return 0;}
int CReferenceSignal::GenDRSSequence()
{return 0;}
int CReferenceSignal::DemodRSmap()
{return 0;}
int CReferenceSignal::SoundingRSmap()
{return 0;}
int CReferenceSignal::CommonRSmap()
{return 0;}
int CReferenceSignal::DedicatedRSmap()
{return 0;}
int CReferenceSignal::GetDemodRS()
{return 0;}
int CReferenceSignal::GetSoundingRS()
{return 0;}
int CReferenceSignal::GetCommonRS()
{return 0;}
int CReferenceSignal::GetDedicatedRS()
{return 0;}
int CReferenceSignal::Run() 
{
	/*if (mRSState[0])
{
	GenURSSequence();
    SoundingRSmap();
}
if(mRSState[1])
{
	GenURSSequence();
    DemodRSmap();
}
if(mRSState[2])
{ 
	GenDRSSequence();
    CommonRSmap(); 
}
if(mRSState[3])
{
    GenDRSSequence();
    DedicatedRSmap();
}*/
return 0;
}
int  CReferenceSignal::DRun()
{
//if (mRSState[0])
//GetSoundingRS();
//if (mRSState[1])
//GetDemodRS();
//if( mRSState[2])
//GetCommonRS();
//if(mRSState[3])
//GetDedicatedRS();
return 0;}
