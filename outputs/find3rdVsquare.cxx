//Make noise waveforms from data
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <memory>
#include <boost/scoped_ptr.hpp>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fftw3.h>
#include <ctime>
#include "time.h" // for time convert

//ROOT Includes
#include "TTreeIndex.h"
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TText.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TPaveStats.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TMath.h"
#include "TF1.h"
#include "TH1D.h"
#include "TChain.h"
#include "TString.h"
#include "TVector3.h"
#include "TLegend.h"
#include "TStyle.h"


//AraSim Includes
#ifndef ___CINT___ //shielded from CINT
#include "Settings.h"
#endif

#include "RawIcrrStationEvent.h"  
#include "UsefulIcrrStationEvent.h"
#include "RawAraStationEvent.h"
#include "RawAtriStationEvent.h"
#include "UsefulAtriStationEvent.h"
#include "AraEventCalibrator.h"
#include "AraGeomTool.h"
#include "FFTtools.h"
#include "AraStationInfo.h"
#include "TPaletteAxis.h"
TStyle* RootStyle();

using namespace std;

//FFT function

int main(int argc, char **argv)
{ //start main

  TChain chain("eventTree"); //this for the events for the exterior loop     
  for(int file=1; file<argc; file++){
    TString fileKey(argv[file]); //a file key                                
    chain.Add(fileKey); //add files to the chain                             
  }
  int channel = atof(argv[1]);
  RawAtriStationEvent *rawEvPtr=0;
  chain.SetBranchAddress("event",&rawEvPtr);
  int numEntries = chain.GetEntries();
  Int_t run_num=0;
  int stationId=0;

  stationId = rawEvPtr->stationId; //assign the statio id number
  
  chain.SetBranchAddress("run",&run_num); //set the branch address

  AraEventCalibrator *calib = AraEventCalibrator::Instance(); //make a calibrator
  TH1* h1 = new TH1I("h1", "h1 title", 200, 0.0, 15000);

  for(int event=0; event<numEntries; event++){//loop over events
    chain.GetEvent(event);
    int evt_num = rawEvPtr->eventNumber;//event number
    UsefulAtriStationEvent *realAtriEvPtr_fullcalib = new UsefulAtriStationEvent(rawEvPtr, AraCalType::kLatestCalib); //make the event
    bool is_RF_trigger = rawEvPtr->isRFTrigger();
    if(is_RF_trigger){//If RF triggered event
      //  cout<<"RF triggered!"<<endl;
      int eventId=rawEvPtr->eventId;
      int numReadoutBlocks=rawEvPtr->numReadoutBlocks;
      uint32_t triggerInfo[MAX_TRIG_BLOCKS];
      uint8_t triggerBlock[MAX_TRIG_BLOCKS];
      for(int trig=0;trig<MAX_TRIG_BLOCKS;trig++) {
	// triggerInfo[trig]=rawEvPtr->triggerInfo[trig];
	triggerBlock[trig]=rawEvPtr->triggerBlock[trig];
	//cout << "Trigger info: " << triggerInfo[trig] << endl;
	//cout << "Trigger block: " << triggerBlock[trig] << endl;
	TGraph *waveform = realAtriEvPtr_fullcalib->getGraphFromRFChan(channel);//channel 2.
	double vsquared = FFTtools::getPeakSqVal(waveform);
	h1->Fill(vsquared);
	
      }
    }//close if RF triggered event
    
  }//end loop over events
  TCanvas *c2 = new TCanvas("","",650,650);
  h1->Draw();
  c2->SaveAs("vsquared.png");
}//end main
