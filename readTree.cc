#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <time.h>
#include "TTreeIndex.h"
#include "TChain.h"
#include "TH1.h"
#include "TF1.h"
#include "TF2.h"
#include "TFile.h"
#include "TRandom.h"
#include "TRandom2.h"
#include "TRandom3.h" 
#include "TTree.h"
#include "TLegend.h"
#include "TLine.h"
#include "TROOT.h"
#include "TPostScript.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TText.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TMath.h"
#include <unistd.h>
#include "TVector3.h"
#include "TRotation.h"
#include "TSpline.h"
//$(FFTWSYS) -llibRootFftwWrapper
//#include <fftw3.h>

using namespace std;

#include "RawIcrrStationEvent.h"  
#include "UsefulIcrrStationEvent.h"
#include "RawAraStationEvent.h"
#include "RawAtriStationEvent.h"
#include "UsefulAtriStationEvent.h"
#include "AraEventCalibrator.h"
//#include "FFTtools.h"


#include "PlottingFns.h"
#include "Constants.h"


class EarthModel; //class




int getPeakBin(TGraph *gr);

double getPeak(TGraph *gr);

RawAtriStationEvent *rawAtriEvPtr;
UsefulAtriStationEvent *realAtriEvPtr;

//int main() {
int main(int argc, char **argv) {    // this is for manual power threshold value


  string readfile;
  if (argc<2) { // no setup file input, use default
    readfile = "outputs/AraOut.root";
  }
  else if (argc == 2) { // read file!!
    readfile = string( argv[1] );
  }
  else { // no mode for argc > 2!
    cout<<"too many info! just use default AraOut.root file!"<<endl;
    readfile = "outputs/AraOut.root";
  }


  //  Settings *settings = new Settings();

  //  Detector *detector=new Detector(settings->DETECTOR); // builds antenna array, 0 for testbed
  //  Detector *detector=0; // builds antenna array, 0 for testbed
 

  
  TFile *AraFile=new TFile(( readfile ).c_str());
  //TFile *AraFile=new TFile((outputdir+"/AraOut.root").c_str());
  TTree *eventTree=(TTree*)AraFile->Get("eventTree");
  eventTree->SetBranchAddress("UsefulAtriStationEvent", &realAtriEvPtr);
  eventTree->GetEvent(0);
 

  int numEntries = eventTree->GetEntries();
  stringstream ss;
  string xLabel, yLabel;
  vector<string> titlesForGraphs;
  for (int i = 0; i < 16; i++){
    ss.str("");
    ss << "Channel " << i;
    titlesForGraphs.push_back(ss.str());
  }
  xLabel = "Time (ns)"; yLabel = "Voltage (mV)";
  vector<TGraph*> grWaveformsRaw = makeGraphsFromRF(realAtriEvPtr, 16, xLabel, yLabel, titlesForGraphs);
  
  ss.str("");

      
}


int getPeakBin(TGraph *gr) 
{
  double x,y;
  gr->GetPoint(0,x,y);
  double peakVal=y;
  int peakBin=0;
  for(int i=1;i<gr->GetN();i++) {
    gr->GetPoint(i,x,y);
    if(peakVal<y) {
      peakVal=y;
      peakBin=i;
    }      
  }
  return peakBin;
}

double getPeak(TGraph *gr)
{
  double x,y;
  gr->GetPoint(0,x,y);
  double peakVal=y*y;
  int peakBin=0;
  for(int i=1;i<gr->GetN();i++) {
    gr->GetPoint(i,x,y);
    if( peakVal<(y*y) ) {
      peakVal=(y*y);
      peakBin=i;
    }      
  }
  //return peakBin;
  return peakVal;
}




