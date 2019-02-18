//////////////////////////////////////////////////////////////////////////////                                                 
/////  AraCuts.cxx                                                       /////                                                  
/////                                                                    /////                                                  
/////  Description:                                                      /////                                                  
/////     Class for making a variety of cuts on Ara Data                 /////                                                  
/////  Author: CGP, adapted from ACG, class structure from RJN            /////                                                 
//////////////////////////////////////////////////////////////////////////////  

#ifndef PLOTTINGFNS_H
#define PLOTTINGFNS_H

                                                
//System includes                                                                                                               
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <cstdlib>
#include <sstream>
#include <functional>
#include <numeric>
#include <deque>
using namespace std;

//AraRoot Includes                                                                                                              
#include "RawIcrrStationEvent.h"
#include "RawAtriStationEvent.h"
#include "UsefulAraStationEvent.h"
#include "UsefulIcrrStationEvent.h"
#include "UsefulAtriStationEvent.h"
#include "AraGeomTool.h"
#include "AraAntennaInfo.h"
//#include "RayTraceCorrelator_test.h"
#include "FFTtools.h"
#include "WaveformFns.h"

//ROOT Includes                                                                                                
#include "TROOT.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
#include "TButton.h"
#include "TGroupButton.h"
#include <TGClient.h>
#include "TStyle.h"
#include "TPostScript.h"
#include "TTree.h"
#include "math.h"
#include "TText.h"
#include "TF1.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TNtuple.h"
#include "TProfile.h"
#include "Math/Interpolator.h"
#include "TImage.h"
#include "TMarker.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TObjArray.h"
#include "TLine.h"



void saveGraph(TGraph* gr1, string filename);
void save2Graphs(TGraph* gr1, TGraph* gr2, string filename);
void save16Graphs(vector<TGraph*> grs, string filename);
void saveNGraphs(vector<TGraph*> grs, string filename, int numGraphs);

void saveHist(TH1* hist, string filename);

void deleteGraphVector(vector<TGraph*> graphs);

vector<TGraph*> makeGraphsFromRF(UsefulAtriStationEvent* realAtriEvPtr, int numGraphs, string xlabel, string ylabel, vector<string> titles);
vector<TGraph*> makeInterpolatedGraphs(vector<TGraph*> graphsIn, double intTimestep, string xlabel, string ylabel, vector<string> titles);


void saveGraph(TGraph* gr1, string filename){

  TCanvas* c1 = new TCanvas("c1", "", 640, 480);
  gr1->Draw("AL");
  c1->SaveAs(filename.c_str());

  delete c1;
};

void saveHist(TH1* hist, string filename){

  TCanvas* c1 = new TCanvas("c1", "", 640, 480);
  hist->Draw();
  c1->SaveAs(filename.c_str());

  delete c1;
};

void saveHistWithFit(TH1* hist, TF1* myFit, string filename){

  TCanvas* cHist = new TCanvas("c1", "", 640, 480);
  hist->Draw();
  myFit->Draw("same");
  cHist->SaveAs(filename.c_str());

  delete cHist;
};


void save2Graphs(TGraph* gr1, TGraph* gr2, string filename){
  
  TCanvas* c1 = new TCanvas("c1", "", 640, 480);
  c1->Divide(1,2);
  c1->cd(1);
  gr1->Draw("AL");
  c1->cd(2);
  gr2->Draw("AL");
  c1->SaveAs(filename.c_str());

  delete c1;
};

void save16Graphs(vector<TGraph*> grs, string filename){
  
  TCanvas* c1 = new TCanvas("c1", "", 1280, 960);
  c1->Divide(4,4);
  //  TGraph* gr;
  for (int i = 0; i < 16; i++){
    c1->cd(i+1);
    //    gr = (TGraph*)grs->At(i);
    
    grs[i]->Draw("AL");
    //    delete gr;
  }
  c1->SaveAs(filename.c_str());
  
  
  delete c1;
}

void saveNGraphs(vector<TGraph*> grs, string filename, int numGraphs){

  int x_dim = (int)floor(sqrt((double)numGraphs));
  int y_dim = (int)ceil(sqrt((double)numGraphs));

  //  cout << x_dim << " : " << y_dim << endl;
  
  TCanvas* c1 = new TCanvas("c1", "", 1280, 960);
  c1->Divide(x_dim, y_dim);
  //  TGraph* gr;
  for (int i = 0; i < numGraphs; i++){
    c1->cd(i+1);
    //    gr = (TGraph*)grs->At(i);
    
    grs[i]->Draw("AL");
    //    delete gr;
  }
  c1->SaveAs(filename.c_str());
  
  
  delete c1;
}

void saveNGraphs_2Tlines(vector<TGraph*> grs, string filename, int numGraphs, vector<double> averages, vector<double> deviations){

  TLine *line_average[numGraphs];
  TLine *line_deviation[numGraphs];

  int x_dim = (int)floor(sqrt((double)numGraphs));
  int y_dim = (int)ceil(sqrt((double)numGraphs));

  //  cout << x_dim << " : " << y_dim << endl;
  
  TCanvas* c1 = new TCanvas("c1", "", 1280, 960);
  c1->Divide(x_dim, y_dim);
  //  TGraph* gr;
  for (int i = 0; i < numGraphs; i++){
    c1->cd(i+1);
    //    gr = (TGraph*)grs->At(i);
    line_average[i] = new TLine(-1000., averages[i], 1000., averages[i]);
    line_deviation[i]  = new TLine(-1000., averages[i]+deviations[i], 1000., averages[i]+deviations[i]);
    
    grs[i]->Draw("AL");
    line_average[i]->Draw("same");
    line_deviation[i]->Draw("same");
    //    delete gr;
  }
  c1->SaveAs(filename.c_str());
  
  for (int i = 0; i < numGraphs; i++){
    delete line_average[i];
    delete line_deviation[i];
  }
  delete c1;
}



vector<TGraph*> makeGraphsFromRF(UsefulAtriStationEvent* realAtriEvPtr, int numGraphs, string xlabel, string ylabel, vector<string> titles){

  //  TObjArray* graphArray = new TObjArray();

  //const int numGraphsConst = numGraphs;
  vector<TGraph*> graphs;
  //  TGraph* graphs[numGraphsConst];

  for (int i = 0; i < numGraphs; i++){
    TGraph* gr = realAtriEvPtr->getGraphFromRFChan(i);
    gr->GetXaxis()->SetTitle(xlabel.c_str());
    gr->GetYaxis()->SetTitle(ylabel.c_str());
    gr->SetTitle(titles[i].c_str());
    graphs.push_back(gr);
  }
  
  return graphs;
}

void deleteGraphVector(vector<TGraph*> graphs){
  for (int i=0; i < graphs.size(); i++){
    delete graphs[i];
  }
  graphs.clear();
}

vector<TGraph*> makeInterpolatedGraphs(vector<TGraph*> graphsIn, double intTimestep, string xlabel, string ylabel, vector<string> titles){

  vector<TGraph*> graphsOut;
  for (int i = 0; i < graphsIn.size(); i++){
    TGraph* gr = FFTtools::getInterpolatedGraph(graphsIn[i], intTimestep);
    gr->GetXaxis()->SetTitle(xlabel.c_str());
    gr->GetYaxis()->SetTitle(ylabel.c_str());
    gr->SetTitle(titles[i].c_str());
    graphsOut.push_back(gr);
  }

  return  graphsOut;
}

vector<TGraph*> makePaddedGraphs(vector<TGraph*> graphsIn, bool isSoftTrigger, string xlabel, string ylabel, vector<string> titles){

  vector<TGraph*> graphsOut;
  int maxPoints = 0;
  for (int i = 0; i < graphsIn.size(); i++){
    int nPoints=graphsIn[i]->GetN();
    if (nPoints > maxPoints){ maxPoints = nPoints;}
  }
  int nPointsNew = (int)pow(2, ceil(log2((double)maxPoints)));
  //  if (isSoftTrigger == false){ nPointsNew = 1024;} else {nPointsNew = 512;}
  if (isSoftTrigger == false){ nPointsNew = 1024;} else {nPointsNew = 1024;}

  for (int i = 0; i < graphsIn.size(); i++){
    TGraph* gr = new TGraph();
    double x, y;
    double x0, x1;
    int graphSize = graphsIn[i]->GetN();
    if (isSoftTrigger == false && graphSize > 1024)
      { graphSize = 1024;} 
    //    if (isSoftTrigger == true && graphSize > 512) 
      //      {graphSize = 512;}    
    if (isSoftTrigger == true && graphSize > 1024) 
      {graphSize = 1024;}    

    for (int point = 0; point < graphSize; point++){
      graphsIn[i]->GetPoint(point, x, y);
      if (point == 0){ x0 = x; }
      if (point == 1){ x1 = x; }
      gr->SetPoint(gr->GetN(), x, y);
    }
    double dt = x1-x0;

    while (gr->GetN() < nPointsNew){
      gr->GetPoint(gr->GetN()-1, x, y);
      double x_new = x+dt;
      gr->SetPoint(gr->GetN(), x_new, 0.);
    }


    gr->GetXaxis()->SetTitle(xlabel.c_str());
    gr->GetYaxis()->SetTitle(ylabel.c_str());
    gr->SetTitle(titles[i].c_str());
    graphsOut.push_back(gr);
  }

  return  graphsOut;
}


vector<TGraph*> makePowerSpectrumGraphs(vector<TGraph*> graphsIn, string xlabel, string ylabel, vector<string> titles){

  vector<TGraph*> graphsOut;
  for (int i = 0; i < graphsIn.size(); i++){
    TGraph* gr = FFTtools::makePowerSpectrumMilliVoltsNanoSeconds(graphsIn[i]);
    gr->GetXaxis()->SetTitle(xlabel.c_str());
    gr->GetYaxis()->SetTitle(ylabel.c_str());
    gr->SetTitle(titles[i].c_str());
    graphsOut.push_back(gr);
  }

  return  graphsOut;
}

vector<TGraph*> makeVoltageSpectrumGraphs(vector<TGraph*> graphsIn, string xlabel, string ylabel, vector<string> titles){

  vector<TGraph*> graphsOut;
  for (int i = 0; i < graphsIn.size(); i++){
    TGraph* gr = makeSpectrum_mVPerRootHz(graphsIn[i]);
    gr->GetXaxis()->SetTitle(xlabel.c_str());
    gr->GetYaxis()->SetTitle(ylabel.c_str());
    gr->SetTitle(titles[i].c_str());
    graphsOut.push_back(gr);
  }

  return  graphsOut;
}


vector<TGraph*> makePhaseGraphs(vector<TGraph*> graphsIn, string xlabel, string ylabel, vector<string> titles){

  vector<TGraph*> graphsOut;
  for (int i = 0; i < graphsIn.size(); i++){
    TGraph* gr = makeRawPhase(graphsIn[i]);
    gr->GetXaxis()->SetTitle(xlabel.c_str());
    gr->GetYaxis()->SetTitle(ylabel.c_str());
    gr->SetTitle(titles[i].c_str());
    graphsOut.push_back(gr);
  }

  return  graphsOut;
}


void saveHistOverlap(vector<TH1D*> hists, string filename){
  gStyle->SetOptStat(0);

  TCanvas* c1 = new TCanvas("c1", "", 640, 480);
  for (int i = 0; i < hists.size(); i++){

    if (i ==0 ){
      hists[i]->Draw();
    }
    else {
      hists[i]->Draw("same");
    }
    
    if (i == 0){
      hists[i]->SetLineColor(kBlack);
    }
    if (i == 1){
      hists[i]->SetLineColor(kRed);
    }
    if (i == 2){
      hists[i]->SetLineColor(kBlue);
    }
    if (i == 3){
      hists[i]->SetLineColor(kGreen);
    }
    if (i == 4){
      hists[i]->SetLineColor(kViolet);
    }
  }

  c1->SaveAs(filename.c_str());

  delete c1;
};

void setHistLineColors(vector<TH1D*> hists);

void setHistLineColors(vector<TH1D*> hists){
  for (int i = 0; i < hists.size(); i++){
    if (i == 0){
      hists[i]->SetLineColor(kBlack);
    }
    if (i == 1){
      hists[i]->SetLineColor(kRed);
    }
    if (i == 2){
      hists[i]->SetLineColor(kBlue);
    }
    if (i == 3){
      hists[i]->SetLineColor(kGreen);
    }
    if (i == 4){
      hists[i]->SetLineColor(kViolet);
    }
  }

};

void saveHistOverlapPolarized(vector<TH1D*> histsV, vector<TH1D*> histsH, string filename, int logOpt=0);
void saveHistOverlapPolarized(vector<TH1D*> histsV, vector<TH1D*> histsH, string filename, int logOpt){

  TCanvas* c1 = new TCanvas("c1", "", 640, 960);

  c1->Divide(1,2);
  c1->cd(1);

  if (logOpt == 1){
    gPad->SetLogx();
  }
  if (logOpt == 2){
    gPad->SetLogy();
  }
  if (logOpt == 3){
    gPad->SetLogx();
    gPad->SetLogy();
  }
  setHistLineColors(histsV);

  for (int i = 0; i < histsV.size(); i++){
    if (i ==0 ){
      histsV[i]->SetMinimum(0.1);
      histsV[i]->Draw();
    }
    else {
      histsV[i]->Draw("same");
    }
  }
  c1->cd(2);
  if (logOpt == 1){
    gPad->SetLogx();
  }
  if (logOpt == 2){
    gPad->SetLogy();
  }
  if (logOpt == 3){
    gPad->SetLogx();
    gPad->SetLogy();
  }
  setHistLineColors(histsH);
  for (int i = 0; i < histsH.size(); i++){
    if (i ==0 ){
      histsH[i]->SetMinimum(0.1);
      histsH[i]->Draw();
    }
    else {
      histsH[i]->Draw("same");
    }
  }

  c1->SaveAs(filename.c_str());

  delete c1;
}


void deleteHistVector(vector<TH1D*> hists){
  for (int i=0; i < hists.size(); i++){
    delete hists[i];
  }
  hists.clear();
}

vector<TH1D*> initializeHistVector(int size, string name, int bins, double lowBound, double upperBound){
  vector<TH1D*> vHist;
  stringstream ss;

  for (int i = 0; i < size; i++){
    ss.str("");
    ss << name << "_" << i;
    TH1D* hist = new TH1D(ss.str().c_str(), "", bins, lowBound,upperBound);
    vHist.push_back(hist);
  }
  return vHist;
}


int getRunNum(char* runfile){
  string file = string(runfile);
  string chRun = "Run";
  size_t foundRun=file.find(chRun);
  string strRunNum = file.substr(foundRun + 3, foundRun + 9);
  int runNum = atoi(strRunNum.c_str());
  return runNum;
}

int getrunNum(char* runfile){
  string file = string(runfile);
  string chRun = "run";
  size_t foundRun=file.find(chRun);
  string strRunNum = file.substr(foundRun + 3, foundRun + 9);
  int runNum = atoi(strRunNum.c_str());
  return runNum;
}

void saveHist2D(TH2* hist, string filename){

  TCanvas* c1 = new TCanvas("c1", "", 1280, 960);
  hist->Draw("colz");
  c1->SaveAs(filename.c_str());

  delete c1;
};

void saveHist2D(TH2* hist, string filename, int logOpt){

  TCanvas* c1 = new TCanvas("c1", "", 1280, 960);
  if (logOpt > 0){
    gPad->SetLogz();
  }
  hist->SetMinimum(0.);
  hist->SetMaximum(0.35);
  gStyle->SetOptStat(0);
  hist->Draw("colz");
  c1->SaveAs(filename.c_str());

  delete c1;
};


vector<TH2D*> initializeHistVector2D(int size, string name, int bins_x, double lowBound_x, double upperBound_x, int bins_y, double lowBound_y, double upperBound_y){
  vector<TH2D*> vHist;
  stringstream ss;

  for (int i = 0; i < size; i++){
    ss.str("");
    ss << name << "_" << i;
    TH2D* hist = new TH2D(ss.str().c_str(), "", bins_x, lowBound_x, upperBound_x, bins_y, lowBound_y, upperBound_y);
    vHist.push_back(hist);
  }
  return vHist;
}

void saveNHists2D(vector<TH2D*> hists, string filename_base){

  stringstream ss;

  for (int i = 0; i < hists.size(); i++){
    TCanvas* c1 = new TCanvas("c1", "", 1280, 960);
    hists[i]->Draw("colz");
    ss.str("");
    ss << filename_base << "_" << i << ".pdf";
    c1->SetLogz();
    c1->SaveAs(ss.str().c_str());
    delete c1;
  }
}

string getProcessedFilename(int stationID, char* outputdir, char* runfile ){
  
  stringstream ss;
  int runNum = getrunNum(runfile);
  ss.str("");
  ss << outputdir << "ProcessedFile/processed_station_" << stationID << "_run_" << runNum << ".root";

  return ss.str();
}

string getProcessedFilename_filter(int stationID, char* outputdir, char* runfile ){
  
  stringstream ss;
  int runNum = getrunNum(runfile);
  ss.str("");
  ss << outputdir << "processed_station_" << stationID << "_run_" << runNum << "_filter.root";

  return ss.str();
}

string getProcessedFilename_recoRadius(int stationID, char* outputdir, char* runfile, int radius){
  
  stringstream ss;
  int runNum = getrunNum(runfile);
  ss.str("");
  ss << outputdir << "processed_station_" << stationID << "_run_" << runNum << "_recoRadius_" << radius << ".root";

  return ss.str();
}


string getNoiseFilename(int stationID, char* outputdir, char* runfile ){
  
  stringstream ss;
  int runNum = getrunNum(runfile);
  ss.str("");
  ss << outputdir << "Noise/processed_station_" << stationID << "_run_" << runNum << ".root";

  return ss.str();
}


string getRunSummaryFilename(int stationID, char* outputdir, char* runfile ){
  
  stringstream ss;
  int runNum = getrunNum(runfile);
  ss.str("");
  ss << outputdir << "RunSummary/run_summary_station_" << stationID << "_run_" << runNum << ".root";

  return ss.str();
}

vector<TGraph*> initializeGraphVector(int size){
  vector<TGraph*> vGr;

  for (int i = 0; i < size; i++){
    TGraph* gr = new TGraph();
    vGr.push_back(gr);
  }

  return vGr;
}


vector<TGraph*> makeDifferenceGraphs(vector<TGraph*> graphsIn, int graphBegin, int graphEnd, string xlabel, string ylabel, vector<string> titles){

  stringstream ss_temp;
  vector<TGraph*> graphsOut;
  for (int i = graphBegin; i < graphEnd; i++){
    for (int j = i+1; j < graphEnd+1; j++){
      
      TGraph* gr = getGraphDifference(graphsIn[i], graphsIn[j]);
      gr->GetXaxis()->SetTitle(xlabel.c_str());
      gr->GetYaxis()->SetTitle(ylabel.c_str());
      ss_temp.str("");
      ss_temp << "Difference between " << titles[i] << " and " << titles[j];
      gr->SetTitle(ss_temp.str().c_str());
      graphsOut.push_back(gr);
    }
  }

  return  graphsOut;
}

vector<TGraph*> makeSumGraphs(vector<TGraph*> graphsIn, int graphBegin, int graphEnd, string xlabel, string ylabel, vector<string> titles){

  stringstream ss_temp;
  vector<TGraph*> graphsOut;
  for (int i = graphBegin; i < graphEnd; i++){
    for (int j = i+1; j < graphEnd+1; j++){
      
      TGraph* gr = getGraphSum(graphsIn[i], graphsIn[j]);
      gr->GetXaxis()->SetTitle(xlabel.c_str());
      gr->GetYaxis()->SetTitle(ylabel.c_str());
      ss_temp.str("");
      ss_temp << "Sum of " << titles[i] << " and " << titles[j];
      gr->SetTitle(ss_temp.str().c_str());
      graphsOut.push_back(gr);
    }
  }

  return  graphsOut;
}


#endif
