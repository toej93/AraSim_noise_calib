////////////////////////////////////////////////////////////////////////////// 
/////  WaveformFns.h                                                     /////  
/////                                                                    /////  
/////  Description:                                                      ///// 
/////     Functions for making a variety of cuts on Ara Data             ///// 
/////  Author: CGP, adapted from ACG, class structure from RJN           ///// 
////////////////////////////////////////////////////////////////////////////// 

#ifndef WAVEFORMFNS_H
#define WAVEFORMFNS_H

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

#include "Constants.h"
//#include "PlottingFns.h"

double getRMS( TGraph *plot, int numPointsToInclude);
double getRMS( double *array, int bin);
void getAbsMaximum(TGraph *plot, double &x_max, double &y_max);
double integrateBinPower( TGraph *plot, int numBinsToIntegrate, vector<double> &integratedBins);
vector<TGraph*> makeIntegratedBinPowerGraphs(vector<TGraph*> graphsInput, int numBinsToIntegrate, string xlabel, string ylabel, vector<string> titles);
void getAbsMaximum(vector<TGraph*> graphs, vector<double> &xs, vector<double> &ys );
vector<vector<vector<vector<int> > > > setupPairs();


void getMaximum(TGraph * plot, double &x_max, double &y_max){
  int nPoints = plot->GetN();

  double x_temp, y_temp;
  double y_good = -9.0E99;
  double x_good = 0;
  int test = 0;

  if (nPoints > 0){
    test = plot->GetPoint(0, x_temp, y_temp);
    x_good = x_temp;
    y_good = y_temp;
  }

  for (int i = 0; i < nPoints; i++){
    test = plot->GetPoint(i, x_temp, y_temp);
    if (y_temp > y_good){
      y_good = y_temp;
      x_good = x_temp;
    }
  }

  x_max = x_good;
  y_max = y_good;

  return;
}



void getMinimum(TGraph * plot, double &x_min, double &y_min){
  int nPoints = plot->GetN();

  double x_temp, y_temp;
  double y_good = 9.0E99;
  double x_good = 0;
  int test = 0;
  
  if (nPoints > 0){
    test = plot->GetPoint(0, x_temp, y_temp);
    x_good = x_temp;
    y_good = y_temp;
  }

  for (int i = 0; i < nPoints; i++){
    test = plot->GetPoint(i, x_temp, y_temp);
    if (y_temp < y_good){
      y_good = y_temp;
      x_good = x_temp;
    }
  }

  x_min = x_good;
  y_min = y_good;

  return;
}




void getAbsMaximum( TGraph *plot, double &x_max, double &y_max )
{
  int nPoints = plot->GetN();

  double x_temp, y_temp;
  double y_good = -9.0E99;
  double x_good = 0;
  int test = 0;

  for (int i = 0; i < nPoints; i++){
    test = plot->GetPoint(i, x_temp, y_temp);
    if (abs(y_temp) > y_good){
      y_good = abs(y_temp);
      x_good = x_temp;
    }
  }

  x_max = x_good;
  y_max = y_good;

  return;
}

void getPeaksAboveThreshold(TGraph * plot, double threshold, vector<double> &x_max, vector<double> &y_max){
  int nPoints = plot->GetN();
  double x_temp, y_temp;
  int test = 0;

  for (int i = 0; i < nPoints; i++){
    test = plot->GetPoint(i, x_temp, y_temp);
    if (y_temp > threshold){
      x_max.push_back(x_temp);
      y_max.push_back(y_temp);
    }
  }

  return;
}




double getRMS(double *array, int bin) {

  double x,y;
  double RMSVal = 0.;
  for(int i=0;i<bin;i++) {
    RMSVal += (array[i]*array[i]);
  }

  return sqrt(RMSVal/((double)bin));

}

double getRMS( TGraph *plot, int numPointsToInclude)
{
  int nPoints = plot->GetN();
  //  Double_t *xVals = plot->GetX();                                                                                                          
  Double_t *yVals = plot->GetY();

  int pointsToAdd;
  if (numPointsToInclude > nPoints || numPointsToInclude == 0){
    pointsToAdd = nPoints;
  } else {
    pointsToAdd = numPointsToInclude;
  }

  double RMS = getRMS(yVals, pointsToAdd);

  return RMS;
}

double getWaveformVariance(double *array, int bin) {

  double x,y;
  double variance = 0.;
  for(int i=0;i<bin;i++) {
    variance += (array[i]*array[i]);
  }

  return variance;

}


double getWaveformVariance( TGraph *plot, int numPointsToInclude)
{
  int nPoints = plot->GetN();
  //  Double_t *xVals = plot->GetX();                                                                                                          
  Double_t *yVals = plot->GetY();

  int pointsToAdd;
  if (numPointsToInclude > nPoints || numPointsToInclude == 0){
    pointsToAdd = nPoints;
  } else {
    pointsToAdd = numPointsToInclude;
  }

  double variance = getWaveformVariance(yVals, pointsToAdd);

  return variance;
}


double integrateBinPower( TGraph *plot, int numBinsToIntegrate, vector<double> &integratedBins)
{
  int nPoints = plot->GetN();
  if (nPoints < numBinsToIntegrate){
    return 0;
  }

  //  Double_t *xVals = plot->GetX();                                                                                                          
  Double_t *yVals = plot->GetY();
  std::deque<double> integrator;
  double sum = 0.;
  for (int i = 0; i < numBinsToIntegrate; i++){
    integrator.push_back(pow(yVals[i], 2));
    sum = accumulate(integrator.begin(), integrator.end(), 0);
  }
  double max = 0.;
  integratedBins.push_back(sum);

  for (int i = 0+numBinsToIntegrate; i < nPoints; i++){

    sum = sum - integrator[0];
    integrator.pop_front();
    integrator.push_back(pow(yVals[i], 2));
    sum = sum + integrator[numBinsToIntegrate-1];

    integratedBins.push_back(sum);

    if (sum > max){
      max = sum;
    }
  }

  return max;
}

double integrateAllPower(TGraph *plot){
  int nPoints = plot->GetN();
  double sum = 0.;
  double *yVals = plot->GetY();
  for (int i = 0; i < nPoints; i++){
    sum = sum +  pow(yVals[i],2);
  }
  return sum;
}

double getAveragePower(TGraph *plot){
  int nPoints = plot->GetN();
  double average = integrateAllPower(plot)/nPoints;
  return average;
}

vector<int> getDividers(TGraph *plot, int divisions){
  int nPoints = plot->GetN();
  vector<int> dividers;
  dividers.push_back(0);
  for (int i = 0; i < divisions; i++){
    dividers.push_back( (i+1)*nPoints/divisions);
  }
  return dividers;  
}

double integrate(TGraph *plot, int firstBin, int lastBin){
  int nPoints = plot->GetN();
  if (lastBin > nPoints){
    cerr << "ERROR: Cannot integrate graph, upper limit not valid" << endl;
    return -1.;
  }
  double sum = 0.;
  double *yVals = plot->GetY();
  for (int i = firstBin; i < lastBin; i++){
    sum = sum +  yVals[i];
  }
  return sum;
}

double getAverage(TGraph *plot, int firstBin, int lastBin){
  double average = integrate(plot, firstBin, lastBin)/(double)(lastBin-firstBin);
  return average;  
}

double getDeviation(TGraph *plot, int firstBin, int lastBin, double average){
  double sumOfSquares = 0.;
  double *yVals = plot->GetY();
  for (int i = firstBin; i < lastBin; i++){
    sumOfSquares = sumOfSquares +  pow((yVals[i]-average), 2);
  }
  double variance = sqrt(sumOfSquares/(double)(lastBin-firstBin));
  return variance;  
}

void getAverageDeviation_Divided(TGraph *plot, int divisions, double &average_out, double &deviation_out){

  int i_lowestAverage = -1;
  int i_secondLowestAverage = -1;
  double lowestAverage = 1E99;
  double secondLowestAverage = 1E99;

  vector<int> dividers = getDividers(plot, divisions);
  vector<double> averages;
  vector<double> variances;
  for (int i = 0; i < divisions; i++){

    double average = getAverage(plot, dividers[i], dividers[i+1]);
    averages.push_back(average);
    double variance = getDeviation(plot, dividers[i], dividers[i+1], average);
    variances.push_back(variance);

    //    cout << plot->GetN() << " : " << dividers[i+1] << endl;
    //    cout << average << " : " << variance << endl;

    if (average < lowestAverage){
      lowestAverage = average;
      i_secondLowestAverage = i_lowestAverage;
      i_lowestAverage = i;
    }
    if (average < secondLowestAverage && average > lowestAverage){
      secondLowestAverage = average;
      i_secondLowestAverage = i;
    }
  }

  /*
  cout << i_lowestAverage << " : " << i_secondLowestAverage << endl;
  cout << averages[i_lowestAverage] << " : " << averages[i_secondLowestAverage] << endl;
  cout<< endl;
  */

  average_out = (averages[i_lowestAverage]+averages[i_secondLowestAverage])/2.;
  deviation_out = (variances[i_lowestAverage]+variances[i_secondLowestAverage])/2.;
}

void getVAveragesDeviations_Divided(vector<TGraph*> vPlots, int divisions, vector<double> &vAverages, vector<double> &vDeviations){
  //vector<double> vDeviations;
  //vector<double> vAverages;
  vAverages.clear();
  vDeviations.clear();

  for (int i = 0; i< vPlots.size(); i++){
    double average;
    double deviation;
    getAverageDeviation_Divided(vPlots[i], divisions, average, deviation);
    vAverages.push_back(average);
    vDeviations.push_back(deviation);
  }
}

double integratePowerFirstNBins(TGraph *plot, int nBins){
  int nPoints = plot->GetN();
  if (nPoints < nBins){
    nBins = nPoints;
  }
  double sum = 0.;
  double *yVals = plot->GetY();
  for (int i = 0; i < nBins; i++){
    sum = sum +  pow(yVals[i],2);
  }
  return sum;
}

double getAveragePowerFirstNBins(TGraph *plot, int nBins){
  int nPoints = plot->GetN();
  if (nPoints < nBins){
    nBins = nPoints;
  }
  double average = integratePowerFirstNBins(plot, nBins)/nBins;
  return average;
}


vector<TGraph*> makeIntegratedBinPowerGraphs(vector<TGraph*> graphsInput, int numBinsToIntegrate, string xlabel, string ylabel, vector<string> titles){

  //  cout << graphsInput.size() << endl;
  vector<TGraph*> graphsOutput;

  for (int i = 0; i < graphsInput.size(); i++){
    vector<double> integratedBins;
    double maxIntPower = integrateBinPower(graphsInput[i], numBinsToIntegrate, integratedBins);
    double *volts = &integratedBins[0];
    TGraph* gIntPower = new TGraph(integratedBins.size(), graphsInput[i]->GetX(), volts);
    gIntPower->GetXaxis()->SetTitle(xlabel.c_str());
    gIntPower->GetYaxis()->SetTitle(ylabel.c_str());
    gIntPower->SetTitle(titles[i].c_str());
    graphsOutput.push_back(gIntPower);
    integratedBins.clear();
    //    delete gIntPower;
  }

  return graphsOutput;
}

vector<TGraph*> makeAvgVoltageFromIntPower(vector<TGraph*> graphsInput, int numBinsToIntegrate, string xlabel, string ylabel, vector<string> titles){

  vector<TGraph*> graphsOutput;
  // graphsOutput.resize(graphsInput.size());
  
  for (int i = 0; i < graphsInput.size(); i++){
    int nPoints = graphsInput[i]->GetN();
    double* xIn = graphsInput[i]->GetX();
    double* yIn = graphsInput[i]->GetY();

    TGraph* newGraph = new TGraph();
    for (int point = 0; point < nPoints; point++){
      double sqrtY;
      if (yIn[point]>0.){
	sqrtY= sqrt(yIn[point]/(double)numBinsToIntegrate);
      } else {
	sqrtY = 0.;
      }
      newGraph->SetPoint(newGraph->GetN(), xIn[point], sqrtY);
    }

    newGraph->GetXaxis()->SetTitle(xlabel.c_str());
    newGraph->GetYaxis()->SetTitle(ylabel.c_str());
    newGraph->SetTitle(titles[i].c_str());
    graphsOutput.push_back(newGraph);
    //    delete gIntPower;
  }

  return graphsOutput;
}




void getAbsMaximum(vector<TGraph*> graphs, vector<double> &xs, vector<double> &ys ){

  //  vector<double> xs;
  // vector<double> ys;

  xs.clear();
  ys.clear();

  double x, y;

  for (int i = 0; i < graphs.size(); i++){
    getAbsMaximum(graphs[i], x, y);
    xs.push_back(x);
    ys.push_back(y);
  }
}

void getRMS(vector<TGraph*> graphs, vector<double> &vRMS, int numPointsToInclude){

  vRMS.clear();
  double RMS;

  for (int i = 0; i < graphs.size(); i++){
    int numPoints_temp = numPointsToInclude;
    if(numPointsToInclude == 0){
      numPoints_temp = graphs[i]->GetN();
    }

    RMS = getRMS(graphs[i], numPoints_temp);
    vRMS.push_back(RMS);
  }
}

void getWaveformVarianceSamples(vector<TGraph*> graphs, int numPointsToInclude, vector<double> &vVariance, vector<double> &vSamples){

  vVariance.clear();
  vSamples.clear();

  double variance;

  for (int i = 0; i < graphs.size(); i++){
    int numPoints_temp = numPointsToInclude;
    if(numPointsToInclude == 0){
      numPoints_temp = graphs[i]->GetN();
    }

    variance = getWaveformVariance(graphs[i], numPoints_temp);
    vVariance.push_back(variance);
    vSamples.push_back(numPoints_temp);
  }
}

void addVarianceSamples(vector<double> vWaveformVariance, vector<int> vSamples, vector<double> vWaveformVarianceTotal, vector<long int> vSamplesTotal)
{
  for (int i = 0; i < vWaveformVariance.size(); i++){
    vWaveformVarianceTotal[i] = vWaveformVarianceTotal[i] + vWaveformVariance[i];
    vSamplesTotal[i] = vSamplesTotal[i] + vSamples[i];
  }
}



vector<double> getTotalPowers(vector<TGraph*> graphs){
  vector<double> powers;
  double power;
  for (int i = 0; i < graphs.size(); i++){
    power = integrateAllPower(graphs[i]);
    powers.push_back(power);
 }
  return powers;
}

vector<double> getAveragePowers(vector<TGraph*> graphs){
  vector<double> powers;
  double power;
  for (int i = 0; i < graphs.size(); i++){
    power = getAveragePower(graphs[i]);
    powers.push_back(power);
 }
  return powers;
}

vector<double> getAveragePowersFirstNBins(vector<TGraph*> graphs, int nBins){
  vector<double> powers;
  double power;
  for (int i = 0; i < graphs.size(); i++){
    power = getAveragePowerFirstNBins(graphs[i], nBins);
    powers.push_back(power);
 }
  return powers;
}


vector<vector<vector<vector<int> > > > setupPairs(int stationId){

  vector<vector<vector<vector<int> > > > Pairs;
  vector<vector<vector<int> > > pairsType;
  vector<vector<int> > pairsPol; 
  vector<int> pair;

  if (stationId == 2){
    //Vertical Pairs
    pair.push_back(0); pair.push_back(4);
    pairsPol.push_back(pair); pair.clear();
    pair.push_back(1); pair.push_back(5);
    pairsPol.push_back(pair); pair.clear();
    pair.push_back(2); pair.push_back(6);
    pairsPol.push_back(pair); pair.clear();
    pair.push_back(3); pair.push_back(7);
    pairsPol.push_back(pair); pair.clear();
    pairsType.push_back(pairsPol);pairsPol.clear();
    
    pair.push_back(8); pair.push_back(12);
    pairsPol.push_back(pair); pair.clear();
    pair.push_back(9); pair.push_back(13);
    pairsPol.push_back(pair); pair.clear();
    pair.push_back(10); pair.push_back(14);
    pairsPol.push_back(pair); pair.clear();
    pair.push_back(11); pair.push_back(15);
    pairsPol.push_back(pair); pair.clear();
    pairsType.push_back(pairsPol);pairsPol.clear();
    
    Pairs.push_back(pairsType);pairsType.clear();
    
    // Side pairs 1
    
    pair.push_back(0); pair.push_back(1);
    pairsPol.push_back(pair); pair.clear();
    pair.push_back(4); pair.push_back(5);
    pairsPol.push_back(pair); pair.clear();
    pair.push_back(2); pair.push_back(3);
    pairsPol.push_back(pair); pair.clear();
    pair.push_back(6); pair.push_back(7);
    pairsPol.push_back(pair); pair.clear();
    pairsType.push_back(pairsPol);pairsPol.clear();
    
    pair.push_back(8); pair.push_back(9);
    pairsPol.push_back(pair); pair.clear();
    pair.push_back(12); pair.push_back(13);
    pairsPol.push_back(pair); pair.clear();
    pair.push_back(10); pair.push_back(11);
    pairsPol.push_back(pair); pair.clear();
    pair.push_back(14); pair.push_back(15);
    pairsPol.push_back(pair); pair.clear();
    pairsType.push_back(pairsPol);pairsPol.clear();
    
    Pairs.push_back(pairsType);pairsType.clear();
    
    // Side pairs 2
    
    pair.push_back(0); pair.push_back(2);
    pairsPol.push_back(pair); pair.clear();
    pair.push_back(4); pair.push_back(6);
    pairsPol.push_back(pair); pair.clear();
    pair.push_back(1); pair.push_back(3);
    pairsPol.push_back(pair); pair.clear();
    pair.push_back(5); pair.push_back(7);
    pairsPol.push_back(pair); pair.clear();
    pairsType.push_back(pairsPol);pairsPol.clear();
    
    pair.push_back(8); pair.push_back(10);
    pairsPol.push_back(pair); pair.clear();
    pair.push_back(12); pair.push_back(14);
    pairsPol.push_back(pair); pair.clear();
    pair.push_back(9); pair.push_back(11);
    pairsPol.push_back(pair); pair.clear();
    pair.push_back(13); pair.push_back(15);
    pairsPol.push_back(pair); pair.clear();
    pairsType.push_back(pairsPol);pairsPol.clear();
    
    Pairs.push_back(pairsType);pairsType.clear();
    
    //Diagonal pairs 1
    
    pair.push_back(0); pair.push_back(3);
    pairsPol.push_back(pair); pair.clear();
    pair.push_back(4); pair.push_back(7);
    pairsPol.push_back(pair); pair.clear();
    pairsType.push_back(pairsPol);pairsPol.clear();
    
    pair.push_back(8); pair.push_back(11);
    pairsPol.push_back(pair); pair.clear();
    pair.push_back(12); pair.push_back(15);
    pairsPol.push_back(pair); pair.clear();
    pairsType.push_back(pairsPol);pairsPol.clear();
    
    Pairs.push_back(pairsType);pairsType.clear();
    
    //Diagonal pairs 2
    
    pair.push_back(1); pair.push_back(2);
    pairsPol.push_back(pair); pair.clear();
    pair.push_back(5); pair.push_back(6);
    pairsPol.push_back(pair); pair.clear();
    pairsType.push_back(pairsPol);pairsPol.clear();
    
    pair.push_back(9); pair.push_back(10);
    pairsPol.push_back(pair); pair.clear();
    pair.push_back(13); pair.push_back(14);
    pairsPol.push_back(pair); pair.clear();
    pairsType.push_back(pairsPol);pairsPol.clear();
    
    Pairs.push_back(pairsType);pairsType.clear();

    //Diagonal pairs 1 - side
    
    pair.push_back(1); pair.push_back(4);
    pairsPol.push_back(pair); pair.clear();
    pair.push_back(3); pair.push_back(6);
    pairsPol.push_back(pair); pair.clear();
    pairsType.push_back(pairsPol);pairsPol.clear();
    
    pair.push_back(9); pair.push_back(12);
    pairsPol.push_back(pair); pair.clear();
    pair.push_back(11); pair.push_back(14);
    pairsPol.push_back(pair); pair.clear();
    pairsType.push_back(pairsPol);pairsPol.clear();
    
    Pairs.push_back(pairsType);pairsType.clear();

    //Diagonal pairs 2 - side
    
    pair.push_back(0); pair.push_back(5);
    pairsPol.push_back(pair); pair.clear();
    pair.push_back(2); pair.push_back(7);
    pairsPol.push_back(pair); pair.clear();
    pairsType.push_back(pairsPol);pairsPol.clear();
    
    pair.push_back(8); pair.push_back(13);
    pairsPol.push_back(pair); pair.clear();
    pair.push_back(10); pair.push_back(15);
    pairsPol.push_back(pair); pair.clear();
    pairsType.push_back(pairsPol);pairsPol.clear();
    
    Pairs.push_back(pairsType);pairsType.clear();


    //Diagonal pairs 3 - side
    
    pair.push_back(1); pair.push_back(7);
    pairsPol.push_back(pair); pair.clear();
    pair.push_back(0); pair.push_back(6);
    pairsPol.push_back(pair); pair.clear();
    pairsType.push_back(pairsPol);pairsPol.clear();
    
    pair.push_back(9); pair.push_back(15);
    pairsPol.push_back(pair); pair.clear();
    pair.push_back(8); pair.push_back(14);
    pairsPol.push_back(pair); pair.clear();
    pairsType.push_back(pairsPol);pairsPol.clear();
    
    Pairs.push_back(pairsType);pairsType.clear();

    //Diagonal pairs 4 - side
    
    pair.push_back(2); pair.push_back(4);
    pairsPol.push_back(pair); pair.clear();
    pair.push_back(3); pair.push_back(5);
    pairsPol.push_back(pair); pair.clear();
    pairsType.push_back(pairsPol);pairsPol.clear();
    
    pair.push_back(10); pair.push_back(12);
    pairsPol.push_back(pair); pair.clear();
    pair.push_back(11); pair.push_back(13);
    pairsPol.push_back(pair); pair.clear();
    pairsType.push_back(pairsPol);pairsPol.clear();
    
    Pairs.push_back(pairsType);pairsType.clear();


  }

  if (stationId == 3){
    //Vertical Pairs
    pair.push_back(1); pair.push_back(5);
    pairsPol.push_back(pair); pair.clear();
    pair.push_back(0); pair.push_back(4);
    pairsPol.push_back(pair); pair.clear();
    pair.push_back(2); pair.push_back(6);
    pairsPol.push_back(pair); pair.clear();
    pair.push_back(3); pair.push_back(7);
    pairsPol.push_back(pair); pair.clear();
    pairsType.push_back(pairsPol);pairsPol.clear();
    
    pair.push_back(9); pair.push_back(13);
    pairsPol.push_back(pair); pair.clear();
    pair.push_back(8); pair.push_back(12);
    pairsPol.push_back(pair); pair.clear();
    pair.push_back(10); pair.push_back(14);
    pairsPol.push_back(pair); pair.clear();
    pair.push_back(11); pair.push_back(15);
    pairsPol.push_back(pair); pair.clear();
    pairsType.push_back(pairsPol);pairsPol.clear();
    
    Pairs.push_back(pairsType);pairsType.clear();
    
    // Side pairs 1
    
    pair.push_back(1); pair.push_back(0);
    pairsPol.push_back(pair); pair.clear();
    pair.push_back(5); pair.push_back(4);
    pairsPol.push_back(pair); pair.clear();
    pair.push_back(2); pair.push_back(3);
    pairsPol.push_back(pair); pair.clear();
    pair.push_back(6); pair.push_back(7);
    pairsPol.push_back(pair); pair.clear();
    pairsType.push_back(pairsPol);pairsPol.clear();
    
    pair.push_back(9); pair.push_back(8);
    pairsPol.push_back(pair); pair.clear();
    pair.push_back(13); pair.push_back(12);
    pairsPol.push_back(pair); pair.clear();
    pair.push_back(10); pair.push_back(11);
    pairsPol.push_back(pair); pair.clear();
    pair.push_back(14); pair.push_back(15);
    pairsPol.push_back(pair); pair.clear();
    pairsType.push_back(pairsPol);pairsPol.clear();
    
    Pairs.push_back(pairsType);pairsType.clear();
    
    // Side pairs 2
    
    pair.push_back(1); pair.push_back(2);
    pairsPol.push_back(pair); pair.clear();
    pair.push_back(5); pair.push_back(6);
    pairsPol.push_back(pair); pair.clear();
    pair.push_back(0); pair.push_back(3);
    pairsPol.push_back(pair); pair.clear();
    pair.push_back(4); pair.push_back(7);
    pairsPol.push_back(pair); pair.clear();
    pairsType.push_back(pairsPol);pairsPol.clear();
    
    pair.push_back(9); pair.push_back(10);
    pairsPol.push_back(pair); pair.clear();
    pair.push_back(13); pair.push_back(14);
    pairsPol.push_back(pair); pair.clear();
    pair.push_back(8); pair.push_back(11);
    pairsPol.push_back(pair); pair.clear();
    pair.push_back(12); pair.push_back(15);
    pairsPol.push_back(pair); pair.clear();
    pairsType.push_back(pairsPol);pairsPol.clear();
    
    Pairs.push_back(pairsType);pairsType.clear();
    
    //Diagonal pairs 1 - through the middle
    
    pair.push_back(1); pair.push_back(3);
    pairsPol.push_back(pair); pair.clear();
    pair.push_back(5); pair.push_back(7);
    pairsPol.push_back(pair); pair.clear();
    pairsType.push_back(pairsPol);pairsPol.clear();
    
    pair.push_back(9); pair.push_back(11);
    pairsPol.push_back(pair); pair.clear();
    pair.push_back(13); pair.push_back(15);
    pairsPol.push_back(pair); pair.clear();
    pairsType.push_back(pairsPol);pairsPol.clear();
    
    Pairs.push_back(pairsType);pairsType.clear();
    
    //Diagonal pairs 2 - through the middle
    
    pair.push_back(0); pair.push_back(2);
    pairsPol.push_back(pair); pair.clear();
    pair.push_back(4); pair.push_back(6);
    pairsPol.push_back(pair); pair.clear();
    pairsType.push_back(pairsPol);pairsPol.clear();
    
    pair.push_back(8); pair.push_back(10);
    pairsPol.push_back(pair); pair.clear();
    pair.push_back(12); pair.push_back(14);
    pairsPol.push_back(pair); pair.clear();
    pairsType.push_back(pairsPol);pairsPol.clear();
    
    Pairs.push_back(pairsType);pairsType.clear();

    //Diagonal pairs 1 - side
    
    pair.push_back(0); pair.push_back(5);
    pairsPol.push_back(pair); pair.clear();
    pair.push_back(3); pair.push_back(6);
    pairsPol.push_back(pair); pair.clear();
    pairsType.push_back(pairsPol);pairsPol.clear();
    
    pair.push_back(8); pair.push_back(13);
    pairsPol.push_back(pair); pair.clear();
    pair.push_back(11); pair.push_back(14);
    pairsPol.push_back(pair); pair.clear();
    pairsType.push_back(pairsPol);pairsPol.clear();
    
    Pairs.push_back(pairsType);pairsType.clear();

    //Diagonal pairs 2 - side
    
    pair.push_back(1); pair.push_back(4);
    pairsPol.push_back(pair); pair.clear();
    pair.push_back(2); pair.push_back(7);
    pairsPol.push_back(pair); pair.clear();
    pairsType.push_back(pairsPol);pairsPol.clear();
    
    pair.push_back(9); pair.push_back(12);
    pairsPol.push_back(pair); pair.clear();
    pair.push_back(10); pair.push_back(15);
    pairsPol.push_back(pair); pair.clear();
    pairsType.push_back(pairsPol);pairsPol.clear();
    
    Pairs.push_back(pairsType);pairsType.clear();


    //Diagonal pairs 3 - side
    
    pair.push_back(0); pair.push_back(7);
    pairsPol.push_back(pair); pair.clear();
    pair.push_back(1); pair.push_back(6);
    pairsPol.push_back(pair); pair.clear();
    pairsType.push_back(pairsPol);pairsPol.clear();
    
    pair.push_back(8); pair.push_back(15);
    pairsPol.push_back(pair); pair.clear();
    pair.push_back(9); pair.push_back(14);
    pairsPol.push_back(pair); pair.clear();
    pairsType.push_back(pairsPol);pairsPol.clear();
    
    Pairs.push_back(pairsType);pairsType.clear();

    //Diagonal pairs 4 - side
    
    pair.push_back(2); pair.push_back(5);
    pairsPol.push_back(pair); pair.clear();
    pair.push_back(3); pair.push_back(4);
    pairsPol.push_back(pair); pair.clear();
    pairsType.push_back(pairsPol);pairsPol.clear();
    
    pair.push_back(10); pair.push_back(13);
    pairsPol.push_back(pair); pair.clear();
    pair.push_back(11); pair.push_back(12);
    pairsPol.push_back(pair); pair.clear();
    pairsType.push_back(pairsPol);pairsPol.clear();
    
    Pairs.push_back(pairsType);pairsType.clear();

  }


  return Pairs;
}

vector<vector<vector<double> > > getHitTimeDifferences(vector<double> hitTimes, vector<double> peakIntPowers, vector<vector<vector<vector<int> > > > Pairs){
  
  vector<vector<vector<double> > > delays;

  delays.resize(Pairs.size());
  for (int i = 0; i < Pairs.size(); i++){
    delays[i].resize(Pairs[i].size());
    for (int j = 0; j < Pairs[i].size(); j++){
      //      delays[i][j].resize(Pairs[i][j].size());
      for (int k = 0; k < Pairs[i][j].size(); k++){
	double delay = hitTimes[Pairs[i][j][k][1]] - hitTimes[Pairs[i][j][k][0]];
	delays[i][j].push_back(delay);
	//	cout << delay << endl;
      }
    }
  }

  return delays;
}


vector<vector<vector<double> > > getHitTimeWeight_LesserPower(vector<double> hitTimes, vector<double> peakIntPowers, vector<vector<vector<vector<int> > > > Pairs){
  
  vector<vector<vector<double> > > vvv;

  vvv.resize(Pairs.size());
  for (int i = 0; i < Pairs.size(); i++){
    vvv[i].resize(Pairs[i].size());
    for (int j = 0; j < Pairs[i].size(); j++){
      //      delays[i][j].resize(Pairs[i][j].size());
      for (int k = 0; k < Pairs[i][j].size(); k++){
	double weight = 0.;
	if ( peakIntPowers[Pairs[i][j][k][1]] > peakIntPowers[Pairs[i][j][k][0]]) {
	  weight = peakIntPowers[Pairs[i][j][k][0]];
	} else {
	  weight = peakIntPowers[Pairs[i][j][k][1]];
	}

	vvv[i][j].push_back(weight);
	
	//	cout << delay << endl;
      }
    }
  }

  return vvv;
}

vector<vector<vector<double> > > getHitTimeWeight_LesserFractionalPower(vector<double> hitTimes, vector<double> peakIntPowers, vector<double> totalPowers, vector<vector<vector<vector<int> > > > Pairs){
  
  vector<vector<vector<double> > > vvv;

  vvv.resize(Pairs.size());
  for (int i = 0; i < Pairs.size(); i++){
    vvv[i].resize(Pairs[i].size());
    for (int j = 0; j < Pairs[i].size(); j++){
      //      delays[i][j].resize(Pairs[i][j].size());
      for (int k = 0; k < Pairs[i][j].size(); k++){
	double weight = 0.;
	int index1 = Pairs[i][j][k][0];
	int index2 = Pairs[i][j][k][1];

	double fractionPower1 = peakIntPowers[index1]/totalPowers[index1];
	double fractionPower2 = peakIntPowers[index2]/totalPowers[index2];

	if ( fractionPower2 > fractionPower1) {
	  weight = fractionPower1;
	} else {
	  weight = fractionPower2;
	}

	vvv[i][j].push_back(weight);
	
	//	cout << delay << endl;
      }
    }
  }

  return vvv;
}

vector<vector<vector<double> > > getHitTimeWeight_LesserRelativePower(vector<double> hitTimes, vector<double> peakIntPowers, vector<double> averagePowers, vector<vector<vector<vector<int> > > > Pairs){
  
  vector<vector<vector<double> > > vvv;

  vvv.resize(Pairs.size());
  for (int i = 0; i < Pairs.size(); i++){
    vvv[i].resize(Pairs[i].size());
    for (int j = 0; j < Pairs[i].size(); j++){
      //      delays[i][j].resize(Pairs[i][j].size());
      for (int k = 0; k < Pairs[i][j].size(); k++){
	double weight = 0.;
	int index1 = Pairs[i][j][k][0];
	int index2 = Pairs[i][j][k][1];

	double relativePower1 = peakIntPowers[index1]/averagePowers[index1];
	double relativePower2 = peakIntPowers[index2]/averagePowers[index2];

	if ( relativePower2 > relativePower1) {
	  weight = relativePower1;
	} else {
	  weight = relativePower2;
	}

	vvv[i][j].push_back(weight);
	
	//	cout << delay << endl;
      }
    }
  }

  return vvv;
}

vector<vector<vector<double> > > getHitTimeWeight_LesserRelativeRMS(vector<double> hitTimes, vector<double> peakIntRMS, vector<double> waveformRMS, vector<vector<vector<vector<int> > > > Pairs){
  
  vector<vector<vector<double> > > vvv;

  vvv.resize(Pairs.size());
  for (int i = 0; i < Pairs.size(); i++){
    vvv[i].resize(Pairs[i].size());
    for (int j = 0; j < Pairs[i].size(); j++){
      //      delays[i][j].resize(Pairs[i][j].size());
      for (int k = 0; k < Pairs[i][j].size(); k++){
	double weight = 0.;
	int index1 = Pairs[i][j][k][0];
	int index2 = Pairs[i][j][k][1];

	double relativeRMS1 = peakIntRMS[index1]/waveformRMS[index1];
	double relativeRMS2 = peakIntRMS[index2]/waveformRMS[index2];

	if ( relativeRMS2 > relativeRMS1) {
	  weight = relativeRMS1;
	} else {
	  weight = relativeRMS2;
	}

	vvv[i][j].push_back(weight);
	
	//	cout << delay << endl;
      }
    }
  }

  return vvv;
}






vector<TH1D*> makeHitTimeDist(vector<vector<vector<double> > > delays, int polarization){
  vector<TH1D*> vhDelays;
  stringstream ss;
  
  const int numHists = delays.size();

  TH1D *hDelays[numHists];

  for (int i = 0; i < delays.size(); i++){
    ss.str("");    
    ss << "hDelays_" << polarization << "_" << i;
    hDelays[i] = new TH1D(ss.str().c_str(), "", 100, -400, 400);
    //  for (int j = 0; j < delays[i].size(); j++){
    for (int j = polarization; j < polarization+1; j++){
      for (int k = 0; k < delays[i][j].size(); k++){
	hDelays[i]->Fill(delays[i][j][k]);
      }
    }
    vhDelays.push_back(hDelays[i]);
    
  }

  return vhDelays;
}

double getMean(vector<double> v){
  double sum = accumulate(v.begin(), v.end(), 0.0);
  double mean = sum / (double)v.size();
  return mean;
}

		      
double getRms(vector<double> v){
  double mean = getMean(v);
  double sumOfSquares = 0.0;
  for (int i = 0; i < v.size(); i++){
    sumOfSquares = sumOfSquares + pow(v[i] - mean, 2);
  }
  double rms = sqrt(sumOfSquares / (double)v.size());
  return rms;
}

vector<vector<double> > getRmsDelays(vector<vector<vector<double> > > delays){
  
  vector<vector<double> > rms_all;
  rms_all.resize(delays.size());
  for (int i = 0; i < delays.size(); i++){
    for (int j = 0; j < delays[i].size(); j++){
      double rms = getRms(delays[i][j]);
      rms_all[i].push_back(rms);
    }
  }

  return rms_all;
}

double getRmsPol(vector<vector<vector<double> > > delays, int polarization){
  
  double totalRms = 0.0;
  for (int i = 0; i < delays.size(); i++){
    for (int j = polarization; j < polarization+1; j++){
      double rms = getRms(delays[i][j]);
      totalRms = totalRms + pow(rms, 2)*delays[i][j].size();
    }
  }
  //  cout << totalRms << endl;
  return totalRms;
}

double getWeightedSum(vector<double> v, vector<double> weights){
  if (v.size() != weights.size()){
    cerr << "ERROR: getWeightedSum(): value vector not same size as weight vector" << endl;
    return 0.0;
  }
  double sum = 0.0;
  for (int i = 0; i < v.size(); i++){
    sum = sum + v[i]*weights[i];
  }
  return sum;
}

double getWeightedMean(vector<double> v, vector<double> weights){
  if (v.size() != weights.size()){
    cerr << "ERROR: getWeightedMean(): value vector not same size as weight vector" << endl;
    return 0.0;
  }
  double sum = getWeightedSum(v, weights);
  double sumOfWeights = accumulate(weights.begin(), weights.end(), 0.0);
  double mean = sum / sumOfWeights;
  return mean;
}

double getWeightedVariance(vector<double> v, vector<double> weights){
  if (v.size() != weights.size()){
    cerr << "ERROR: getWeightedVariance(): value vector not same size as weight vector" << endl;
    return 0.0;
  }

double mean = getWeightedMean(v, weights);
  double sumOfSquares = 0.0;
  //  double sumOfWeights = accumulate(weights.begin(), weights.end(), 0.0);
  for (int i = 0; i < v.size(); i++){
    sumOfSquares = sumOfSquares + weights[i]*pow((v[i] - mean), 2);
  }
  return sumOfSquares;
}

double getWeightedRms(vector<double> v, vector<double> weights){
  double mean = getWeightedMean(v, weights);
  double sumOfSquares = 0.0;
  double sumOfWeights = accumulate(weights.begin(), weights.end(), 0.0);
  for (int i = 0; i < v.size(); i++){
    sumOfSquares = sumOfSquares + weights[i]*pow((v[i] - mean), 2);
  }
  double rms = sqrt(sumOfSquares / sumOfWeights);
  return rms;
}

double getWeightedRmsPol(vector<vector<vector<double> > > delays, vector<vector<vector<double> > > weights, int polarization){
  
  double totalVariance = 0.0;
  double totalWeight = 0.0;
  for (int i = 0; i < delays.size(); i++){
    for (int j = polarization; j < polarization+1; j++){
      //double sum = getWeightedSum(delays[i][j], weights[i][j]);
      double sumOfWeights = accumulate(weights[i][j].begin(), weights[i][j].end(), 0.0);
      totalVariance = totalVariance + getWeightedVariance(delays[i][j], weights[i][j]);
      totalWeight = totalWeight + sumOfWeights;
    }
  }
  //  cout << totalRms << endl;
  return sqrt(totalVariance/totalWeight);
}

vector<vector<vector<double> > > getVVVToNPower(vector<vector<vector<double> > > vvvIn, double power){

  vector<vector<vector<double> > > vvvOut;

  vvvOut.resize(vvvIn.size());
  for (int i = 0; i < vvvIn.size(); i++){
    vvvOut[i].resize(vvvIn[i].size());
    for (int j = 0; j < vvvIn[i].size(); j++){
      for (int k = 0; k < vvvIn[i][j].size(); k++){
	double weight = pow(vvvIn[i][j][k], power);
	vvvOut[i][j].push_back(weight);
	//	cout << vvvIn[i][j][k] << " : " << weight << endl;
	//	cout << delay << endl;
      }
    }
  }
  return vvvOut;
}



vector<vector<vector<double> > > getHitTimeMagnitude_Lesser(vector<double> hitTimes, vector<double> magnitude, vector<vector<vector<vector<int> > > > Pairs){
  
  vector<vector<vector<double> > > pairMagnitudes;

  pairMagnitudes.resize(Pairs.size());
  for (int i = 0; i < Pairs.size(); i++){
    pairMagnitudes[i].resize(Pairs[i].size());
    for (int j = 0; j < Pairs[i].size(); j++){
      //      delays[i][j].resize(Pairs[i][j].size());
      for (int k = 0; k < Pairs[i][j].size(); k++){
	double pairMagnitude;
	if (magnitude[Pairs[i][j][k][1]] > magnitude[Pairs[i][j][k][0]]){
	  pairMagnitude = magnitude[Pairs[i][j][k][0]];
	}
	else {
	  pairMagnitude = magnitude[Pairs[i][j][k][1]];
	}
	pairMagnitudes[i][j].push_back(pairMagnitude);
	//	cout << delay << endl;
      }
    }
  }

  return pairMagnitudes;
}



bool isUsablePattern(vector<double> magnitudes, double threshold){
  bool isUsable = false;
  int counter = 0;
  for (int i = 0; i < magnitudes.size(); i++){
    if (magnitudes[i] >= threshold){
      counter++;
    }
  }
  if (counter > 1){
    isUsable = true;
  }

  return isUsable;
}



double getWeightedSum_Threshold(vector<double> v, vector<double> weights, vector<double> magnitudes, double threshold){
  if (v.size() != weights.size()){
    cerr << "ERROR: getWeightedSum(): value vector not same size as weight vector" << endl;
    return 0.0;
  }
  double sum = 0.0;
  for (int i = 0; i < v.size(); i++){
    if (magnitudes[i] > threshold){
      sum = sum + v[i]*weights[i];
    }
  }
  return sum;
}



double accumulateWeights_Threshold(vector<double> weights, vector<double> magnitudes, double threshold){
  double sum = 0.0;
  for (int i = 0; i < weights.size(); i++){
    if (magnitudes[i] > threshold){
      sum = sum + weights[i];
    }
  }
  return sum;
}

double getWeightedMean_Threshold(vector<double> v, vector<double> weights, vector<double> magnitudes, double threshold){
  if (v.size() != weights.size()){
    cerr << "ERROR: getWeightedMean(): value vector not same size as weight vector" << endl;
    return 0.0;
  }
  double sum = getWeightedSum_Threshold(v, weights, magnitudes, threshold);
  double sumOfWeights = accumulateWeights_Threshold(weights, magnitudes, threshold);
  //  double sumOfWeights = accumulate(weights.begin(), weights.end(), 0.0);
  double mean = sum / sumOfWeights;
  return mean;
}

void getWeightedVariance_Threshold(vector<double> v, vector<double> weights, vector<double> magnitudes, double threshold, double &sumOfSquares_out, double &sumOfWeights_out){
  vector<double> variance_weight;
 if (v.size() != weights.size()){
    cerr << "ERROR: getWeightedVariance(): value vector not same size as weight vector" << endl;
    return;
  }

  double mean = getWeightedMean_Threshold(v, weights, magnitudes, threshold);
  double sumOfSquares = 0.0;
  double sumOfWeights = 0.0;
  
  //  double sumOfWeights = accumulate(weights.begin(), weights.end(), 0.0);
  for (int i = 0; i < v.size(); i++){
    if (magnitudes[i] > threshold){
      sumOfSquares = sumOfSquares + weights[i]*pow((v[i] - mean), 2);
      sumOfWeights = sumOfWeights + weights[i];
    }
  }  

  sumOfSquares_out = sumOfSquares;
  sumOfWeights_out = sumOfWeights;

  return;
}

double getWeightedRmsPol_Threshold(int polarization, vector<vector<vector<double> > > delays, vector<vector<vector<double> > > weights, vector<vector<vector<double> > > magnitudes,  double threshold ){
  
  double totalVariance = 0.0;
  double totalWeight = 0.0;
  for (int i = 0; i < delays.size(); i++){
    for (int j = polarization; j < polarization+1; j++){
      //double sum = getWeightedSum(delays[i][j], weights[i][j]);
      bool isUsable = isUsablePattern(magnitudes[i][j], threshold);
      if (isUsable == true){
	double sumOfSquares;
	double sumOfWeights;
	getWeightedVariance_Threshold(delays[i][j], weights[i][j], magnitudes[i][j], threshold, sumOfSquares, sumOfWeights);
	totalVariance = totalVariance + sumOfSquares;
	totalWeight = totalWeight + sumOfWeights;
      }
    }
  }
  //  cout << totalRms << endl;
  if (totalWeight > 0.){
    return sqrt(totalVariance/totalWeight);
  }
  else {
    return 999.;
  }
}

double getLeastDifferenceFromTetrad(vector<vector<vector<double> > > delays)
{
  double leastDelayDifferenceFrom4 = 1.E9;
  for (int i = 0; i < delays.size(); i++){
    for (int j = 0; j < delays[i].size(); j++){
      for (int k = 0; k < delays[i][j].size(); k++){
	if (leastDelayDifferenceFrom4 > abs(delays[i][j][k])) {
	  leastDelayDifferenceFrom4 = abs(delays[i][j][k]);
	  
	}
      }
    }
  }
  return leastDelayDifferenceFrom4;
 
}

int getHitsUnderTimeFromTetrad(int polarization, double timeLimit, vector<vector<vector<double> > > delays)
{
  int hits = 0;
  for (int i = 0; i < delays.size(); i++){
    //    for (int j = 0; j < delays[i].size(); j++){
      for (int k = 0; k < delays[i][polarization].size(); k++){
	if (timeLimit > abs(delays[i][polarization][k]) ) {
	  //	  cout << timeLimit << delays[i][j][k] <<endl;
	  //	  cout << i << " : " << k << " : " << delays[i][polarization][k] << endl;
	  hits++;
	  
	}
      }
      //    }
  }
  return hits;
 
}

int getHitsUnderTimeFromTetrad(int polarization, double timeLimit, vector<vector<vector<double> > > delays, vector<vector<vector<double> > > weights)
{
  int hits = 0;
  for (int i = 0; i < delays.size(); i++){
    //    for (int j = 0; j < delays[i].size(); j++){
      for (int k = 0; k < delays[i][polarization].size(); k++){
	if (timeLimit > abs(delays[i][polarization][k]) ) {
	  //	  cout << timeLimit << delays[i][j][k] <<endl;
	  //	  cout << i << " : " << k << " : " << delays[i][polarization][k] << endl;
	  hits++;
	  
	}
      }
      //    }
  }
  return hits;
 
}


vector<vector<vector<double> > > getHitTimeDifferencesFromTetrad(vector<double> hitTimes, vector<double> peakIntPowers, vector<vector<vector<vector<int> > > > Pairs){
  
  vector<vector<vector<double> > > delays;

  delays.resize(Pairs.size());
  for (int i = 0; i < Pairs.size(); i++){
    delays[i].resize(Pairs[i].size());
    for (int j = 0; j < Pairs[i].size(); j++){
      //      delays[i][j].resize(Pairs[i][j].size());
      for (int k = 0; k < Pairs[i][j].size()-1; k++){
	for (int k2 = k+1; k2 < Pairs[i][j].size(); k2++){
	  double delay1 = hitTimes[Pairs[i][j][k][1]] - hitTimes[Pairs[i][j][k][0]];
	  double delay2 = hitTimes[Pairs[i][j][k2][1]] - hitTimes[Pairs[i][j][k2][0]];

	  if (j==0){
	    cout << Pairs[i][j][k][0] << " : " << Pairs[i][j][k][1]  << " : " << Pairs[i][j][k2][0] << " : " << Pairs[i][j][k2][1] << endl;
	  cout << i << " : " << j << " : " << k << " : " << k2 << " : " << delay1 << " : " << delay2 << " : " << delay2-delay1 << endl;
	  }
	  
	  delays[i][j].push_back(delay2-delay1);
	  //	cout << delay << endl;
	}
      }
    }
  }

  return delays;
}



/*
void getVAveragesDeviations_Divided(vector<TGraph*> vPlots, int divisions, vector<double> &vAverages, vector<double> &vDeviations){
  //vector<double> vDeviations;
  //vector<double> vAverages;
  vAverages.clear();
  vDeviations.clear();

  for (int i = 0; i< vPlots.size(); i++){
    double average;
    double deviation;
    getAverageDeviation_Divided(vPlots[i], divisions, average, deviation);
    vAverages.push_back(average);
    vDeviations.push_back(deviation);
  }
}
*/


vector<double> getVSigmas(vector<double> peaks, vector<double> averages, vector<double> variances){
  vector<double> vSigmas;

  for (int i = 0; i < peaks.size(); i++){
    double sigma = (peaks[i] - averages[i])/variances[i];
    vSigmas.push_back(sigma);
  }
  return vSigmas;

}

vector<vector<vector<double> > > getVVVSigmas(vector<double> vSigmas, vector<vector<vector<vector<int> > > > Pairs){
  
  vector<vector<vector<double> > > vvv;

  vvv.resize(Pairs.size());
  for (int i = 0; i < Pairs.size(); i++){
    vvv[i].resize(Pairs[i].size());
    for (int j = 0; j < Pairs[i].size(); j++){
      //      delays[i][j].resize(Pairs[i][j].size());
      for (int k = 0; k < Pairs[i][j].size(); k++){
	double weight = 0.;
	if ( vSigmas[Pairs[i][j][k][1]] > vSigmas[Pairs[i][j][k][0]]) {
	  weight = vSigmas[Pairs[i][j][k][0]];
	} else {
	  weight = vSigmas[Pairs[i][j][k][1]];
	}

	vvv[i][j].push_back(weight);
	
	//	cout << delay << endl;
      }
    }
  }

  return vvv;
}

double getLesser(double input1, double input2){
  if (input1 == input1 && input2 == input2){
    if (input1 > input2){
      return input2;
    } else {
      return input1;
    }
  } else {
    cerr << "getLesser: Comparison cannot be made! One of the arguments fails self-consistency!" << endl;
    return -1;
  }

}

void getMinimum(int length, double *array, int &index, double &minimum){ 
  double min = 1.E100;
  int minIndex = -1;
  for (int i = 0; i < length; i++){
    //    cout << array[i] << " : ";

    if (array[i] < min && array[i] != -1000){
      min = array[i];
      minIndex = i;
    }
  }
  //  cout << endl;
  index = minIndex;
  minimum = min;
}

double* getCorrelation_NoNorm(int length, double *oldY1, double *oldY2)
{

  //    cout << "Here in getCorrelation" << endl;                                                                                                                                    
  FFTWComplex *theFFT1=FFTtools::doFFT(length,oldY1);
  FFTWComplex *theFFT2=FFTtools::doFFT(length,oldY2);


  int newLength=(length/2)+1;
  //     cout << "newLength " << newLength << endl;                                                                                                                                  
  FFTWComplex *tempStep = new FFTWComplex [newLength];
  int no2=length>>1;
  for(int i=0;i<newLength;i++) {
    double reFFT1=theFFT1[i].re;
    double imFFT1=theFFT1[i].im;
    double reFFT2=theFFT2[i].re;
    double imFFT2=theFFT2[i].im;

    //Real part of output                                                                                                                                                      
    tempStep[i].re=(reFFT1*reFFT2+imFFT1*imFFT2);
    //Imaginary part of output                                                                                                                                                 
    tempStep[i].im=(imFFT1*reFFT2-reFFT1*imFFT2);
  }
  //    cout << "finished messing around" << endl;                                                                                                                                   
  double *theOutput=FFTtools::doInvFFT(length,tempStep);
  //    cout << "got inverse" << endl;                                                                                                                                               
  delete [] theFFT1;
  delete [] theFFT2;
  delete [] tempStep;
  return theOutput;

}

// test get Correlation graph with normalization factor from overlapping waveform's sqrt(power)                                                                                    
TGraph *getCorrelationGraph_WFweight(TGraph *gr1, TGraph *gr2) {
  //Now we'll extend this up to a power of 2                                                                                                                                      
  int length=gr1->GetN();
  int length2=gr2->GetN();

  //  cout << "Lengths: " << length << " : " << length2 << endl;

  int N=int(TMath::Power(2,int(TMath::Log2(length))+2));
  if(N<length2)
    N=int(TMath::Power(2,int(TMath::Log2(length2))+2));

  //Will really assume that N's are equal for now                                                                                                                                
  int firstRealSamp=(N-length)/2;

  double *oldY1 = new double [N];
  double *oldY2 = new double [N];

  double x,y;
  Double_t x2,y2;
  gr1->GetPoint(1,x2,y2);
  gr1->GetPoint(0,x,y);
  double deltaT=x2-x;
  double firstX=x;

  gr2->GetPoint(0,x2,y2);
  double waveOffset=firstX-x2;

  int OffsetBin = (int)(waveOffset/deltaT);


  //    gr1->GetPoint(N/2,x2,y2);                                                                                                                                                
  //    double offset=x-x2;                                                                                                                                                      
  //    std::cout << length << "\t" << length2 << "\n";                                                                                                                          

  for(int i=0;i<N;i++) {

    if(i<firstRealSamp || i>=firstRealSamp+length)
      y=0;
    else {
      gr1->GetPoint(i-firstRealSamp,x,y);
    }
    oldY1[i]=y;

    if(i<firstRealSamp || i>=firstRealSamp+length2)
      y=0;
    else {
      gr2->GetPoint(i-firstRealSamp,x,y);
    }
    oldY2[i]=y;

  }


  //    offset+=waveOffset;                                                                                                                                                      

  double *xVals = new double [N];
  double *yVals = new double [N];

  int dBin;
  double Norm1 = 0.;
  double Norm2 = 0.;

  //double *corVals=FFTtools::getCorrelation(N,oldY1,oldY2);                                                                         

  double *corVals=getCorrelation_NoNorm(N,oldY1,oldY2);
  /*
  for (int i = 0; i < N; i++){
    cout << corVals[i] << endl;
  }
  */
  for(int i=0;i<N;i++) {

    Norm1 = 0.;
    Norm2 = 0.;

    if(i<N/2) {
      //Positive                                                                                                                                                               
      xVals[i+(N/2)]=(i*deltaT)+waveOffset;
      //yVals[i+(N/2)]=corVals[i];                                                                                                                                             
      dBin = i+OffsetBin;

      if (dBin<0) {
	for (int i=-(dBin); i<N; i++) {
	  Norm1 += oldY1[i]*oldY1[i];
	}
	for (int i=0; i<N+(dBin); i++) {
	  Norm2 += oldY2[i]*oldY2[i];
	}
      }
      else { // dBin >= 0                                                                                                                                                      
	for (int i=0; i<N-(dBin); i++) {
	  Norm1 += oldY1[i]*oldY1[i];
	}
	for (int i=(dBin); i<N; i++) {
	  Norm2 += oldY2[i]*oldY2[i];
	}
      }

      //cout<<"Norm1 : "<<Norm1<<", Norm2 : "<<Norm2<<endl;                                                                                                                    

      if ( Norm1>0. && Norm2>0. )
	yVals[i+(N/2)]=corVals[i] / ( sqrt(Norm1)*sqrt(Norm2) );
      else
	yVals[i+(N/2)]=corVals[i];
      
    }
    else {
      //Negative                                                                                                                                                               
      xVals[i-(N/2)]=((i-N)*deltaT)+waveOffset;
      
      dBin = i-N+OffsetBin;
      
      if (dBin<0) {
	for (int i=-(dBin); i<N; i++) {
	  Norm1 += oldY1[i]*oldY1[i];
	}
	for (int i=0; i<N+(dBin); i++) {
	  Norm2 += oldY2[i]*oldY2[i];
	}
      }
      else { // dBin >= 0                                                                                                                                                      
	for (int i=0; i<N-(dBin); i++) {
	  Norm1 += oldY1[i]*oldY1[i];
	}
	for (int i=(dBin); i<N; i++) {
	  Norm2 += oldY2[i]*oldY2[i];
	}
      }

      //cout<<"Norm1 : "<<Norm1<<", Norm2 : "<<Norm2<<endl;                                                                                                                    

      if ( Norm1>0. && Norm2>0. )
	yVals[i-(N/2)]=corVals[i] / ( sqrt(Norm1)*sqrt(Norm2) );
      else
	yVals[i-(N/2)]=corVals[i];


    }
    //    cout << i<<"/"<<N << " : " << yVals[i-N/2] << endl;
  }

  
  TGraph *grCor = new TGraph(N,xVals,yVals);
  
  /*
  double x_test,y_test;
  for (int i = 0; i < N; i++){
    grCor->GetPoint(i, x_test, y_test);
    cout << i << "/" << N << " : " << x_test << " : " << y_test << endl;
  }
  */
  delete [] oldY1;
  delete [] oldY2;
  delete [] xVals;
  delete [] yVals;
  delete [] corVals;



  return grCor;
}

vector<TGraph*> getCorrelationGraphs_wBestTimes(vector<TGraph*> grIn, vector<vector< int > > pairs, vector<double> &bestTimes, vector<double> &bestCorrs){
  
  vector<TGraph*> grCorrs;
  double x_max, y_max;
  
  for (int pair = 0; pair < pairs.size(); pair++){
    //    cout << pair << " : " << pairs[pair][0] << " : " << pairs[pair][1] <<endl;
    TGraph* grCorr = getCorrelationGraph_WFweight(grIn[pairs[pair][0]], grIn[pairs[pair][1]]);
    //    TGraph* grCorr = FFTtools::getCorrelationGraph(grIn[pairs[pair][0]], grIn[pairs[pair][1]]);
    
    getMaximum(grCorr, x_max, y_max);
    //    cout << x_max << " : " << y_max << endl;
    bestTimes.push_back(x_max);
    bestCorrs.push_back(y_max);
    grCorrs.push_back(grCorr);
  }

  return grCorrs;
}

vector<TGraph*> getCorrelationGraphs(vector<TGraph*> grIn, vector<vector< int > > pairs ){
  
  vector<TGraph*> grCorrs;
  double x_max, y_max;
  
  for (int pair = 0; pair < pairs.size(); pair++){
    //    cout << pair << " : " << pairs[pair][0] << " : " << pairs[pair][1] <<endl;
    TGraph* grCorr = getCorrelationGraph_WFweight(grIn[pairs[pair][0]], grIn[pairs[pair][1]]);
    //    TGraph* grCorr = FFTtools::getCorrelationGraph(grIn[pairs[pair][0]], grIn[pairs[pair][1]]);
    
    //    getMaximum(grCorr, x_max, y_max);
    //    cout << x_max << " : " << y_max << endl;
    //    bestTimes.push_back(x_max);
    //    bestCorrs.push_back(y_max);
    grCorrs.push_back(grCorr);
  }

  return grCorrs;
}


vector<TGraph*> getHilbertGraphs_wBestTimes(vector<TGraph*> grIn, vector<double> &bestTimes, vector<double> &bestCorrs){

  vector<TGraph*> grHilberts;
  double x_max, y_max;

  for (int i = 0; i < grIn.size(); i++){
    TGraph* grHilbert = FFTtools::getHilbertEnvelope(grIn[i]);
    getMaximum(grHilbert, x_max, y_max);
    bestTimes.push_back(x_max);
    bestCorrs.push_back(y_max);
    grHilberts.push_back(grHilbert);
  }

  return grHilberts;

}

vector<TGraph*> getHilbertCorrelationGraphs_wBestTimes(vector<TGraph*> grIn, vector<vector< int > > pairs, vector<double> &bestTimes, vector<double> &bestCorrs){

  vector<TGraph*> grHilberts;
  double x_max, y_max;

  for (int pair = 0; pair < pairs.size(); pair++){
    //    cout << pair << " : " << pairs[pair][0] << " : " << pairs[pair][1] <<endl;
    TGraph* grCorr = getCorrelationGraph_WFweight(grIn[pairs[pair][0]], grIn[pairs[pair][1]]);
    TGraph* grHilbert = FFTtools::getHilbertEnvelope(grCorr);

    getMaximum(grHilbert, x_max, y_max);
    bestTimes.push_back(x_max);
    bestCorrs.push_back(y_max);
    grHilberts.push_back(grHilbert);
    delete grCorr;
  }

  return grHilberts;
}



void setupCorrelationPairs(int StationID, vector<vector<int> > &pairs_V, vector<vector<int> > &pairs_H){

  int chan_V[] = {0,1,2,3,4,5,6,7};
  int chan_H[] = {8,9,10,11,12,13,14,15};

  for (int i = 0; i < 8; i++){

    for (int j = i+1; j<8; j++){
      vector<int> pair;
      pair.push_back(chan_V[i]);
      pair.push_back(chan_V[j]);
      pairs_V.push_back(pair);
    }

    for (int j = i+1; j<8; j++){
      vector<int> pair;
      pair.push_back(chan_H[i]);
      pair.push_back(chan_H[j]);
      pairs_H.push_back(pair);
    }
  }
}


vector<double> getRms_Corr(vector<double> BestTimes, vector<vector<int> > pairs_from_correlation, vector<vector<vector<vector<int> > > > pairs_from_faces, vector<vector<double> > ant_loc){

  //  vector< vector<double> > DelaysForFaces;
  vector< vector<double> > DelaysForFaces_Norm;

  for (int polarization = 0; polarization < 2; polarization++){
  for (int faceType = 0; faceType < pairs_from_faces.size(); faceType++){
    vector<double> delays_temp_norm;
    for (int i = 0; i < pairs_from_correlation.size(); i++){
      for (int j = 0; j < pairs_from_faces[faceType][polarization].size(); j++){
	if (pairs_from_correlation[i][0] == pairs_from_faces[faceType][polarization][j][0] && pairs_from_correlation[i][1] == pairs_from_faces[faceType][polarization][j][1]){
	  double baseline_length = sqrt( pow( (ant_loc[pairs_from_correlation[i][0]][0] - ant_loc[pairs_from_correlation[i][1]][0]),2)+
                                         pow( (ant_loc[pairs_from_correlation[i][0]][1] - ant_loc[pairs_from_correlation[i][1]][1]),2)+
                                         pow( (ant_loc[pairs_from_correlation[i][0]][2] - ant_loc[pairs_from_correlation[i][1]][2]),2) );
	  delays_temp_norm.push_back(BestTimes[i]/baseline_length*c*nsTos/n_ice);

	  //	  delays_temp.push_back(bestTimes_V_Hilbert[i]);
	}
	if (pairs_from_correlation[i][1] == pairs_from_faces[faceType][polarization][j][0] && pairs_from_correlation[i][0] == pairs_from_faces[faceType][polarization][j][1]){
	  //	  delays_temp.push_back(bestTimes_V_Hilbert[i]);
	  double baseline_length = sqrt( pow( (ant_loc[pairs_from_correlation[i][0]][0] - ant_loc[pairs_from_correlation[i][1]][0]),2)+
                                         pow( (ant_loc[pairs_from_correlation[i][0]][1] - ant_loc[pairs_from_correlation[i][1]][1]),2)+
                                         pow( (ant_loc[pairs_from_correlation[i][0]][2] - ant_loc[pairs_from_correlation[i][1]][2]),2) );
	  delays_temp_norm.push_back(BestTimes[i]/baseline_length*c*nsTos/n_ice);

	}
      }
    }
    //    DelaysForFaces.push_back(delays_temp);
    DelaysForFaces_Norm.push_back(delays_temp_norm);
  }
  }

  vector<double> RmsValues;
 
  for (int faceType = 0; faceType < DelaysForFaces_Norm.size(); faceType++){
    RmsValues.push_back(getRms(DelaysForFaces_Norm[faceType]));
    //    cout << event << " : " << RmsValues[faceType] << " : " << RmsValues_Norm[faceType] << endl;
  }

  return RmsValues;
}

vector<double> getRms_Corr_Threshold(vector<double> BestTimes, vector<vector<int> > pairs_from_correlation, vector<vector<vector<vector<int> > > > pairs_from_faces, vector<vector<double> > ant_loc, vector<double> BestCorrs, double CorrThreshold){

  //  vector< vector<double> > DelaysForFaces;
  vector< vector<double> > DelaysForFaces_Norm;

  double baseline_length = 0.;

  for (int polarization = 0; polarization < 2; polarization++){  
  for (int faceType = 0; faceType < pairs_from_faces.size(); faceType++){
    vector<double> delays_temp_norm;
    for (int i = 0; i < pairs_from_correlation.size(); i++){
      for (int j = 0; j < pairs_from_faces[faceType][polarization].size(); j++){
	if (BestCorrs[i] > CorrThreshold){
	  if (pairs_from_correlation[i][0] == pairs_from_faces[faceType][polarization][j][0] && pairs_from_correlation[i][1] == pairs_from_faces[faceType][polarization][j][1]){
	    baseline_length = sqrt( pow( (ant_loc[pairs_from_correlation[i][0]][0] - ant_loc[pairs_from_correlation[i][1]][0]),2)+
				    pow( (ant_loc[pairs_from_correlation[i][0]][1] - ant_loc[pairs_from_correlation[i][1]][1]),2)+
				    pow( (ant_loc[pairs_from_correlation[i][0]][2] - ant_loc[pairs_from_correlation[i][1]][2]),2) );
	    delays_temp_norm.push_back(BestTimes[i]/baseline_length*c*nsTos/n_ice);
	  }
	  //	  delays_temp.push_back(bestTimes_V_Hilbert[i]);
	  
	  if (pairs_from_correlation[i][1] == pairs_from_faces[faceType][polarization][j][0] && pairs_from_correlation[i][0] == pairs_from_faces[faceType][polarization][j][1]){
	    //	  delays_temp.push_back(bestTimes_V_Hilbert[i]);
	    baseline_length = sqrt( pow( (ant_loc[pairs_from_correlation[i][0]][0] - ant_loc[pairs_from_correlation[i][1]][0]),2)+
				    pow( (ant_loc[pairs_from_correlation[i][0]][1] - ant_loc[pairs_from_correlation[i][1]][1]),2)+
				    pow( (ant_loc[pairs_from_correlation[i][0]][2] - ant_loc[pairs_from_correlation[i][1]][2]),2) );
	    delays_temp_norm.push_back(BestTimes[i]/baseline_length*c*nsTos/n_ice);
	    
	  }
	}
      }
    }
  //    DelaysForFaces.push_back(delays_temp);
  DelaysForFaces_Norm.push_back(delays_temp_norm);
  }
  }

  vector<double> RmsValues;
 
 
  for (int faceType = 0; faceType < DelaysForFaces_Norm.size(); faceType++){
    if (DelaysForFaces_Norm[faceType].size() > 1){
      RmsValues.push_back(getRms(DelaysForFaces_Norm[faceType]));
    //    cout << event << " : " << RmsValues[faceType] << " : " << RmsValues_Norm[faceType] << endl;
    } else {
      RmsValues.push_back(2000.);
    }

  }
  
  return RmsValues;
}



double  getPolarizationRatio(vector<TGraph*> waveforms, vector<int> polarizations){

  if (waveforms.size() != polarizations.size()){
    cerr << "waveforms vector not the same size as the polarizations vector!" << endl;
    return -100.;
  }
  double totalPower[2] = {0., 0.};
  int totalSamples[2] = {0, 0};

  for (int channel = 0; channel < waveforms.size(); channel++){
    double power = integrateAllPower(waveforms[channel]);
    int samples = waveforms[channel]->GetN();
    if (polarizations[channel] == 0){
      totalPower[0] = totalPower[0]+power;
      totalSamples[0] = totalSamples[0]+samples;
    }
    if (polarizations[channel] == 1){
      totalPower[1] = totalPower[1]+power;
      totalSamples[1] = totalSamples[1]+samples;
    }
  }
  double averagePower[2];  
  averagePower[0] = totalPower[0]/(double)totalSamples[0];
  averagePower[1] = totalPower[1]/(double)totalSamples[1];
  
  double ratio = averagePower[0]/averagePower[1];

  return ratio;
}

void getThirdVPeakOverRMS(vector<double> vPeakOverRms, vector<int> polarizations, vector<double> &ThirdVpeakOverRms){

  ThirdVpeakOverRms.resize(3);

  vector<double> vPeakOverRms_V;
  vector<double> vPeakOverRms_H;
  vector<double> vPeakOverRms_Total;

  for (int chan = 0; chan < vPeakOverRms.size(); chan++){
    if (polarizations[chan] == 0){
      vPeakOverRms_V.push_back(vPeakOverRms[chan]);
    }
    if (polarizations[chan] == 1){
      vPeakOverRms_H.push_back(vPeakOverRms[chan]);
    }
    vPeakOverRms_Total.push_back(vPeakOverRms[chan]);
  }
  
  sort(vPeakOverRms_V.begin(), vPeakOverRms_V.end());
  sort(vPeakOverRms_H.begin(), vPeakOverRms_H.end());
  sort(vPeakOverRms_Total.begin(), vPeakOverRms_Total.end());
  
  ThirdVpeakOverRms[0] = vPeakOverRms_V[(int)vPeakOverRms_V.size()-3];
  ThirdVpeakOverRms[1] = vPeakOverRms_H[(int)vPeakOverRms_H.size()-3];
  ThirdVpeakOverRms[2] = vPeakOverRms_Total[(int)vPeakOverRms_Total.size()-3];
}

double getRms_Remove1(vector<double> v){
  double lowestRMS = 1.0E5;
  /*
    for(int i = 0; i < v.size(); i++){
      cout << v[i] << " : ";
    }
    cout << endl;
  */
  for (int removed = 0; removed < v.size(); removed++){
    vector<double> v_subset;
    for(int i = 0; i < v.size(); i++){
      if (i != removed){
	v_subset.push_back(v[i]);
      } 
    }
    double mean = getMean(v_subset);
    double sumOfSquares = 0.0;
    for (int i = 0; i < v_subset.size(); i++){
      //  cout << v_subset[i] << " : ";
      sumOfSquares = sumOfSquares + pow(v_subset[i] - mean, 2);
    }
    //    cout << endl;
    double rms = sqrt(sumOfSquares / (double)v_subset.size());
    if (rms < lowestRMS){
      lowestRMS = rms;
    }
  }


    return lowestRMS;
}

double getRms_Remove1(vector<double> v, int &dropped_pair){
  double lowestRMS = 1.0E5;
  /*
    for(int i = 0; i < v.size(); i++){
      cout << v[i] << " : ";
    }
    cout << endl;
  */
  for (int removed = 0; removed < v.size(); removed++){
    vector<double> v_subset;
    for(int i = 0; i < v.size(); i++){
      if (i != removed){
	v_subset.push_back(v[i]);
      } 
    }
    double mean = getMean(v_subset);
    double sumOfSquares = 0.0;
    for (int i = 0; i < v_subset.size(); i++){
      //  cout << v_subset[i] << " : ";
      sumOfSquares = sumOfSquares + pow(v_subset[i] - mean, 2);
    }
    //    cout << endl;
    double rms = sqrt(sumOfSquares / (double)v_subset.size());
    if (rms < lowestRMS){
      lowestRMS = rms;
      dropped_pair = removed;
    }
  }


    return lowestRMS;
}


vector<double> getRms_Corr_Regular_Threshold_Remove1(vector<double> BestTimes, vector<vector<int> > pairs_from_correlation, vector<vector<vector<vector<int> > > > pairs_from_faces, vector<vector<double> > ant_loc, vector<double> BestCorrs, double CorrThreshold){

  //  vector< vector<double> > DelaysForFaces;
  vector< vector<double> > DelaysForFaces_Norm;

  double baseline_length = 0.;
  for (int polarization = 0; polarization < 2; polarization++){
  for (int faceType = 0; faceType < pairs_from_faces.size(); faceType++){
    vector<double> delays_temp_norm;
    for (int i = 0; i < pairs_from_correlation.size(); i++){
      for (int j = 0; j < pairs_from_faces[faceType][polarization].size(); j++){
	if (BestCorrs[i] > CorrThreshold){
	  if (pairs_from_correlation[i][0] == pairs_from_faces[faceType][polarization][j][0] && pairs_from_correlation[i][1] == pairs_from_faces[faceType][polarization][j][1]){
	    baseline_length = sqrt( pow( (ant_loc[pairs_from_correlation[i][0]][0] - ant_loc[pairs_from_correlation[i][1]][0]),2)+
				    pow( (ant_loc[pairs_from_correlation[i][0]][1] - ant_loc[pairs_from_correlation[i][1]][1]),2)+
				    pow( (ant_loc[pairs_from_correlation[i][0]][2] - ant_loc[pairs_from_correlation[i][1]][2]),2) );
	    delays_temp_norm.push_back(BestTimes[i]/baseline_length*c*nsTos/n_ice);
	  }
	  //	  delays_temp.push_back(bestTimes_V_Hilbert[i]);
	  
	  if (pairs_from_correlation[i][1] == pairs_from_faces[faceType][polarization][j][0] && pairs_from_correlation[i][0] == pairs_from_faces[faceType][polarization][j][1]){
	    //	  delays_temp.push_back(bestTimes_V_Hilbert[i]);
	    baseline_length = sqrt( pow( (ant_loc[pairs_from_correlation[i][0]][0] - ant_loc[pairs_from_correlation[i][1]][0]),2)+
				    pow( (ant_loc[pairs_from_correlation[i][0]][1] - ant_loc[pairs_from_correlation[i][1]][1]),2)+
				    pow( (ant_loc[pairs_from_correlation[i][0]][2] - ant_loc[pairs_from_correlation[i][1]][2]),2) );
	    delays_temp_norm.push_back(BestTimes[i]/baseline_length*c*nsTos/n_ice);
	    
	  }
	}
      }
    }
  //    DelaysForFaces.push_back(delays_temp);
  DelaysForFaces_Norm.push_back(delays_temp_norm);
  }
  }
  vector<double> RmsValues;
 
 
  for (int faceType = 0; faceType < DelaysForFaces_Norm.size(); faceType++){
    if (DelaysForFaces_Norm[faceType].size() > 2){
      RmsValues.push_back(getRms_Remove1(DelaysForFaces_Norm[faceType]));
    //    cout << event << " : " << RmsValues[faceType] << " : " << RmsValues_Norm[faceType] << endl;
    } else {
      RmsValues.push_back(2000.);
    }

  }
  
  return RmsValues;
}

vector<double> getRms_Corr_Regular_Threshold_Remove1(vector<double> BestTimes, vector<vector<int> > pairs_from_correlation, vector<vector<vector<vector<int> > > > pairs_from_faces, vector<vector<double> > ant_loc, vector<double> BestCorrs, double CorrThreshold, vector<int> &dropped_pairs){

  dropped_pairs.clear();
  //  vector< vector<double> > DelaysForFaces;
  vector< vector<double> > DelaysForFaces_Norm;
  

  double baseline_length = 0.;
  for (int polarization = 0; polarization < 2; polarization++){
  for (int faceType = 0; faceType < pairs_from_faces.size(); faceType++){
    vector<double> delays_temp_norm;
    for (int i = 0; i < pairs_from_correlation.size(); i++){
      for (int j = 0; j < pairs_from_faces[faceType][polarization].size(); j++){
	if (BestCorrs[i] > CorrThreshold){
	  if (pairs_from_correlation[i][0] == pairs_from_faces[faceType][polarization][j][0] && pairs_from_correlation[i][1] == pairs_from_faces[faceType][polarization][j][1]){
	    baseline_length = sqrt( pow( (ant_loc[pairs_from_correlation[i][0]][0] - ant_loc[pairs_from_correlation[i][1]][0]),2)+
				    pow( (ant_loc[pairs_from_correlation[i][0]][1] - ant_loc[pairs_from_correlation[i][1]][1]),2)+
				    pow( (ant_loc[pairs_from_correlation[i][0]][2] - ant_loc[pairs_from_correlation[i][1]][2]),2) );
	    delays_temp_norm.push_back(BestTimes[i]/baseline_length*c*nsTos/n_ice);
	  }
	  //	  delays_temp.push_back(bestTimes_V_Hilbert[i]);
	  
	  if (pairs_from_correlation[i][1] == pairs_from_faces[faceType][polarization][j][0] && pairs_from_correlation[i][0] == pairs_from_faces[faceType][polarization][j][1]){
	    //	  delays_temp.push_back(bestTimes_V_Hilbert[i]);
	    baseline_length = sqrt( pow( (ant_loc[pairs_from_correlation[i][0]][0] - ant_loc[pairs_from_correlation[i][1]][0]),2)+
				    pow( (ant_loc[pairs_from_correlation[i][0]][1] - ant_loc[pairs_from_correlation[i][1]][1]),2)+
				    pow( (ant_loc[pairs_from_correlation[i][0]][2] - ant_loc[pairs_from_correlation[i][1]][2]),2) );
	    delays_temp_norm.push_back(BestTimes[i]/baseline_length*c*nsTos/n_ice);
	    
	  }
	}
      }
    }
  //    DelaysForFaces.push_back(delays_temp);
  DelaysForFaces_Norm.push_back(delays_temp_norm);
  }
  }
  vector<double> RmsValues;
 
 
  for (int faceType = 0; faceType < DelaysForFaces_Norm.size(); faceType++){
    int dropped_pair;
    if (DelaysForFaces_Norm[faceType].size() > 2){
      RmsValues.push_back(getRms_Remove1(DelaysForFaces_Norm[faceType], dropped_pair));
    //    cout << event << " : " << RmsValues[faceType] << " : " << RmsValues_Norm[faceType] << endl;
    } else {
      RmsValues.push_back(2000.);
      dropped_pair = -1;
    }
    dropped_pairs.push_back(dropped_pair);
  }
  
  return RmsValues;
}

vector<vector<vector<vector<int> > > > setupFaces(int stationId){

  vector<vector<vector<vector<int> > > > Pairs; //Pairs_pol_face_pair_ind;
  vector<vector<vector<int> > > faces;
  vector<vector<int> > face;
  vector<int> pair;

  if (stationId == 2){

//VPol
    //Face 1 - side face
    pair.push_back(0); pair.push_back(4);
    face.push_back(pair); pair.clear();
    pair.push_back(1); pair.push_back(5);
    face.push_back(pair); pair.clear();
    faces.push_back(face);face.clear();

    //Face 2 - side face
    pair.push_back(1); pair.push_back(5);
    face.push_back(pair); pair.clear();
    pair.push_back(3); pair.push_back(7);
    face.push_back(pair); pair.clear();
    faces.push_back(face);face.clear();

    //Face 3 - side face
    pair.push_back(3); pair.push_back(7);
    face.push_back(pair); pair.clear();
    pair.push_back(2); pair.push_back(6);
    face.push_back(pair); pair.clear();
    faces.push_back(face);face.clear();

    //Face 4 - side face
    pair.push_back(2); pair.push_back(6);
    face.push_back(pair); pair.clear();
    pair.push_back(0); pair.push_back(4);
    face.push_back(pair); pair.clear();
    faces.push_back(face);face.clear();

    //Face 5 - top face
    pair.push_back(0); pair.push_back(1);
    face.push_back(pair); pair.clear();
    pair.push_back(2); pair.push_back(3);
    face.push_back(pair); pair.clear();
    faces.push_back(face);face.clear();

    //Face 6 - bottom face
    pair.push_back(4); pair.push_back(5);
    face.push_back(pair); pair.clear();
    pair.push_back(6); pair.push_back(7);
    face.push_back(pair); pair.clear();
    faces.push_back(face);face.clear();

    //Face 7 - cross face
    pair.push_back(0); pair.push_back(4);
    face.push_back(pair); pair.clear();
    pair.push_back(3); pair.push_back(7);
    face.push_back(pair); pair.clear();
    faces.push_back(face);face.clear();

    //Face 8 - cross face
    pair.push_back(1); pair.push_back(5);
    face.push_back(pair); pair.clear();
    pair.push_back(2); pair.push_back(6);
    face.push_back(pair); pair.clear();
    faces.push_back(face);face.clear();

    //Face 9 - slant face
    pair.push_back(0); pair.push_back(5);
    face.push_back(pair); pair.clear();
    pair.push_back(2); pair.push_back(7);
    face.push_back(pair); pair.clear();
    faces.push_back(face);face.clear();

    //Face 10 - slant face
    pair.push_back(4); pair.push_back(1);
    face.push_back(pair); pair.clear();
    pair.push_back(6); pair.push_back(3);
    face.push_back(pair); pair.clear();
    faces.push_back(face);face.clear();

    //Face 11 - slant face
    pair.push_back(0); pair.push_back(6);
    face.push_back(pair); pair.clear();
    pair.push_back(1); pair.push_back(7);
    face.push_back(pair); pair.clear();
    faces.push_back(face);face.clear();

    //Face 12 - slant face
    pair.push_back(2); pair.push_back(4);
    face.push_back(pair); pair.clear();
    pair.push_back(3); pair.push_back(5);
    face.push_back(pair); pair.clear();
    faces.push_back(face);face.clear();

    Pairs.push_back(faces);faces.clear();
    
//HPol
    //Face 1 - side face
    pair.push_back(8); pair.push_back(12);
    face.push_back(pair); pair.clear();
    pair.push_back(9); pair.push_back(13);
    face.push_back(pair); pair.clear();
    faces.push_back(face);face.clear();

    //Face 2 - side face
    pair.push_back(9); pair.push_back(13);
    face.push_back(pair); pair.clear();
    pair.push_back(11); pair.push_back(15);
    face.push_back(pair); pair.clear();
    faces.push_back(face);face.clear();

    //Face 3 - side face
    pair.push_back(11); pair.push_back(15);
    face.push_back(pair); pair.clear();
    pair.push_back(10); pair.push_back(14);
    face.push_back(pair); pair.clear();
    faces.push_back(face);face.clear();

    //Face 4 - side face
    pair.push_back(10); pair.push_back(14);
    face.push_back(pair); pair.clear();
    pair.push_back(8); pair.push_back(12);
    face.push_back(pair); pair.clear();
    faces.push_back(face);face.clear();

    //Face 5 - top face
    pair.push_back(8); pair.push_back(9);
    face.push_back(pair); pair.clear();
    pair.push_back(10); pair.push_back(11);
    face.push_back(pair); pair.clear();
    faces.push_back(face);face.clear();

    //Face 6 - bottom face
    pair.push_back(12); pair.push_back(13);
    face.push_back(pair); pair.clear();
    pair.push_back(14); pair.push_back(15);
    face.push_back(pair); pair.clear();
    faces.push_back(face);face.clear();

    //Face 7 - cross face
    pair.push_back(8); pair.push_back(12);
    face.push_back(pair); pair.clear();
    pair.push_back(11); pair.push_back(15);
    face.push_back(pair); pair.clear();
    faces.push_back(face);face.clear();

    //Face 8 - cross face
    pair.push_back(9); pair.push_back(13);
    face.push_back(pair); pair.clear();
    pair.push_back(10); pair.push_back(14);
    face.push_back(pair); pair.clear();
    faces.push_back(face);face.clear();

    //Face 9 - slant face
    pair.push_back(8); pair.push_back(13);
    face.push_back(pair); pair.clear();
    pair.push_back(10); pair.push_back(15);
    face.push_back(pair); pair.clear();
    faces.push_back(face);face.clear();

    //Face 10 - slant face
    pair.push_back(12); pair.push_back(9);
    face.push_back(pair); pair.clear();
    pair.push_back(14); pair.push_back(11);
    face.push_back(pair); pair.clear();
    faces.push_back(face);face.clear();

    //Face 11 - slant face
    pair.push_back(8); pair.push_back(14);
    face.push_back(pair); pair.clear();
    pair.push_back(9); pair.push_back(15);
    face.push_back(pair); pair.clear();
    faces.push_back(face);face.clear();

    //Face 12 - slant face
    pair.push_back(10); pair.push_back(12);
    face.push_back(pair); pair.clear();
    pair.push_back(11); pair.push_back(13);
    face.push_back(pair); pair.clear();
    faces.push_back(face);face.clear();
    Pairs.push_back(faces);faces.clear();

  }

  if (stationId == 3){

//VPol
    //Face 1 - side face
    pair.push_back(1); pair.push_back(5);
    face.push_back(pair); pair.clear();
    pair.push_back(0); pair.push_back(4);
    face.push_back(pair); pair.clear();
    faces.push_back(face);face.clear();

    //Face 2 - side face
    pair.push_back(0); pair.push_back(4);
    face.push_back(pair); pair.clear();
    pair.push_back(3); pair.push_back(7);
    face.push_back(pair); pair.clear();
    faces.push_back(face);face.clear();

    //Face 3 - side face
    pair.push_back(3); pair.push_back(7);
    face.push_back(pair); pair.clear();
    pair.push_back(2); pair.push_back(6);
    face.push_back(pair); pair.clear();
    faces.push_back(face);face.clear();

    //Face 4 - side face
    pair.push_back(2); pair.push_back(6);
    face.push_back(pair); pair.clear();
    pair.push_back(1); pair.push_back(5);
    face.push_back(pair); pair.clear();
    faces.push_back(face);face.clear();

    //Face 5 - top face
    pair.push_back(1); pair.push_back(0);
    face.push_back(pair); pair.clear();
    pair.push_back(2); pair.push_back(3);
    face.push_back(pair); pair.clear();
    faces.push_back(face);face.clear();

    //Face 6 - bottom face
    pair.push_back(5); pair.push_back(4);
    face.push_back(pair); pair.clear();
    pair.push_back(6); pair.push_back(7);
    face.push_back(pair); pair.clear();
    faces.push_back(face);face.clear();

    //Face 7 - cross face
    pair.push_back(1); pair.push_back(5);
    face.push_back(pair); pair.clear();
    pair.push_back(3); pair.push_back(7);
    face.push_back(pair); pair.clear();
    faces.push_back(face);face.clear();

    //Face 8 - cross face
    pair.push_back(0); pair.push_back(4);
    face.push_back(pair); pair.clear();
    pair.push_back(2); pair.push_back(6);
    face.push_back(pair); pair.clear();
    faces.push_back(face);face.clear();

    //Face 9 - slant face
    pair.push_back(1); pair.push_back(4);
    face.push_back(pair); pair.clear();
    pair.push_back(2); pair.push_back(7);
    face.push_back(pair); pair.clear();
    faces.push_back(face);face.clear();

    //Face 10 - slant face
    pair.push_back(5); pair.push_back(0);
    face.push_back(pair); pair.clear();
    pair.push_back(6); pair.push_back(3);
    face.push_back(pair); pair.clear();
    faces.push_back(face);face.clear();

    //Face 11 - slant face
    pair.push_back(1); pair.push_back(6);
    face.push_back(pair); pair.clear();
    pair.push_back(0); pair.push_back(7);
    face.push_back(pair); pair.clear();
    faces.push_back(face);face.clear();

    //Face 12 - slant face
    pair.push_back(2); pair.push_back(5);
    face.push_back(pair); pair.clear();
    pair.push_back(3); pair.push_back(4);
    face.push_back(pair); pair.clear();
    faces.push_back(face);face.clear();

    Pairs.push_back(faces);faces.clear();
    
//HPol
    //Face 1 - side face
    pair.push_back(9); pair.push_back(13);
    face.push_back(pair); pair.clear();
    pair.push_back(8); pair.push_back(12);
    face.push_back(pair); pair.clear();
    faces.push_back(face);face.clear();

    //Face 2 - side face
    pair.push_back(8); pair.push_back(12);
    face.push_back(pair); pair.clear();
    pair.push_back(11); pair.push_back(15);
    face.push_back(pair); pair.clear();
    faces.push_back(face);face.clear();

    //Face 3 - side face
    pair.push_back(11); pair.push_back(15);
    face.push_back(pair); pair.clear();
    pair.push_back(10); pair.push_back(14);
    face.push_back(pair); pair.clear();
    faces.push_back(face);face.clear();

    //Face 4 - side face
    pair.push_back(10); pair.push_back(14);
    face.push_back(pair); pair.clear();
    pair.push_back(9); pair.push_back(13);
    face.push_back(pair); pair.clear();
    faces.push_back(face);face.clear();

    //Face 5 - top face
    pair.push_back(9); pair.push_back(8);
    face.push_back(pair); pair.clear();
    pair.push_back(10); pair.push_back(11);
    face.push_back(pair); pair.clear();
    faces.push_back(face);face.clear();

    //Face 6 - bottom face
    pair.push_back(13); pair.push_back(12);
    face.push_back(pair); pair.clear();
    pair.push_back(14); pair.push_back(15);
    face.push_back(pair); pair.clear();
    faces.push_back(face);face.clear();

    //Face 7 - cross face
    pair.push_back(9); pair.push_back(13);
    face.push_back(pair); pair.clear();
    pair.push_back(11); pair.push_back(15);
    face.push_back(pair); pair.clear();
    faces.push_back(face);face.clear();

    //Face 8 - cross face
    pair.push_back(8); pair.push_back(12);
    face.push_back(pair); pair.clear();
    pair.push_back(10); pair.push_back(14);
    face.push_back(pair); pair.clear();
    faces.push_back(face);face.clear();

    //Face 9 - slant face
    pair.push_back(9); pair.push_back(12);
    face.push_back(pair); pair.clear();
    pair.push_back(10); pair.push_back(15);
    face.push_back(pair); pair.clear();
    faces.push_back(face);face.clear();

    //Face 10 - slant face
    pair.push_back(13); pair.push_back(8);
    face.push_back(pair); pair.clear();
    pair.push_back(14); pair.push_back(11);
    face.push_back(pair); pair.clear();
    faces.push_back(face);face.clear();

    //Face 11 - slant face
    pair.push_back(9); pair.push_back(14);
    face.push_back(pair); pair.clear();
    pair.push_back(8); pair.push_back(15);
    face.push_back(pair); pair.clear();
    faces.push_back(face);face.clear();

    //Face 12 - slant face
    pair.push_back(10); pair.push_back(13);
    face.push_back(pair); pair.clear();
    pair.push_back(11); pair.push_back(12);
    face.push_back(pair); pair.clear();
    faces.push_back(face);face.clear();
    Pairs.push_back(faces);faces.clear();

  }

  return Pairs;
}

vector<double> getRms_Faces(
			    vector<double> BestTimes,
			    int polarization,
			    vector<vector<int> > pairs_from_correlation, 
			    vector<vector<vector<vector<int> > > > pairs_pol_faces_pair_ind, 
			    vector<vector<double> > ant_loc){

  vector<double> RmsValues;
  
  double baseline_length = 0.;
  //  for (int polarization = 0; polarization < 2; polarization++){
    for (int faceType = 0; faceType < pairs_pol_faces_pair_ind[polarization].size(); faceType++){
      vector<double> delays_temp_norm_standard;
      vector<double> delays_temp_norm_cross;
      for (int i = 0; i < pairs_from_correlation.size(); i++){
	
	//standard pairs
	if (pairs_from_correlation[i][0] == pairs_pol_faces_pair_ind[polarization][faceType][0][0] && pairs_from_correlation[i][1] == pairs_pol_faces_pair_ind[polarization][faceType][0][1]){
	  baseline_length = sqrt( pow( (ant_loc[pairs_from_correlation[i][0]][0] - ant_loc[pairs_from_correlation[i][1]][0]),2)+
				  pow( (ant_loc[pairs_from_correlation[i][0]][1] - ant_loc[pairs_from_correlation[i][1]][1]),2)+
				  pow( (ant_loc[pairs_from_correlation[i][0]][2] - ant_loc[pairs_from_correlation[i][1]][2]),2) );
	  delays_temp_norm_standard.push_back(BestTimes[i]/baseline_length*c*nsTos/n_ice);
	}
	
	if (pairs_from_correlation[i][1] == pairs_pol_faces_pair_ind[polarization][faceType][0][0] && pairs_from_correlation[i][0] == pairs_pol_faces_pair_ind[polarization][faceType][0][1]){
	  
	  baseline_length = sqrt( pow( (ant_loc[pairs_from_correlation[i][0]][0] - ant_loc[pairs_from_correlation[i][1]][0]),2)+
				  pow( (ant_loc[pairs_from_correlation[i][0]][1] - ant_loc[pairs_from_correlation[i][1]][1]),2)+
				  pow( (ant_loc[pairs_from_correlation[i][0]][2] - ant_loc[pairs_from_correlation[i][1]][2]),2) );
	  delays_temp_norm_standard.push_back(BestTimes[i]/baseline_length*c*nsTos/n_ice * -1.);
	}

	if (pairs_from_correlation[i][0] == pairs_pol_faces_pair_ind[polarization][faceType][1][0] && pairs_from_correlation[i][1] == pairs_pol_faces_pair_ind[polarization][faceType][1][1]){
	  baseline_length = sqrt( pow( (ant_loc[pairs_from_correlation[i][0]][0] - ant_loc[pairs_from_correlation[i][1]][0]),2)+
				  pow( (ant_loc[pairs_from_correlation[i][0]][1] - ant_loc[pairs_from_correlation[i][1]][1]),2)+
				  pow( (ant_loc[pairs_from_correlation[i][0]][2] - ant_loc[pairs_from_correlation[i][1]][2]),2) );
	  delays_temp_norm_standard.push_back(BestTimes[i]/baseline_length*c*nsTos/n_ice);
	}
	
	if (pairs_from_correlation[i][1] == pairs_pol_faces_pair_ind[polarization][faceType][1][0] && pairs_from_correlation[i][0] == pairs_pol_faces_pair_ind[polarization][faceType][1][1]){
	  
	  baseline_length = sqrt( pow( (ant_loc[pairs_from_correlation[i][0]][0] - ant_loc[pairs_from_correlation[i][1]][0]),2)+
				  pow( (ant_loc[pairs_from_correlation[i][0]][1] - ant_loc[pairs_from_correlation[i][1]][1]),2)+
				  pow( (ant_loc[pairs_from_correlation[i][0]][2] - ant_loc[pairs_from_correlation[i][1]][2]),2) );
	  delays_temp_norm_standard.push_back(BestTimes[i]/baseline_length*c*nsTos/n_ice * -1.);
	}
	

	//cross pairs
	if (pairs_from_correlation[i][0] == pairs_pol_faces_pair_ind[polarization][faceType][0][0] && pairs_from_correlation[i][1] == pairs_pol_faces_pair_ind[polarization][faceType][1][0]){
	  baseline_length = sqrt( pow( (ant_loc[pairs_from_correlation[i][0]][0] - ant_loc[pairs_from_correlation[i][1]][0]),2)+
				  pow( (ant_loc[pairs_from_correlation[i][0]][1] - ant_loc[pairs_from_correlation[i][1]][1]),2)+
				  pow( (ant_loc[pairs_from_correlation[i][0]][2] - ant_loc[pairs_from_correlation[i][1]][2]),2) );
	  delays_temp_norm_cross.push_back(BestTimes[i]/baseline_length*c*nsTos/n_ice);
	}
	
	if (pairs_from_correlation[i][1] == pairs_pol_faces_pair_ind[polarization][faceType][0][0] && pairs_from_correlation[i][0] == pairs_pol_faces_pair_ind[polarization][faceType][1][0]){
	  
	  baseline_length = sqrt( pow( (ant_loc[pairs_from_correlation[i][0]][0] - ant_loc[pairs_from_correlation[i][1]][0]),2)+
				  pow( (ant_loc[pairs_from_correlation[i][0]][1] - ant_loc[pairs_from_correlation[i][1]][1]),2)+
				  pow( (ant_loc[pairs_from_correlation[i][0]][2] - ant_loc[pairs_from_correlation[i][1]][2]),2) );
	  delays_temp_norm_cross.push_back(BestTimes[i]/baseline_length*c*nsTos/n_ice * -1.);
	}

	if (pairs_from_correlation[i][0] == pairs_pol_faces_pair_ind[polarization][faceType][0][1] && pairs_from_correlation[i][1] == pairs_pol_faces_pair_ind[polarization][faceType][1][1]){
	  baseline_length = sqrt( pow( (ant_loc[pairs_from_correlation[i][0]][0] - ant_loc[pairs_from_correlation[i][1]][0]),2)+
				  pow( (ant_loc[pairs_from_correlation[i][0]][1] - ant_loc[pairs_from_correlation[i][1]][1]),2)+
				  pow( (ant_loc[pairs_from_correlation[i][0]][2] - ant_loc[pairs_from_correlation[i][1]][2]),2) );
	  delays_temp_norm_cross.push_back(BestTimes[i]/baseline_length*c*nsTos/n_ice);
	}
	
	if (pairs_from_correlation[i][1] == pairs_pol_faces_pair_ind[polarization][faceType][0][1] && pairs_from_correlation[i][0] == pairs_pol_faces_pair_ind[polarization][faceType][1][1]){
	  
	  baseline_length = sqrt( pow( (ant_loc[pairs_from_correlation[i][0]][0] - ant_loc[pairs_from_correlation[i][1]][0]),2)+
				  pow( (ant_loc[pairs_from_correlation[i][0]][1] - ant_loc[pairs_from_correlation[i][1]][1]),2)+
				  pow( (ant_loc[pairs_from_correlation[i][0]][2] - ant_loc[pairs_from_correlation[i][1]][2]),2) );
	  delays_temp_norm_cross.push_back(BestTimes[i]/baseline_length*c*nsTos/n_ice * -1.);
	}
	

      }
      //    DelaysForFaces.push_back(delays_temp);
      double error_standard = getRms(delays_temp_norm_standard);
      double error_cross = getRms(delays_temp_norm_cross);

      /*
      for (int i = 0; i < delays_temp_norm_standard.size(); i++){
	cout << delays_temp_norm_standard[i] << " : " << delays_temp_norm_cross[i] << endl;
      }
      */
      double error_combined = sqrt( pow(error_standard,2.) + pow(error_cross,2.));
      
      RmsValues.push_back(error_combined);
      //      cout << polarization << " : " << faceType << " : " << error_combined << " : " << error_standard << " : " << error_cross << endl;


    }
    //  }

  return RmsValues;
}

vector<double> getRms_Faces(
			    vector<double> BestTimes,
			    // vector<double> BestSignal,
			    int polarization,
			    vector<vector<vector<vector<int> > > > pairs_pol_faces_pair_ind, 
			    vector<vector<double> > ant_loc){

  vector<double> RmsValues;
  double dt1, dt2;
  double baseline_length1, baseline_length2;
  

  // double baseline_length = 0.;
  //  for (int polarization = 0; polarization < 2; polarization++){
    for (int faceType = 0; faceType < pairs_pol_faces_pair_ind[polarization].size(); faceType++){
      vector<double> delays_temp_norm_standard;
      vector<double> delays_temp_norm_cross;
    
      // standard pairs
      dt1 = BestTimes[pairs_pol_faces_pair_ind[polarization][faceType][0][0]] - BestTimes[pairs_pol_faces_pair_ind[polarization][faceType][0][1]];
      baseline_length1 = sqrt( 
			      pow( (ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][0][0]][0] - ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][0][1]][0]),2)+
			      pow( (ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][0][0]][1] - ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][0][1]][1]),2)+
			      pow( (ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][0][0]][2] - ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][0][1]][2]),2)
			       );
      delays_temp_norm_standard.push_back(dt1/baseline_length1*c*nsTos/n_ice);
      
      
      dt2 = BestTimes[pairs_pol_faces_pair_ind[polarization][faceType][1][0]] - BestTimes[pairs_pol_faces_pair_ind[polarization][faceType][1][1]];
      baseline_length2 = sqrt( 
			      pow( (ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][1][0]][0] - ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][1][1]][0]),2)+
			      pow( (ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][1][0]][1] - ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][1][1]][1]),2)+
			      pow( (ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][1][0]][2] - ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][1][1]][2]),2)
			       );
      delays_temp_norm_standard.push_back(dt2/baseline_length2*c*nsTos/n_ice);
      
      
      //cross pairs
      dt1 = BestTimes[pairs_pol_faces_pair_ind[polarization][faceType][0][0]] - BestTimes[pairs_pol_faces_pair_ind[polarization][faceType][1][0]];
      baseline_length1 = sqrt( 
			      pow( (ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][0][0]][0] - ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][1][0]][0]),2)+
			      pow( (ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][0][0]][1] - ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][1][0]][1]),2)+
			      pow( (ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][0][0]][2] - ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][1][0]][2]),2)
			       );
      delays_temp_norm_cross.push_back(dt1/baseline_length1*c*nsTos/n_ice);
      
      dt2 = BestTimes[pairs_pol_faces_pair_ind[polarization][faceType][0][1]] - BestTimes[pairs_pol_faces_pair_ind[polarization][faceType][1][1]];
      baseline_length1 = sqrt( 
			      pow( (ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][0][1]][0] - ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][1][1]][0]),2)+
			      pow( (ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][0][1]][1] - ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][1][1]][1]),2)+
			      pow( (ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][0][1]][2] - ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][1][1]][2]),2)
			       );
      delays_temp_norm_cross.push_back(dt2/baseline_length2*c*nsTos/n_ice);

    //    DelaysForFaces.push_back(delays_temp);
    double error_standard = getRms(delays_temp_norm_standard);
    double error_cross = getRms(delays_temp_norm_cross);
    
    /*
      for (int i = 0; i < delays_temp_norm_standard.size(); i++){
      cout << delays_temp_norm_standard[i] << " : " << delays_temp_norm_cross[i] << endl;
      }
    */
    double error_combined = sqrt( pow(error_standard,2.) + pow(error_cross,2.));
    
    RmsValues.push_back(error_combined);
    //      cout << polarization << " : " << faceType << " : " << error_combined << " : " << error_standard << " : " << error_cross << endl;

    }
  
    return RmsValues;
}




bool isCutOnFaceError(int nFaces, double* faceCuts, vector<double> faceErrors)
{
  bool isCut = true;
  for ( int i = 0; i < nFaces; i++){
    if (faceErrors[i] < faceCuts[i]){
      isCut = false;
    }
  }

  return isCut;
}


vector<double> getRms_Faces_Thresh(
				   vector<double> BestTimes, // vector of arrival times at channels
				   vector<double> BestSignals, // vector of signal strengths to compare with threshold found at arrival times
				   double threshold, // threshold of signal strengths needed to be considered of sufficient strength to be included in face
				   int polarization, // polarization of antennas to check
				   vector<vector<vector<vector<int> > > > pairs_pol_faces_pair_ind, // vector of antenna pairs for faces
				   vector<vector<double> > ant_loc // vector of antenna locations 
				   ){
  
  vector<double> RmsValues;
  double dt1, dt2;
  double baseline_length1, baseline_length2;
  

  // double baseline_length = 0.;
  //  for (int polarization = 0; polarization < 2; polarization++){
    for (int faceType = 0; faceType < pairs_pol_faces_pair_ind[polarization].size(); faceType++){
      double error_combined = 0.;

      if (BestSignals[pairs_pol_faces_pair_ind[polarization][faceType][0][0]] > threshold &&
	  BestSignals[pairs_pol_faces_pair_ind[polarization][faceType][0][1]] > threshold &&
	  BestSignals[pairs_pol_faces_pair_ind[polarization][faceType][1][0]] > threshold &&
	  BestSignals[pairs_pol_faces_pair_ind[polarization][faceType][1][1]] > threshold){
      
	vector<double> delays_temp_norm_standard;
	vector<double> delays_temp_norm_cross;
	
	// standard pairs
	dt1 = BestTimes[pairs_pol_faces_pair_ind[polarization][faceType][0][0]] - BestTimes[pairs_pol_faces_pair_ind[polarization][faceType][0][1]];
	baseline_length1 = sqrt( 
				pow( (ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][0][0]][0] - ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][0][1]][0]),2)+
				pow( (ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][0][0]][1] - ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][0][1]][1]),2)+
				pow( (ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][0][0]][2] - ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][0][1]][2]),2)
				 );
	delays_temp_norm_standard.push_back(dt1/baseline_length1*c*nsTos/n_ice);
	
	
	dt2 = BestTimes[pairs_pol_faces_pair_ind[polarization][faceType][1][0]] - BestTimes[pairs_pol_faces_pair_ind[polarization][faceType][1][1]];
	baseline_length2 = sqrt( 
				pow( (ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][1][0]][0] - ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][1][1]][0]),2)+
				pow( (ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][1][0]][1] - ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][1][1]][1]),2)+
				pow( (ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][1][0]][2] - ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][1][1]][2]),2)
				 );
	delays_temp_norm_standard.push_back(dt2/baseline_length2*c*nsTos/n_ice);
	
	
	//cross pairs
	dt1 = BestTimes[pairs_pol_faces_pair_ind[polarization][faceType][0][0]] - BestTimes[pairs_pol_faces_pair_ind[polarization][faceType][1][0]];
	baseline_length1 = sqrt( 
				pow( (ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][0][0]][0] - ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][1][0]][0]),2)+
				pow( (ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][0][0]][1] - ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][1][0]][1]),2)+
				pow( (ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][0][0]][2] - ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][1][0]][2]),2)
				 );
	delays_temp_norm_cross.push_back(dt1/baseline_length1*c*nsTos/n_ice);
	
	dt2 = BestTimes[pairs_pol_faces_pair_ind[polarization][faceType][0][1]] - BestTimes[pairs_pol_faces_pair_ind[polarization][faceType][1][1]];
	baseline_length1 = sqrt( 
				pow( (ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][0][1]][0] - ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][1][1]][0]),2)+
				pow( (ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][0][1]][1] - ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][1][1]][1]),2)+
				pow( (ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][0][1]][2] - ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][1][1]][2]),2)
				 );
	delays_temp_norm_cross.push_back(dt2/baseline_length2*c*nsTos/n_ice);
	
	//    DelaysForFaces.push_back(delays_temp);
	double error_standard = getRms(delays_temp_norm_standard);
	double error_cross = getRms(delays_temp_norm_cross);
	
	/*
	  for (int i = 0; i < delays_temp_norm_standard.size(); i++){
	  cout << delays_temp_norm_standard[i] << " : " << delays_temp_norm_cross[i] << endl;
	  }
	*/
	error_combined = sqrt( pow(error_standard,2.) + pow(error_cross,2.));
      }
      else {
	error_combined = 1000.;
      }
      RmsValues.push_back(error_combined);
    //      cout << polarization << " : " << faceType << " : " << error_combined << " : " << error_standard << " : " << error_cross << endl;

    }
  
    return RmsValues;
}

void getAbsMaximum_N(TGraph* plot, int nPeaks, double timeApart, vector<double> &xs, vector<double> &ys)
{
  xs.clear();
  ys.clear();

  int nPoints = plot->GetN();
  if (nPoints < nPeaks){
    cerr << "Number of points in waveform is fewer than the number of peaks requested. Decreasing number of peaks requested to match number points." << endl;
    nPeaks = nPoints;
  } 

  double x_temp, y_temp;
  double y_good, x_good;
  int test;
  double y_upper;

    for (int iPeak = 0; iPeak < nPeaks; iPeak++){
      y_good = -9.0E99;
      y_upper = 9.0E99;
      if (iPeak > 0){
	y_upper = ys[iPeak-1];
      }
      for (int i = 0; i < nPoints; i++){
	test = plot->GetPoint(i, x_temp, y_temp);
	if (iPeak == 0){
	  if (y_temp > y_good){
	    x_good = x_temp;
	    y_good = y_temp;
	  }
	} 
	else {
	  for (int iTimeTest = 0; iTimeTest < xs.size(); iTimeTest++){
	    if (y_temp > y_good && y_temp < y_upper && abs(x_temp - xs[iTimeTest]) > timeApart){
	      x_good = x_temp;
	      y_good = y_temp;
	    }
	  }
	}
      }
      xs.push_back(x_good);
      ys.push_back(y_good);
      //cout << iPeak << " : " << xs[iPeak] << " : " << ys[iPeak] << endl;
    }

  return;
}



void getAbsMaximum_N(vector<TGraph*> graphs, int nPeaks, double timeApart, vector<vector<double> > &xs, vector<vector<double> > &ys){

  //  vector<double> xs;
  // vector<double> ys;

  xs.clear();
  ys.clear();

  vector<double> xs_temp;
  vector<double> ys_temp;

  for (int i = 0; i < graphs.size(); i++){
    getAbsMaximum_N(graphs[i], nPeaks, timeApart, xs_temp, ys_temp);
    xs.push_back(xs_temp);
    ys.push_back(ys_temp);
  }
  
  //    cout << ys.size() << endl;
  //  cout << ys[0].size() << endl;


}

vector<double> getRms_Faces_Thresh_N(
				   vector<vector<double> > BestTimes, // vector of arrival times at channels
				   vector<vector<double> > BestSignals, // vector of signal strengths to compare with threshold found at arrival times
				   double threshold, // threshold of signal strengths needed to be considered of sufficient strength to be included in face
				   int polarization, // polarization of antennas to check
				   vector<vector<vector<vector<int> > > > pairs_pol_faces_pair_ind, // vector of antenna pairs for faces
				   vector<vector<double> > ant_loc // vector of antenna locations 
				   ){
  
  vector<double> RmsValues;
  double dt1, dt2;
  double baseline_length1, baseline_length2;
  

  // double baseline_length = 0.;
  //  for (int polarization = 0; polarization < 2; polarization++){
    for (int faceType = 0; faceType < pairs_pol_faces_pair_ind[polarization].size(); faceType++){
      double error_combined = 0.;
      
      vector<double> test_face;

      //      cout << BestSignals[0].size() << endl;

      for (int i1 = 0; i1< BestSignals[0].size(); i1++){
      for (int i2 = 0; i2< BestSignals[0].size(); i2++){
      for (int i3 = 0; i3< BestSignals[0].size(); i3++){
      for (int i4 = 0; i4< BestSignals[0].size(); i4++){
	

      if (BestSignals[pairs_pol_faces_pair_ind[polarization][faceType][0][0]][i1] > threshold &&
	  BestSignals[pairs_pol_faces_pair_ind[polarization][faceType][0][1]][i2] > threshold &&
	  BestSignals[pairs_pol_faces_pair_ind[polarization][faceType][1][0]][i3] > threshold &&
	  BestSignals[pairs_pol_faces_pair_ind[polarization][faceType][1][1]][i4] > threshold){
      
	vector<double> delays_temp_norm_standard;
	vector<double> delays_temp_norm_cross;
	
	// standard pairs
	dt1 = BestTimes[pairs_pol_faces_pair_ind[polarization][faceType][0][0]][i1] - BestTimes[pairs_pol_faces_pair_ind[polarization][faceType][0][1]][i2];
	baseline_length1 = sqrt( 
				pow( (ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][0][0]][0] - ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][0][1]][0]),2)+
				pow( (ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][0][0]][1] - ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][0][1]][1]),2)+
				pow( (ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][0][0]][2] - ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][0][1]][2]),2)
				 );
	delays_temp_norm_standard.push_back(dt1/baseline_length1*c*nsTos/n_ice);
	
	
	dt2 = BestTimes[pairs_pol_faces_pair_ind[polarization][faceType][1][0]][i3] - BestTimes[pairs_pol_faces_pair_ind[polarization][faceType][1][1]][i4];
	baseline_length2 = sqrt( 
				pow( (ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][1][0]][0] - ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][1][1]][0]),2)+
				pow( (ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][1][0]][1] - ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][1][1]][1]),2)+
				pow( (ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][1][0]][2] - ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][1][1]][2]),2)
				 );
	delays_temp_norm_standard.push_back(dt2/baseline_length2*c*nsTos/n_ice);
	
	
	//cross pairs
	dt1 = BestTimes[pairs_pol_faces_pair_ind[polarization][faceType][0][0]][i1] - BestTimes[pairs_pol_faces_pair_ind[polarization][faceType][1][0]][i3];
	baseline_length1 = sqrt( 
				pow( (ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][0][0]][0] - ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][1][0]][0]),2)+
				pow( (ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][0][0]][1] - ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][1][0]][1]),2)+
				pow( (ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][0][0]][2] - ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][1][0]][2]),2)
				 );
	delays_temp_norm_cross.push_back(dt1/baseline_length1*c*nsTos/n_ice);
	
	dt2 = BestTimes[pairs_pol_faces_pair_ind[polarization][faceType][0][1]][i2] - BestTimes[pairs_pol_faces_pair_ind[polarization][faceType][1][1]][i4];
	baseline_length1 = sqrt( 
				pow( (ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][0][1]][0] - ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][1][1]][0]),2)+
				pow( (ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][0][1]][1] - ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][1][1]][1]),2)+
				pow( (ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][0][1]][2] - ant_loc[pairs_pol_faces_pair_ind[polarization][faceType][1][1]][2]),2)
				 );
	delays_temp_norm_cross.push_back(dt2/baseline_length2*c*nsTos/n_ice);
	
	//    DelaysForFaces.push_back(delays_temp);
	double error_standard = getRms(delays_temp_norm_standard);
	double error_cross = getRms(delays_temp_norm_cross);
	
	/*
	  for (int i = 0; i < delays_temp_norm_standard.size(); i++){
	  cout << delays_temp_norm_standard[i] << " : " << delays_temp_norm_cross[i] << endl;
	  }
	*/
	error_combined = sqrt( pow(error_standard,2.) + pow(error_cross,2.));
      }
      else {
	error_combined = 1000.;
      }
      test_face.push_back(error_combined);
      }
      }
      }
      }

      /*
      cout << faceType << " : " << test_face.size() << " : ";
      for (int i = 0; i < test_face.size(); i++){
	cout << test_face[i] << " : ";
	
      }
      cout << endl;
      */

      sort(test_face.begin(), test_face.end());

      RmsValues.push_back(test_face[0]);
    //      cout << polarization << " : " << faceType << " : " << error_combined << " : " << error_standard << " : " << error_cross << endl;

    }
  
    return RmsValues;
}

double getPhase(FFTWComplex &theNum){
  return atan2(theNum.im, theNum.re);
}

TGraph *makeRawPhase(TGraph *grWave) {

  double *oldY = grWave->GetY();
  double *oldX = grWave->GetX();
  double deltaT=oldX[1]-oldX[0];
  int length=grWave->GetN();
  FFTWComplex *theFFT=FFTtools::doFFT(length,oldY);

  int newLength=(length/2)+1;
  double *newY = new double [newLength];
  double *newX = new double [newLength];

  double deltaF=1/(deltaT*length);
  //    double fMax = 1/(2*deltaT);  // In GHz                                                                                                                       

  double tempF=0;
  for(int i=0;i<newLength;i++) {
    float phase=getPhase(theFFT[i]);
    //    if(i>0 && i<newLength-1) power*=2; //account for symmetry                                                                                                        
    newX[i]=tempF;
    newY[i]=phase;
    tempF+=deltaF;
  }


  TGraph *grPhase = new TGraph(newLength,newX,newY);
  delete [] theFFT;
  delete [] newY;
  delete [] newX;
  return grPhase;

}


TGraph *getGraphDifference(TGraph *gr1, TGraph* gr2){

  TGraph* grOut = new TGraph();
  
  if (gr1->GetN() != gr2->GetN()){
    cerr << "Graphs not the same size!" << endl;
    return 0;
  } else {

    double *x1 = gr1->GetX();
    double *x2 = gr2->GetX();
    double *y1 = gr1->GetY();
    double *y2 = gr2->GetY();
    
    for (int i = 0; i< gr1->GetN(); i++){
      grOut->SetPoint(grOut->GetN(), x1[i], y1[i]-y2[i]);
    }
  }

  return grOut;

}

TGraph *getGraphSum(TGraph *gr1, TGraph* gr2){

  TGraph* grOut = new TGraph();
  
  if (gr1->GetN() != gr2->GetN()){
    cerr << "Graphs not the same size!" << endl;
    return 0;
  } else {

    double *x1 = gr1->GetX();
    double *x2 = gr2->GetX();
    double *y1 = gr1->GetY();
    double *y2 = gr2->GetY();
    
    for (int i = 0; i< gr1->GetN(); i++){
      grOut->SetPoint(grOut->GetN(), x1[i], y1[i]+y2[i]);
    }
  }

  return grOut;

}


TGraph *makeSpectrum_mVPerRootHz(TGraph *grWave){

  double *oldY = grWave->GetY();
  double *oldX = grWave->GetX();
  double deltaT=(oldX[1]-oldX[0])*1.e-9;
  int length=grWave->GetN();
  FFTWComplex *theFFT=FFTtools::doFFT(length,oldY);

  int newLength=(length/2)+1;

  double *newY = new double [newLength];
  double *newX = new double [newLength];

  //    double fMax = 1/(2*deltaT);  // In Hz                                                       
  double deltaF=1/(deltaT*length); //Hz                                                             
  deltaF*=1e-6; //MHz                                                                                
  //    std::cout << fMax << "\t" << deltaF << "\t" << deltaT << "\t" << length << std::endl;       

  double tempF=0;
  for(int i=0;i<newLength;i++) {
    newY[i] = FFTtools::getAbs(theFFT[i])*1.e-3;

			       /*
    float power=pow(FFTtools::getAbs(theFFT[i]),2);
    if(i>0 && i<newLength-1) power*=2; //account for symmetry                                       
    power*=deltaT/(length); //For time-integral squared amplitude                                   
    //    power*=(1e3*1e3)/1e9;                                                                     
    power/=deltaF;//Just to normalise bin-widths                                                    
    //Ends up the same as dt^2, need to integrate the power (multiply by df)                        
    //to get a meaningful number out.                                                               

    power=sqrt(power);

    //if (power>0 ) power=10*TMath::Log10(power);                                                   
    //else power=-1000; //no reason                                                                 

    */
    newX[i]=tempF;
			       //    newY[i]=power;
    tempF+=deltaF;
  }

  TGraph *grPower = new TGraph(newLength,newX,newY);
  delete [] theFFT;
  delete [] newY;
  delete [] newX;
  return grPower;

}

/*
void getWavefrontRMS(vector<TGraph*> graphs_raw, vector<double> waveformRMS, 
		       int stationID, vector< vector<double> > ant_loc,
		       double threshold,
		       double &wavefrontRMS_V, double &wavefrontRMS_H,
		       double interpolationTimeStep = 0.625,
		       double intergrationTime = 5.0,
		       int numSearchPeaks=2,
		       double peakSeparation=5.0
){
  //  double interpolationTimeStep = 0.625;
  //  int numBinsToIntegrate = (int)(5./interpolationTimeStep);
  int nGraphs = graphs_raw.size();
  if (graphs_raw.size() != waveformRMS.size()){ 
    cerr << "Graph vector and RMS vector not the same size!" << endl;
    return;}




  stringstream ss;
  string xLabel, yLabel;
  vector<string> titlesForGraphs;
  for (int i = 0; i < nGraphs; i++){
    ss.str("");
    ss << "Channel " << i;
    titlesForGraphs.push_back(ss.str());
  }

  xLabel = "Time (ns)"; yLabel = "Voltage (mV)";
  vector<TGraph*> grWaveformInt = makeInterpolatedGraphs(graphs_raw, interpolationTimeStep, yLabel, xLabel, titlesForGraphs);

  vector<double> vVPeak;
  vector<double> vVPeakTimes;
  
  getAbsMaximum(grWaveformsInt, vVPeakTimes, vVPeak);

  xLabel = "Time (ns)"; yLabel = "Integrated Power (arb units)";
  vector<TGraph*> grIntPower = makeIntegratedBinPowerGraphs(grWaveformsInt, numBinsToIntegrate, xLabel, yLabel, titlesForGraphs);


  vector<vector<double> > vvHitTimes;
  vector<vector<double> > vvPeakIntPowers;
  getAbsMaximum_N(grIntPower, numSearchPeaks, peakSeparation, vvHitTimes, vvPeakIntPowers);

  vector<vector<double> > vvRMS_10overRMS;
  vvRMS_10overRMS.resize(nGraphs);
  for (int i = 0; i < 16; i++){
    vvRMS_10overRMS[i].resize(vvPeakIntPowers[i].size());
    for (int j = 0; j < vvPeakIntPowers[i].size(); j++){
      vvRMS_10overRMS[i][j] = sqrt(vvPeakIntPowers[i][j]/numBinsToIntegrate)/waveformRMS[i];
    }
  }

  vector<vector<vector<vector<int> > > > faces = setupFaces(2);
  vector<double> faceRmsAllForReco_sorted;

  vector<double> rms_faces_V = getRms_Faces_Thresh_N(vvHitTimes, vvRMS_10overRMS, threshold, 0, faces, ant_loc);
  vector<double> rms_faces_H = getRms_Faces_Thresh_N(vvHitTimes, vvRMS_10overRMS, threshold, 1, faces, ant_loc);

  sort(rms_faces_V.begin(), rms_faces_V.end());
  sort(rms_faces_H.begin(), rms_faces_H.end());

  wavefrontRMS_V = log10(rms_faces_V[0]);
  wavefrontRMS_H = log10(rms_faces_H[0]);

  return;
}


*/
#endif
