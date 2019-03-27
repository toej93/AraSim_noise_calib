#include "inputParameters.h"

#ifndef OUTPUT_OBJS
#define OUTPUT_OBJS


// Event summary parameters
double TSQualParam;
bool isCalpulser;
bool isSoftTrigger;
bool isSimulation = false;
double thirdVPeakOverRMS[3];

// Waveform - channel summary parameters
double VPeak[16];
double waveformRMS[16];
double waveformRMS_50ns[16];
double VPeakOverRMS[16];
double VPeakTimes[16];
double avgPeakPower_5ns[16];
double peakPowerTimes[16];
int waveformLength[16];

// Polarization parameters
double polarizationRatio;

// Face RMS
//double rms_faces[2][12];
double rms_pol_thresh_face[2][thresholdSteps][numFaces];
//double thirdLowestFaceRms[3];

// Geometry parameters
double detectorCenter[3];

// Simulation parameters
int flavor;
int nu_nubar;
double energy;
double posnu[3];
double weight_out;
double viewAngle[16][2];
double viewAngleAvg[2];

//Reconstruction parameters
bool runReconstruction;
bool runCorrelation;

double minChi2Vals[numPols][numRadii];
double minChi2Phi[numPols][numRadii];
double minChi2Theta[numPols][numRadii];

double peakCorr[2][numRadii];
int peakTheta[2][numRadii];
int peakPhi[2][numRadii];
double minCorr[2][numRadii];
double meanCorr[2][numRadii];
double rmsCorr[2][numRadii];
double peakSigma[2][numRadii];

double peakCorr_single[2];
int peakTheta_single[2];
int peakPhi_single[2];
double minCorr_single[2];
double meanCorr_single[2];
double rmsCorr_single[2];
double peakSigma_single[2];

double peakCorr_single_select[2];
int peakTheta_single_select[2];
int peakPhi_single_select[2];
double minCorr_single_select[2];
double meanCorr_single_select[2];
double rmsCorr_single_select[2];
double peakSigma_single_select[2];


// Abandoned terms

//double peakCorrelationTimes[2][28];
//double peakCorrelationValues[2][28];

//double RMS_threshold_relativeRMS_weightedbyrelativepower_runAvg[2][100];
//double RMS_threshold_relativeRMS_weightedbyrelativepowersquared_runAvg[2][100];

//double CosThetaRMS_Corr_Threshold[2][30][9];
//double CosThetaRMS_Corr_Threshold_Remove1[2][30][9];
//int removedPairs[2][30][9];

#endif


//double RMS_relativePower[2];
//double RMS_relativePower2[2];

//double RMS_relativePower_50[2];
//double RMS_relativePower2_50[2];

//double RMS_relativePower_100[2];
//double RMS_relativePower2_100[2];

//double RMS_relativePower_200[2];
//double RMS_relativePower2_200[2];

//double RMS_threshold_sigma[2][100];
//double RMS_threshold_sigma_volts[2][100];

//double RMS_relativeRMS_100[2];
//double RMS_threshold_relativeRMS[2][100];

//double RMS_power[2];
//double RMS_powersquared[2];
//double RMS_fractional[2];

//double RMS_threshold_relativeRMS_weightedbypower[2][100];
//double RMS_threshold_relativeRMS_weightedbypowersquared[2][100];
//double RMS_threshold_relativeRMS_weightedbyrelativepower_100[2][100];
//double RMS_threshold_relativeRMS_weightedbyrelativepowersquared_100[2][100];
