#ifndef INPUT_PARAMS
#define INPUT_PARAMS

// Analysis Run Settings                                                                                                                        
int calpulserRunMode = 0;
int filterMode = 0;

double interpolationTimeStep = 0.6;
const int first50ns = (int)(50./interpolationTimeStep);
int numBinsToIntegrate = (int)(5./interpolationTimeStep);

double thresholdMin = 2.0;
double thresholdStep = 0.1;
const int thresholdSteps = 15;
int thresholdSteps_v = thresholdSteps;
const double thresholdEnd = thresholdMin + thresholdStep*(double)thresholdSteps;

int numSearchPeaks = 2;
const int numFaces = 12;
int numFaces_v = numFaces;

const int nGraphs = 16;

double recoFilterThreshold[] = {2.5, 2.5};
int recoFilterThresholdBin[] = {5, 5};
double recoFilterWavefrontRMSCut[] = {-1.0, -1.0};

const int numRadii = 35;
int numRadii_v = numRadii;
const int numPols = 2;
int numPols_v = numPols;
const int numPolMaps = 2;
const int sol_num = 0;

double radii[] = {16, 19, 22, 26,
			30, 35, 41, 47, 55, 65, 75, 88, 102, 119, 139, 162, 189, 221, 257,
			300, 350, 408, 475, 554, 646, 754, 879, 1024, 1194, 1392, 1623, 1892, 2207, 2573,
			3000};

//double radii[] = {30, 35, 41, 47, 55, 
//			65, 75, 88, 102, 119, 
//			139, 162, 189, 221, 257, 
//			300};

double radNum[] = {0, 1, 2, 3,
			  4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18,
			  19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33,
			  34};
//  double radNum[] = {0, 1, 2, 3, 4,
//			   5, 6, 7, 8, 9,
//			   10, 11, 12, 13, 14,
//			   15};

const double angularBinSize = 1.;
const int RTTestMode = 4;



#endif

