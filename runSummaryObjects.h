#ifndef RUM_SUMMARY_OBJS
#define RUM_SUMMARY_OBJS

int runNumber;

int numEvents;
int numRFTriggers;
int numSoftTriggers;
int numCalpulsers;

double RMS_SoftTrigger[16];
double RMS_RFTrigger[16];
double RMS_Calpulser[16];
double RMS_All[16];

#endif
