/********************************************************************
* eventSimDict.h
* CAUTION: DON'T CHANGE THIS FILE. THIS FILE IS AUTOMATICALLY GENERATED
*          FROM HEADER FILES LISTED IN G__setup_cpp_environmentXXX().
*          CHANGE THOSE HEADER FILES AND REGENERATE THIS FILE.
********************************************************************/
#ifdef __CINT__
#error eventSimDict.h/C is only for compilation. Abort cint.
#endif
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define G__ANSIHEADER
#define G__DICTIONARY
#define G__PRIVATE_GVALUE
#include "G__ci.h"
#include "FastAllocString.h"
extern "C" {
extern void G__cpp_setup_tagtableeventSimDict();
extern void G__cpp_setup_inheritanceeventSimDict();
extern void G__cpp_setup_typetableeventSimDict();
extern void G__cpp_setup_memvareventSimDict();
extern void G__cpp_setup_globaleventSimDict();
extern void G__cpp_setup_memfunceventSimDict();
extern void G__cpp_setup_funceventSimDict();
extern void G__set_cpp_environmenteventSimDict();
}


#include "TObject.h"
#include "TMemberInspector.h"
#include "Trigger.h"
#include "Detector.h"
#include "Settings.h"
#include "Spectra.h"
#include "IceModel.h"
#include "Primaries.h"
#include "Report.h"
#include "Event.h"
#include "RawIcrrStationHeader.h"
#include "RawIcrrStationEvent.h"
#include "RawAraStationEvent.h"
#include "FullIcrrHkEvent.h"
#include "AraEventCalibrator.h"
#include "IcrrTriggerMonitor.h"
#include "IcrrHkData.h"
#include "AraRawIcrrRFChannel.h"
#include "AraAntennaInfo.h"
#include "AraGeomTool.h"
#include "UsefulAraStationEvent.h"
#include "UsefulIcrrStationEvent.h"
#include "araIcrrStructures.h"
#include "AtriEventHkData.h"
#include "AtriSensorHkData.h"
#include "araAtriStructures.h"
#include "araSoft.h"
#include "araIcrrDefines.h"
#include "RawAtriSimpleStationEvent.h"
#include "RawAtriStationBlock.h"
#include "RawAraGenericHeader.h"
#include "RawAtriStationEvent.h"
#include "UsefulAtriStationEvent.h"
#include "AraStationInfo.h"
#include <algorithm>
namespace std { }
using namespace std;

#ifndef G__MEMFUNCBODY
#endif

extern G__linked_taginfo G__eventSimDictLN_TClass;
extern G__linked_taginfo G__eventSimDictLN_TBuffer;
extern G__linked_taginfo G__eventSimDictLN_TMemberInspector;
extern G__linked_taginfo G__eventSimDictLN_string;
extern G__linked_taginfo G__eventSimDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR;
extern G__linked_taginfo G__eventSimDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__eventSimDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR;
extern G__linked_taginfo G__eventSimDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__eventSimDictLN_Detector;
extern G__linked_taginfo G__eventSimDictLN_Settings;
extern G__linked_taginfo G__eventSimDictLN_Report;
extern G__linked_taginfo G__eventSimDictLN_Trigger;
extern G__linked_taginfo G__eventSimDictLN_vectorlEdoublecOallocatorlEdoublegRsPgR;
extern G__linked_taginfo G__eventSimDictLN_reverse_iteratorlEvectorlEdoublecOallocatorlEdoublegRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__eventSimDictLN_vectorlEvectorlEdoublecOallocatorlEdoublegRsPgRcOallocatorlEvectorlEdoublecOallocatorlEdoublegRsPgRsPgRsPgR;
extern G__linked_taginfo G__eventSimDictLN_reverse_iteratorlEvectorlEvectorlEdoublecOallocatorlEdoublegRsPgRcOallocatorlEvectorlEdoublecOallocatorlEdoublegRsPgRsPgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__eventSimDictLN_vectorlEvectorlEvectorlEdoublecOallocatorlEdoublegRsPgRcOallocatorlEvectorlEdoublecOallocatorlEdoublegRsPgRsPgRsPgRcOallocatorlEvectorlEvectorlEdoublecOallocatorlEdoublegRsPgRcOallocatorlEvectorlEdoublecOallocatorlEdoublegRsPgRsPgRsPgRsPgRsPgR;
extern G__linked_taginfo G__eventSimDictLN_reverse_iteratorlEvectorlEvectorlEvectorlEdoublecOallocatorlEdoublegRsPgRcOallocatorlEvectorlEdoublecOallocatorlEdoublegRsPgRsPgRsPgRcOallocatorlEvectorlEvectorlEdoublecOallocatorlEdoublegRsPgRcOallocatorlEvectorlEdoublecOallocatorlEdoublegRsPgRsPgRsPgRsPgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__eventSimDictLN_Vector;
extern G__linked_taginfo G__eventSimDictLN_Position;
extern G__linked_taginfo G__eventSimDictLN_IceModel;
extern G__linked_taginfo G__eventSimDictLN_TF1;
extern G__linked_taginfo G__eventSimDictLN_Parameters;
extern G__linked_taginfo G__eventSimDictLN_Surface_antenna;
extern G__linked_taginfo G__eventSimDictLN_Antenna;
extern G__linked_taginfo G__eventSimDictLN_Antenna_string;
extern G__linked_taginfo G__eventSimDictLN_vectorlEAntennacOallocatorlEAntennagRsPgR;
extern G__linked_taginfo G__eventSimDictLN_reverse_iteratorlEvectorlEAntennacOallocatorlEAntennagRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__eventSimDictLN_ARA_station;
extern G__linked_taginfo G__eventSimDictLN_vectorlEAntenna_stringcOallocatorlEAntenna_stringgRsPgR;
extern G__linked_taginfo G__eventSimDictLN_reverse_iteratorlEvectorlEAntenna_stringcOallocatorlEAntenna_stringgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__eventSimDictLN_vectorlESurface_antennacOallocatorlESurface_antennagRsPgR;
extern G__linked_taginfo G__eventSimDictLN_reverse_iteratorlEvectorlESurface_antennacOallocatorlESurface_antennagRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__eventSimDictLN_InstalledStation;
extern G__linked_taginfo G__eventSimDictLN_vectorlEintcOallocatorlEintgRsPgR;
extern G__linked_taginfo G__eventSimDictLN_reverse_iteratorlEvectorlEintcOallocatorlEintgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__eventSimDictLN_vectorlEvectorlEintcOallocatorlEintgRsPgRcOallocatorlEvectorlEintcOallocatorlEintgRsPgRsPgRsPgR;
extern G__linked_taginfo G__eventSimDictLN_reverse_iteratorlEvectorlEvectorlEintcOallocatorlEintgRsPgRcOallocatorlEvectorlEintcOallocatorlEintgRsPgRsPgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__eventSimDictLN_IdealStation;
extern G__linked_taginfo G__eventSimDictLN_vectorlEARA_stationcOallocatorlEARA_stationgRsPgR;
extern G__linked_taginfo G__eventSimDictLN_reverse_iteratorlEvectorlEARA_stationcOallocatorlEARA_stationgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__eventSimDictLN_vectorlEInstalledStationcOallocatorlEInstalledStationgRsPgR;
extern G__linked_taginfo G__eventSimDictLN_reverse_iteratorlEvectorlEInstalledStationcOallocatorlEInstalledStationgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__eventSimDictLN_vectorlEIdealStationcOallocatorlEIdealStationgRsPgR;
extern G__linked_taginfo G__eventSimDictLN_reverse_iteratorlEvectorlEIdealStationcOallocatorlEIdealStationgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__eventSimDictLN_TVectorTlEfloatgR;
extern G__linked_taginfo G__eventSimDictLN_TVectorTlEdoublegR;
extern G__linked_taginfo G__eventSimDictLN_TGraph;
extern G__linked_taginfo G__eventSimDictLN_TRandom3;
extern G__linked_taginfo G__eventSimDictLN_Spectra;
extern G__linked_taginfo G__eventSimDictLN_Primaries;
extern G__linked_taginfo G__eventSimDictLN_Secondaries;
extern G__linked_taginfo G__eventSimDictLN_EarthModel;
extern G__linked_taginfo G__eventSimDictLN_Interaction;
extern G__linked_taginfo G__eventSimDictLN_TMatrixTBaselEfloatgR;
extern G__linked_taginfo G__eventSimDictLN_TMatrixTBaselEdoublegR;
extern G__linked_taginfo G__eventSimDictLN_TH2D;
extern G__linked_taginfo G__eventSimDictLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR;
extern G__linked_taginfo G__eventSimDictLN_TF2;
extern G__linked_taginfo G__eventSimDictLN_pairlEunsignedsPintcOintgR;
extern G__linked_taginfo G__eventSimDictLN_vectorlEpairlEunsignedsPintcOintgRcOallocatorlEpairlEunsignedsPintcOintgRsPgRsPgR;
extern G__linked_taginfo G__eventSimDictLN_reverse_iteratorlEvectorlEpairlEunsignedsPintcOintgRcOallocatorlEpairlEunsignedsPintcOintgRsPgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__eventSimDictLN_TCanvas;
extern G__linked_taginfo G__eventSimDictLN_TF3;
extern G__linked_taginfo G__eventSimDictLN_Signal;
extern G__linked_taginfo G__eventSimDictLN_RaySolver;
extern G__linked_taginfo G__eventSimDictLN_Y;
extern G__linked_taginfo G__eventSimDictLN_UsefulIcrrStationEvent;
extern G__linked_taginfo G__eventSimDictLN_UsefulAtriStationEvent;
extern G__linked_taginfo G__eventSimDictLN_Event;
extern G__linked_taginfo G__eventSimDictLN_Surface_antenna_r;
extern G__linked_taginfo G__eventSimDictLN_Antenna_r;
extern G__linked_taginfo G__eventSimDictLN_vectorlEPositioncOallocatorlEPositiongRsPgR;
extern G__linked_taginfo G__eventSimDictLN_reverse_iteratorlEvectorlEPositioncOallocatorlEPositiongRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__eventSimDictLN_String_r;
extern G__linked_taginfo G__eventSimDictLN_vectorlEAntenna_rcOallocatorlEAntenna_rgRsPgR;
extern G__linked_taginfo G__eventSimDictLN_reverse_iteratorlEvectorlEAntenna_rcOallocatorlEAntenna_rgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__eventSimDictLN_Station_r;
extern G__linked_taginfo G__eventSimDictLN_vectorlEString_rcOallocatorlEString_rgRsPgR;
extern G__linked_taginfo G__eventSimDictLN_reverse_iteratorlEvectorlEString_rcOallocatorlEString_rgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__eventSimDictLN_vectorlESurface_antenna_rcOallocatorlESurface_antenna_rgRsPgR;
extern G__linked_taginfo G__eventSimDictLN_reverse_iteratorlEvectorlESurface_antenna_rcOallocatorlESurface_antenna_rgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__eventSimDictLN_vectorlETGraphmUcOallocatorlETGraphmUgRsPgR;
extern G__linked_taginfo G__eventSimDictLN_reverse_iteratorlEvectorlETGraphmUcOallocatorlETGraphmUgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__eventSimDictLN_vectorlEStation_rcOallocatorlEStation_rgRsPgR;
extern G__linked_taginfo G__eventSimDictLN_reverse_iteratorlEvectorlEStation_rcOallocatorlEStation_rgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__eventSimDictLN_vectorlEInteractioncOallocatorlEInteractiongRsPgR;
extern G__linked_taginfo G__eventSimDictLN_reverse_iteratorlEvectorlEInteractioncOallocatorlEInteractiongRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__eventSimDictLN_AraCalType;
extern G__linked_taginfo G__eventSimDictLN_AraCalTypecLcLEAraCalType;
extern G__linked_taginfo G__eventSimDictLN_AraAntType;
extern G__linked_taginfo G__eventSimDictLN_AraAntTypecLcLEAraAntType;
extern G__linked_taginfo G__eventSimDictLN_AraAntPol;
extern G__linked_taginfo G__eventSimDictLN_AraAntPolcLcLEAraAntPol;
extern G__linked_taginfo G__eventSimDictLN_AraDaqChanType;
extern G__linked_taginfo G__eventSimDictLN_AraDaqChanTypecLcLEAraDaqChanType;
extern G__linked_taginfo G__eventSimDictLN_AraLabChip;
extern G__linked_taginfo G__eventSimDictLN_AraLabChipcLcLEAraLabChip;
extern G__linked_taginfo G__eventSimDictLN_AraAntDir;
extern G__linked_taginfo G__eventSimDictLN_AraAntDircLcLEAraAntDir;
extern G__linked_taginfo G__eventSimDictLN_TElementActionTlEfloatgR;
extern G__linked_taginfo G__eventSimDictLN_TElementPosActionTlEfloatgR;
extern G__linked_taginfo G__eventSimDictLN_TMatrixTlEfloatgR;
extern G__linked_taginfo G__eventSimDictLN_TMatrixTRow_constlEfloatgR;
extern G__linked_taginfo G__eventSimDictLN_TMatrixTRowlEfloatgR;
extern G__linked_taginfo G__eventSimDictLN_TMatrixTDiag_constlEfloatgR;
extern G__linked_taginfo G__eventSimDictLN_TMatrixTColumn_constlEfloatgR;
extern G__linked_taginfo G__eventSimDictLN_TMatrixTFlat_constlEfloatgR;
extern G__linked_taginfo G__eventSimDictLN_TMatrixTSub_constlEfloatgR;
extern G__linked_taginfo G__eventSimDictLN_TMatrixTSparseRow_constlEfloatgR;
extern G__linked_taginfo G__eventSimDictLN_TMatrixTSparseDiag_constlEfloatgR;
extern G__linked_taginfo G__eventSimDictLN_TMatrixTColumnlEfloatgR;
extern G__linked_taginfo G__eventSimDictLN_TMatrixTDiaglEfloatgR;
extern G__linked_taginfo G__eventSimDictLN_TMatrixTFlatlEfloatgR;
extern G__linked_taginfo G__eventSimDictLN_TMatrixTSublEfloatgR;
extern G__linked_taginfo G__eventSimDictLN_TMatrixTSparseRowlEfloatgR;
extern G__linked_taginfo G__eventSimDictLN_TMatrixTSparseDiaglEfloatgR;
extern G__linked_taginfo G__eventSimDictLN_vectorlEAraAntennaInfocOallocatorlEAraAntennaInfogRsPgR;
extern G__linked_taginfo G__eventSimDictLN_reverse_iteratorlEvectorlEAraAntennaInfocOallocatorlEAraAntennaInfogRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__eventSimDictLN_vectorlEAraCalAntennaInfocOallocatorlEAraCalAntennaInfogRsPgR;
extern G__linked_taginfo G__eventSimDictLN_reverse_iteratorlEvectorlEAraCalAntennaInfocOallocatorlEAraCalAntennaInfogRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__eventSimDictLN_vectorlEunsignedsPshortcOallocatorlEunsignedsPshortgRsPgR;
extern G__linked_taginfo G__eventSimDictLN_reverse_iteratorlEvectorlEunsignedsPshortcOallocatorlEunsignedsPshortgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__eventSimDictLN_vectorlEvectorlEunsignedsPshortcOallocatorlEunsignedsPshortgRsPgRcOallocatorlEvectorlEunsignedsPshortcOallocatorlEunsignedsPshortgRsPgRsPgRsPgR;
extern G__linked_taginfo G__eventSimDictLN_reverse_iteratorlEvectorlEvectorlEunsignedsPshortcOallocatorlEunsignedsPshortgRsPgRcOallocatorlEvectorlEunsignedsPshortcOallocatorlEunsignedsPshortgRsPgRsPgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__eventSimDictLN_vectorlERawAtriStationBlockcOallocatorlERawAtriStationBlockgRsPgR;
extern G__linked_taginfo G__eventSimDictLN_reverse_iteratorlEvectorlERawAtriStationBlockcOallocatorlERawAtriStationBlockgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__eventSimDictLN_maplEintcOvectorlEdoublecOallocatorlEdoublegRsPgRcOlesslEintgRcOallocatorlEpairlEconstsPintcOvectorlEdoublecOallocatorlEdoublegRsPgRsPgRsPgRsPgR;

/* STUB derived class for protected member access */
