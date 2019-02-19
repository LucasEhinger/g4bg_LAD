#ifndef g4rcIO_HH
#define g4rcIO_HH

#include "TROOT.h"
#include "TObject.h"
#include "G4Run.hh"
#include "g4rctypes.hh"

#include "G4String.hh"

#include "g4rcEvent.hh"

class TFile;
class TTree;

class g4rcDetectorHit;
class g4rcUniformScattering;

#define __IO_MAXHIT 10000
#define __FILENAMELEN 255

// Units for output
#define __E_UNIT GeV
#define __L_UNIT m
#define __T_UNIT ns
#define __ANG_UNIT rad
#define __ASYMM_SCALE 1e-9 // ppb

class g4rcIO {
    public:
	 g4rcIO();
	~g4rcIO();

	void SetFilename(G4String  fn);
	G4String GetFilename(){return fFilename;}

	void FillTree();
	void Flush();
	void WriteTree();

	void WriteRun();

	void InitializeTree();
	Double_t fHRSangle;	

    private:
	TFile *fFile;
	TTree *fTree;

	char fFilename[__FILENAMELEN];

	g4rcUniformScattering* fUS;

	//  Interfaces and buffers to the tree
	//  This is potentially going to get very long...

	// Event data
    public:
	void SetEventData(g4rcEvent *);
	void SetScatteringData();
	void SetUniformScatteringProcess(g4rcUniformScattering* us) { fUS = us; }
    private:
	Double_t fEvPart_X;
	Double_t fEvPart_Y;
	Double_t fEvPart_P;
	Int_t fEvPart_PID;

	Double_t fEbeam;
	Double_t fE0;
	Double_t fEp;
	Double_t fEpre;
	Double_t fEpost;
	Double_t fEobs;

	Double_t fTh0_HCS;
	Double_t fPh0_HCS;
	Double_t fTh0_TCS;
	Double_t fPh0_TCS;
	
	Double_t fThObs_HCS;
	Double_t fPhObs_HCS;
	Double_t fThObs_TCS;
	Double_t fPhObs_TCS;
	
	Double_t fQ2obs;
	Double_t fxBobs;
	Double_t fQ2born;
	Double_t fxBborn;

	bool fDetectedElectron;

	//  DetectorHit
    public:
	void AddDetectorHit(g4rcDetectorHit *);
    private:
	Int_t fNDetHit;
	Int_t fDetHit_det[__IO_MAXHIT];

	Int_t fDetHit_trid[__IO_MAXHIT];
	Int_t fDetHit_pid[__IO_MAXHIT];
	Int_t fDetHit_mtrid[__IO_MAXHIT];

	Double_t fDetHit_X[__IO_MAXHIT];
	Double_t fDetHit_Y[__IO_MAXHIT];
	Double_t fDetHit_Z[__IO_MAXHIT];

	Double_t fDetHit_Px[__IO_MAXHIT];
	Double_t fDetHit_Py[__IO_MAXHIT];
	Double_t fDetHit_Pz[__IO_MAXHIT];
	Double_t fDetHit_P[__IO_MAXHIT];
	Double_t fDetHit_E[__IO_MAXHIT];
	Double_t fDetHit_M[__IO_MAXHIT];


};

#endif//g4rcIO_HH
