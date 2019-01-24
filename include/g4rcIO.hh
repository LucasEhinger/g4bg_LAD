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
class g4rcScintDetectorHit;
class g4rcUniformScattering;
class g4rcCrossSection;

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

    private:
	TFile *fFile;
	TTree *fTree;

	char fFilename[__FILENAMELEN];

	g4rcUniformScattering* fUS;
	g4rcCrossSection* fXS;

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
	Double_t fEvPart_Z;
	Double_t fEvPart_P;
	Double_t fEvPart_Px;
	Double_t fEvPart_Py;
	Double_t fEvPart_Pz;
	Int_t fEvPart_PID;

	Double_t fEpre;
	Double_t fE0;
	Double_t fEp;
	Double_t fEpost;
	Double_t fTheta;
	Double_t fThTarg;
	Double_t fPhTarg;
	
	Double_t fQ2obs;
	Double_t fxBobs;
	Double_t fQ2true;
	Double_t fxBtrue;


	//  DetectorHit
    public:
	void AddDetectorHit(g4rcDetectorHit *);
    private:
	Int_t fNDetHit;
	Int_t fDetHit_det[__IO_MAXHIT];
	Int_t fDetHit_id[__IO_MAXHIT];

	Int_t fDetHit_trid[__IO_MAXHIT];
	Int_t fDetHit_pid[__IO_MAXHIT];
	Int_t fDetHit_gen[__IO_MAXHIT];
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

	Double_t fDetHit_Vx[__IO_MAXHIT];
	Double_t fDetHit_Vy[__IO_MAXHIT];
	Double_t fDetHit_Vz[__IO_MAXHIT];
	Double_t fDetHit_Vdx[__IO_MAXHIT];
	Double_t fDetHit_Vdy[__IO_MAXHIT];
	Double_t fDetHit_Vdz[__IO_MAXHIT];

	//  ScintDetectorHit
    public:
	void AddScintDetectorHit(g4rcScintDetectorHit *);
    private:
	Int_t fNScintDetHit;
	Int_t fScintDetHit_det[__IO_MAXHIT];
	Int_t fScintDetHit_id[__IO_MAXHIT];

	Double_t fScintDetHit_edep[__IO_MAXHIT];

};

#endif//g4rcIO_HH
