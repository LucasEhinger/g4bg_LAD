#include "g4rcIO.hh"

#include "G4SystemOfUnits.hh"
#include "G4AffineTransform.hh"

#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>

#include "G4ParticleDefinition.hh"

#include "g4rcDetectorHit.hh"
#include "g4rcEvent.hh"
#include "g4rcUniformScattering.hh"

#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <assert.h>


g4rcIO::g4rcIO(){
    fTree = NULL;
    InitializeTree();
    // Default filename
     	strcpy(fFilename, "g4rcout.root");
	fDetectedElectron = false;
	fHRSangle = 17.5*deg;
	fFile = NULL;
	fUS = NULL;
}

g4rcIO::~g4rcIO(){
    if( fTree ){ delete fTree; }
    fTree = NULL;
    if( fFile ){ delete fFile; }
    fFile = NULL;
}

void g4rcIO::SetFilename(G4String fn){
    G4cout << "Setting output file to " << fn << G4endl;
    strcpy(fFilename, fn.data());
}

void g4rcIO::InitializeTree(){
    if( fFile ){
	fFile->Close();
	delete fFile;
    }

    fFile = new TFile(fFilename, "RECREATE");

    if( fTree ){ delete fTree; }

    fTree = new TTree("T", "TRITIUM radiative effects simulation");

    // Event information
    fTree->Branch("ev.pid",   &fEvPart_PID, "ev.pid/I");
    fTree->Branch("ev.vx",    &fEvPart_X,   "ev.vx/D");
    fTree->Branch("ev.vy",    &fEvPart_Y,   "ev.vy/D");
    fTree->Branch("ev.p",     &fEvPart_P,   "ev.p/D");

	fTree->Branch("Ebeam",	&fEbeam,	"Ebeam/D");
	fTree->Branch("Epre",	&fEpre,		"Epre/D");
	fTree->Branch("E0",	&fE0,		"E0/D");
	fTree->Branch("Ep",	&fEp,		"Ep/D");
	fTree->Branch("Epost",	&fEpost,	"Epost/D");
	fTree->Branch("Edet",	&fEobs,		"Eded/D");		

	fTree->Branch("thHall.born",	&fTh0_HCS,	"thHall.born/D");	
	fTree->Branch("phHall.born",	&fPh0_HCS,	"phHall.born/D");	
	fTree->Branch("thTarg.born",	&fTh0_TCS,	"thHall.born/D");	
	fTree->Branch("phTarg.born",	&fPh0_TCS,	"phHall.born/D");	
	
	fTree->Branch("thHall.obs",	&fThObs_HCS,	"thHall.obs/D");	
	fTree->Branch("phHall.obs",	&fPhObs_HCS,	"phHall.obs/D");	
	fTree->Branch("thTarg.obs",	&fThObs_TCS,	"thHall.obs/D");	
	fTree->Branch("phTarg.obs",	&fPhObs_TCS,	"phHall.obs/D");	

	fTree->Branch("Q2.obs",		&fQ2obs,	"Q2.obs/D");
	fTree->Branch("xBj.obs",	&fxBobs,	"xBj.obs/D");
	fTree->Branch("Q2.born",	&fQ2born,	"Q2.born/D");
	fTree->Branch("xBj.born",	&fxBborn,	"xBj.born/D");

    // DetectorHit
    fTree->Branch("hit.n",    &fNDetHit,     "hit.n/I");
    fTree->Branch("hit.det",  &fDetHit_det,  "hit.det[hit.n]/I");

    fTree->Branch("hit.pid",  &fDetHit_pid,   "hit.pid[hit.n]/I");
    fTree->Branch("hit.trid", &fDetHit_trid,  "hit.trid[hit.n]/I");
    fTree->Branch("hit.mtrid",&fDetHit_mtrid, "hit.mtrid[hit.n]/I");

    fTree->Branch("hit.x",    &fDetHit_X,   "hit.x[hit.n]/D");
    fTree->Branch("hit.y",    &fDetHit_Y,   "hit.y[hit.n]/D");
    fTree->Branch("hit.z",    &fDetHit_Z,   "hit.z[hit.n]/D");

    fTree->Branch("hit.px",   &fDetHit_Px,   "hit.px[hit.n]/D");
    fTree->Branch("hit.py",   &fDetHit_Py,   "hit.py[hit.n]/D");
    fTree->Branch("hit.pz",   &fDetHit_Pz,   "hit.pz[hit.n]/D");

    fTree->Branch("hit.p",    &fDetHit_P,   "hit.p[hit.n]/D");
    fTree->Branch("hit.e",    &fDetHit_E,   "hit.e[hit.n]/D");
    fTree->Branch("hit.m",    &fDetHit_M,   "hit.m[hit.n]/D");

    return;
}

void g4rcIO::FillTree(){
    if( !fTree ){ 
	fprintf(stderr, "Error %s: %s line %d - Trying to fill non-existant tree\n", __PRETTY_FUNCTION__, __FILE__, __LINE__ );
	return; 
    }

    fTree->Fill();
}

void g4rcIO::Flush(){
    	//  Set arrays to 0
    	fNDetHit = 0;
	fDetectedElectron = false;
}

void g4rcIO::WriteTree(){
    assert( fFile );
    assert( fTree );

    if( !fFile->IsOpen() ){
	G4cerr << "ERROR: " << __FILE__ << " line " << __LINE__ << ": TFile not open" << G4endl;
	exit(1);
    }

    G4cout << "Writing output to " << fFile->GetName() << "... ";

    fFile->cd();

    fTree->Write("T", TObject::kOverwrite);

    fTree->ResetBranchAddresses();
    delete fTree;
    fTree = NULL;

    fFile->Close();

    delete fFile;
    fFile = NULL;

    G4cout << "written" << G4endl;

    return;
}

///////////////////////////////////////////////////////////////////////////////
// Interfaces to output section ///////////////////////////////////////////////

// Event Data

void g4rcIO::SetEventData(g4rcEvent *ev){
    fEvPart_PID = ev->fPartType[0]->GetPDGEncoding();

    fEvPart_X = ev->fPartPos[0].x()/__L_UNIT;
    fEvPart_Y = ev->fPartPos[0].y()/__L_UNIT;

    fEvPart_P = ev->fPartMom[0].mag()/__E_UNIT;

    return;
}


void g4rcIO::SetScatteringData() {


	G4double Mp = 0.938272;
	
	fEbeam = 10.589;

/*
 * The following variables (if they exist) are filled by AddDetectorHit:
 * 
 * 	fEobs
 * 	fThObs_HCS
 * 	fPhObs_TCS	
*/

	// Get energies (convert to GeV)
	fEpre = fUS->fEpre/GeV;
	fE0 = fUS->fE0/GeV;	
	fEp = fUS->fEp/GeV;
	fEpost = fUS->fEpost/GeV;

	G4double nuBorn = fE0 - fEp;

	// Get born event angles (hall coordinates)
	fTh0_HCS = fUS->fTheta;
	fPh0_HCS = fUS->fPhi;	

	// Get or calculate born event x and Q2
	fQ2born = fUS->fQ2born/(GeV*GeV);
	fxBborn = fQ2born/(2.*Mp*nuBorn);

	// Calculate born event angles in TRANSPORT coordinates
	G4ThreeVector p0_HCS = G4ThreeVector(sin(fTh0_HCS)*cos(fPh0_HCS), sin(fTh0_HCS)*sin(fPh0_HCS), cos(fPh0_HCS));

	G4RotationMatrix rotateTarg;
	rotateTarg.rotateY(-fHRSangle);
	rotateTarg.rotateZ(90.*deg);		
	G4AffineTransform target_transform;
	target_transform.SetNetRotation(rotateTarg);
	target_transform.Invert();

	G4ThreeVector p0_TCS = target_transform.TransformPoint(p0_HCS);

	fTh0_TCS = p0_TCS.x()/p0_TCS.z();
	fPh0_TCS = p0_TCS.y()/p0_TCS.z();

	// If the electron entered the spectrometer, calculate the observed quantities
	if(fDetectedElectron) {
		G4double nuObs = fEbeam - fEobs;
		fQ2obs = 2.*fEbeam*fEobs*(1. - cos(fThObs_HCS));
		fxBobs = fQ2obs/(2.*Mp*nuObs);
		G4ThreeVector pObs_HCS = G4ThreeVector(sin(fThObs_HCS)*cos(fPhObs_HCS), sin(fThObs_HCS)*sin(fPhObs_HCS), cos(fPhObs_HCS));
		G4ThreeVector pObs_TCS = target_transform.TransformPoint(pObs_HCS);
		fThObs_TCS = pObs_TCS.x()/pObs_TCS.z();
		fPhObs_TCS = pObs_TCS.y()/pObs_TCS.z();	
	} else {

		fQ2obs = fxBobs = fThObs_TCS = fPhObs_TCS = fThObs_HCS = fPhObs_HCS = fEobs = -333.;

	}

	
	
}

// DetectorHit

void g4rcIO::AddDetectorHit(g4rcDetectorHit *hit){
    int n = fNDetHit;
    //printf("%d hits in detector\n", fNDetHit );

    if( n >= __IO_MAXHIT ){
//	G4cerr << "WARNING: " << __PRETTY_FUNCTION__ << " line " << __LINE__ << ":  Buffer size exceeded!" << G4endl;
	return;
    }

    fDetHit_det[n]  = hit->fDetID;

    fDetHit_trid[n] = hit->fTrID;
    fDetHit_mtrid[n]= hit->fmTrID;
    fDetHit_pid[n]  = hit->fPID;

    fDetHit_X[n]  = hit->f3X.x()/__L_UNIT;
    fDetHit_Y[n]  = hit->f3X.y()/__L_UNIT;
    fDetHit_Z[n]  = hit->f3X.z()/__L_UNIT;

    fDetHit_Px[n]  = hit->f3P.x()/__E_UNIT;
    fDetHit_Py[n]  = hit->f3P.y()/__E_UNIT;
    fDetHit_Pz[n]  = hit->f3P.z()/__E_UNIT;

    fDetHit_P[n]  = hit->fP/__E_UNIT;
    fDetHit_E[n]  = hit->fE/__E_UNIT;
    fDetHit_M[n]  = hit->fM/__E_UNIT;

	if(hit->fTrID == 1) {
		fEobs = hit->fE/GeV;
		fThObs_HCS = hit->f3P.theta();
		fPhObs_HCS = hit->f3P.phi();
		fDetectedElectron = true;
	}

    fNDetHit++;

    return;
}


