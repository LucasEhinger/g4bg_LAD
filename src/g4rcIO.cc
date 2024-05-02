#include "g4rcIO.hh"

#include "G4SystemOfUnits.hh"
#include "G4AffineTransform.hh"

#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>

#include "G4ParticleDefinition.hh"

#include "g4rcDetectorHit.hh"
#include "g4rcEvent.hh"

#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <assert.h>


g4rcIO::g4rcIO(){
    fTree = NULL;
    InitializeTree();
    // Default filename
     	strcpy(fFilename, "g4rcout.root");
	fHRSangle = 17.5*deg;
	fFile = NULL;
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
    fTree->Branch("ev.theta", &fEvPart_Theta, "ev.theta/D");
    fTree->Branch("ev.phi",   &fEvPart_Phi,   "ev.phi/D");

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

    fTree->Branch("hit.edep", &fDetHit_Edep, "hit.edep[hit.n]/D");
    fTree->Branch("hit.t",    &fDetHit_t,    "hit.t[hit.n]/D");

    return;
}

bool g4rcIO::HasProtonHit(){
    for( int i = 0; i < fNDetHit; i++ ){
    if( fDetHit_pid[i] == 2212 ){
        return true;
    }
    }

    return false;
}
void g4rcIO::FillTree(){
    if( !fTree ){ 
	fprintf(stderr, "Error %s: %s line %d - Trying to fill non-existant tree\n", __PRETTY_FUNCTION__, __FILE__, __LINE__ );
	return; 
    }

    if( g4rcIO::HasProtonHit() ){
        fTree->Fill();
    }
}

void g4rcIO::Flush(){
    	//  Set arrays to 0
    	fNDetHit = 0;
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
    fEvPart_Theta = ev->fPartMom[0].theta();
    fEvPart_Phi = ev->fPartMom[0].phi();

    fEvPart_P = ev->fPartMom[0].mag()/__E_UNIT;

    return;
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

    fDetHit_Edep[n] = hit->fEdep/__E_UNIT;
    fDetHit_t[n] = hit->fTime/__T_UNIT;

    fNDetHit++;

    return;
}


