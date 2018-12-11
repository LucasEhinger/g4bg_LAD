#include "g4rcEventAction.hh"
#include "g4rcDetectorHit.hh"
#include "g4rcScintDetectorHit.hh"
#include "g4rcSteppingAction.hh"
#include "g4rcUniformScatteringConstructor.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
//#include "Randomize.hh"

#include "CLHEP/Random/RandFlat.h"


#include "g4rcIO.hh"


g4rcEventAction::g4rcEventAction() {
	fUS = NULL;
}

g4rcEventAction::~g4rcEventAction(){
}


void g4rcEventAction::BeginOfEventAction(const G4Event*ev) {
    // Pretty ongoing status with flush
    if( (ev->GetEventID() % 1000) == 0 ){
	printf("Event %8d\r", ev->GetEventID() );
	fflush(stdout);
    }

	if(!fUS) {	
		fUS = fUSC->GetUniformScatteringProcess();
		fIO->SetUniformScatteringProcess(fUS);
		fUS->SetCentralScatteringAngle(fHRSAngle);
	}


	G4double vz = CLHEP::RandFlat::shoot(-10.*cm, 10.*cm);
	fUS->SetVertexZ(vz);
	fUS->fHasScattered = false;

    return;
}

void g4rcEventAction::EndOfEventAction(const G4Event* evt ) {
  //G4SDManager   *SDman = G4SDManager::GetSDMpointer();
  G4HCofThisEvent *HCE = evt->GetHCofThisEvent();

  G4VHitsCollection *thiscol;

	fIO->SetScatteringData();

if(HCE) {
  // Traverse all hit collections, sort by output type
  for( int hcidx = 0; hcidx < HCE->GetCapacity(); hcidx++ ){
      thiscol = HCE->GetHC(hcidx);
      if(thiscol){ // This is NULL if nothing is stored
	  // Dyanmic cast to test types, process however see fit and feed to IO
	  
	  ////  Detector Hits ///////////////////////////////////
	  if( g4rcDetectorHitsCollection *thiscast = 
		  dynamic_cast<g4rcDetectorHitsCollection *>(thiscol)){
	      for( unsigned int hidx = 0; hidx < thiscast->GetSize(); hidx++ ){
		  fIO->AddDetectorHit(
			  (g4rcDetectorHit *) thiscast->GetHit(hidx) );
	      }
	  }
	  
	  ////  Scint Generic Detector Hits ///////////////////////////////////
	  if( g4rcScintDetectorHitsCollection *thiscast =
		  dynamic_cast<g4rcScintDetectorHitsCollection *>(thiscol)){
	      for( unsigned int hidx = 0; hidx < thiscast->GetSize(); hidx++ ){
		  fIO->AddScintDetectorHit(
			  (g4rcScintDetectorHit *) thiscast->GetHit(hidx) );
	      }
	  }


      }
  }
}


  // Fill tree and reset buffers
  fIO->FillTree();
  fIO->Flush();

  return;
}



