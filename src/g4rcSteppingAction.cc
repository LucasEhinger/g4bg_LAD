#include "g4rcSteppingAction.hh"

#include <math.h>

#include "G4SystemOfUnits.hh"
#include "G4SteppingManager.hh"
#include "G4Track.hh"

g4rcSteppingAction::g4rcSteppingAction()
:drawFlag(false)
{

}

void g4rcSteppingAction::UserSteppingAction(const G4Step *aStep) {


	
	G4Track* aTrack = aStep->GetTrack();

	if(aTrack->GetTrackID() != 1 && aTrack->GetDefinition()->GetPDGEncoding() != 2212) {
		aTrack->SetTrackStatus(fStopAndKill);
	}
	if((aTrack->GetTrackID() != 1) 
	&& (aTrack->GetMomentumDirection().phi() < -20.*deg || aTrack->GetMomentumDirection().phi() > 20.*deg || 
		aTrack->GetMomentumDirection().theta() < 90.*deg || aTrack->GetMomentumDirection().theta() > 170.*deg)){
		aTrack->SetTrackStatus(fStopAndKill);
	}

/*
	G4Track* aTrack = aStep->GetTrack();

	if(aTrack->GetTrackID() != 1 || aTrack->GetTotalEnergy() < 0.1*GeV) {
		aTrack->SetTrackStatus(fStopAndKill);
	}
*/
}


