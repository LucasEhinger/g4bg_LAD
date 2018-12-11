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

	if(aTrack->GetTrackID() != 1 || aTrack->GetTotalEnergy() < 1.*GeV) {
		aTrack->SetTrackStatus(fStopAndKill);
	}

}


