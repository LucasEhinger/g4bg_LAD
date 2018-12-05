#include "g4rcSteppingAction.hh"

#include <math.h>

#include "G4SystemOfUnits.hh"
#include "G4SteppingManager.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "G4UserLimits.hh"

g4rcSteppingAction::g4rcSteppingAction()
:drawFlag(false)
{

}

void g4rcSteppingAction::UserSteppingAction(const G4Step *aStep) {

/*
	double vertex_z = 0.*cm;
	G4Track* aTrack = aStep->GetTrack();
	int stepN = aTrack->GetCurrentStepNumber();
	double stepL = aStep->GetStepLength();
	
	if(aTrack->GetParentID() == 0) {

		double z = aTrack->GetPosition().z();
		G4LogicalVolume* vol_log = aTrack->GetVolume()->GetLogicalVolume();
		G4String vol_name = vol_log->GetName();
		G4UserLimits* fLimits = vol_log->GetUserLimits();

		if (vol_name == "gas_logical") {
			Scatter(aTrack);
		}

		if(vol_name == "gas_logical" && !has_scattered) {

			if (fabs(z - vertex_z) < 1.e-5) {
				Scatter(aTrack);
				fLimits->SetMaxAllowedStep(5.*cm);
			} else {
				double delta_z = vertex_z - z;
				fLimits->SetMaxAllowedStep(delta_z);
			}

		}	
	}
*/

}

void g4rcSteppingAction::Scatter(G4Track* fTrack) {

	double theta = 45.*deg;
	fTrack->SetMomentumDirection(G4ThreeVector(sin(theta),0.,cos(theta)));
	has_scattered = true;			

}

