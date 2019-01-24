#include "g4rcMessenger.hh"

#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADouble.hh"

#include "g4rcDetectorConstruction.hh"
#include "g4rcIO.hh"
#include "g4rcEventAction.hh"
#include "g4rcPrimaryGeneratorAction.hh"
#include "g4rcSteppingAction.hh"
#include "g4rcCrossSection.hh"

#include "G4UImanager.hh"
#include "G4UIdirectory.hh"
#include "G4RunManager.hh"

#include "G4VPhysicalVolume.hh"

#include <iostream>

g4rcMessenger::g4rcMessenger(){
    /*  Initialize all the things it talks to to NULL */

	fIO           = NULL;
	fdetcon       = NULL;
	fevact        = NULL;
	fprigen       = NULL;
	fStepAct      = NULL;
	fXS	      = NULL;

	fRemollDir = new G4UIdirectory("/g4rc/");
	fRemollDir->SetGuidance("UI commands of this code");

	fileCmd = new G4UIcmdWithAString("/g4rc/filename",this);
	fileCmd->SetGuidance("Output filename");
	fileCmd->SetParameterName("filename", false);

	seedCmd = new G4UIcmdWithAnInteger("/g4rc/seed",this);
	seedCmd->SetGuidance("Set random engine seed");
	seedCmd->SetParameterName("seed", false);

	angleCmd = new G4UIcmdWithADoubleAndUnit("/g4rc/angle",this);
	angleCmd->SetGuidance("Set HRS angle");
	angleCmd->SetParameterName("angle",false);	

	targCmd = new G4UIcmdWithAString("/g4rc/target",this);
	targCmd->SetGuidance("Set target");
	targCmd->SetParameterName("targ",false);	

	hrsCmd = new G4UIcmdWithAString("/g4rc/hrs",this);
	hrsCmd->SetGuidance("Set HRS arm");
	hrsCmd->SetParameterName("hrs",false);	


}

		
g4rcMessenger::~g4rcMessenger(){
}


void g4rcMessenger::SetNewValue(G4UIcommand* cmd, G4String newValue){
	if( cmd == fileCmd ){
		fIO->SetFilename(newValue);
	}

	if( cmd == seedCmd ){
		G4int seed = seedCmd->GetNewIntValue(newValue);
		CLHEP::HepRandom::setTheSeed(seed);
	}

	if (cmd == angleCmd ) {
		G4double angle = angleCmd->GetNewDoubleValue(newValue);
		fdetcon->fHRSAngle = angle;
		fevact->fHRSAngle = angle;
		fIO->fHRSangle = angle;	
	}

	if (cmd == targCmd) {
		fdetcon->SetTarget(newValue);
		fXS->SetTarget(newValue);
	}

	if (cmd == hrsCmd) {
		fdetcon->fHRS = newValue;
	}
		
}












