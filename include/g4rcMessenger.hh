#ifndef g4rcMessenger_HH
#define g4rcMessenger_HH

#include "globals.hh"
#include "g4rctypes.hh"
#include "G4UImessenger.hh"
#include "G4UIcommand.hh"
#include "G4VModularPhysicsList.hh"

/*!
 *   Global messenger class
 */

class g4rcIO;
class g4rcDetectorConstruction;
class g4rcEventAction;
class g4rcPrimaryGeneratorAction;
class g4rcSteppingAction;

class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithADouble;
class G4UIcmdWithAString;
class G4UIcmdWithABool;
class G4UIdirectory;

class g4rcMessenger : public G4UImessenger {
    public:
       	g4rcMessenger();
       	~g4rcMessenger();

	void SetIO( g4rcIO *io ){ fIO = io; }
	void SetPriGen( g4rcPrimaryGeneratorAction *pg ){ fprigen = pg; }
	void SetDetCon( g4rcDetectorConstruction *dc ){ fdetcon= dc; }
	void SetEvAct( g4rcEventAction *ev ){ fevact = ev; }
	void SetStepAct( g4rcSteppingAction *st ){ fStepAct = st; }

	void SetNewValue(G4UIcommand* cmd, G4String newValue);

    private:
	g4rcIO *fIO;
	g4rcDetectorConstruction *fdetcon;
	g4rcEventAction *fevact;
	g4rcPrimaryGeneratorAction *fprigen;
	g4rcSteppingAction *fStepAct;

	G4UIdirectory *fRemollDir;
	
	G4UIcmdWithAnInteger *seedCmd;
	G4UIcmdWithAString   *fileCmd;
	//
	//G4UIcmdWithAnInteger *fStandModeCmd;
	G4UIcmdWithAnInteger *fDetModeCmd;
	G4UIcmdWithAnInteger *fStandModeCmd;
	G4UIcmdWithAnInteger *fSourceModeCmd;
	G4UIcmdWithADouble *fQuartzPolishCmd;
	G4UIcmdWithADoubleAndUnit *fQuartzThick1Cmd;
	G4UIcmdWithADoubleAndUnit *fQuartzThick2Cmd;
	G4UIcmdWithABool* fFirstDetCmd;
	G4UIcmdWithABool* fSecDetCmd;
	G4UIcmdWithADoubleAndUnit* fSecDetZCmd;
	G4UIcmdWithADoubleAndUnit *fDetAngleCmd;
	// POSSCAN
	G4UIcmdWithADoubleAndUnit *fDetPosXCmd;	
	G4UIcmdWithADoubleAndUnit *fDetPosYCmd;	
	//
	G4UIcmdWithADoubleAndUnit *fXminCmd;
	G4UIcmdWithADoubleAndUnit *fXmaxCmd;
	G4UIcmdWithADoubleAndUnit *fYminCmd;
	G4UIcmdWithADoubleAndUnit *fYmaxCmd;
	G4UIcmdWithADoubleAndUnit *fEminCmd;
	G4UIcmdWithADoubleAndUnit *fEmaxCmd;

	G4UIcmdWithADoubleAndUnit *fthetaMinCmd;
	G4UIcmdWithADoubleAndUnit *fthetaMaxCmd;

//	G4UIcmdWithADoubleAndUnit *fThetaCmd;
//	G4UIcmdWithADoubleAndUnit *fPhiCmd;
	G4UIcmdWithADoubleAndUnit *fZCmd;

};

#endif//g4rcMessenger_HH


