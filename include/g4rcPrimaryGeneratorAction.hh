
#ifndef g4rcPrimaryGeneratorAction_h
#define g4rcPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4String.hh"
#include "THn.h"

class G4ParticleGun;
class G4Event;
class g4rcIO;
class g4rcEvent;

class g4rcPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {

public:
	g4rcPrimaryGeneratorAction();
	~g4rcPrimaryGeneratorAction();

public:
	void GeneratePrimaries(G4Event* anEvent);
	G4ParticleGun* GetParticleGun();
	void SetIO( g4rcIO *io ){ fIO = io; }

private:
	G4ParticleGun* fParticleGun;
	G4int fSourceMode;

	g4rcEvent *fDefaultEvent;
	g4rcIO *fIO;

public:
	G4double fXmin, fXmax, fYmin, fYmax;
	G4double fZ;
	G4double fEbeam;

};

#endif


