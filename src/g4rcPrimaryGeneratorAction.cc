#include "g4rcPrimaryGeneratorAction.hh"


#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "TFile.h"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "g4rcIO.hh"
#include "g4rcEvent.hh"
#include "g4rctypes.hh"
#include "globals.hh"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"

#include "G4SystemOfUnits.hh"

g4rcPrimaryGeneratorAction::g4rcPrimaryGeneratorAction() {

  	G4int n_particle = 1;

	fXmin = -1.*mm;
	fXmax =  1.*mm;
	fYmin = -1.*mm;
	fYmax =  1.*mm;
	
	fZ = -35.*cm;

	fEbeam = 10.9*GeV;

  	fParticleGun = new G4ParticleGun(n_particle);
  	fDefaultEvent = new g4rcEvent();

}


g4rcPrimaryGeneratorAction::~g4rcPrimaryGeneratorAction() {

	delete fParticleGun;
	delete fDefaultEvent;

}


void g4rcPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {

    /*  Generate event, set IO data */

    // Use default, static single generator
    // Update this just in case things changed
    // from the command user interface
    fDefaultEvent->Reset();

    // Set data //////////////////////////////////
    // Magic happens here

	double xPos, yPos, zPos;	

	xPos = CLHEP::RandFlat::shoot(fXmin, fXmax );
	yPos = CLHEP::RandFlat::shoot(fYmin, fYmax );
	
	zPos = fZ;

	double mass = 0.511*MeV;
	double p = sqrt(fEbeam*fEbeam - mass*mass);

	fDefaultEvent->ProduceNewParticle(G4ThreeVector(xPos, yPos, zPos), G4ThreeVector(0., 0., p), "e-");

    	// Register and create event
	double kinE = sqrt(fDefaultEvent->fPartMom[0].mag()*fDefaultEvent->fPartMom[0].mag()
	    + fDefaultEvent->fPartType[0]->GetPDGMass()*fDefaultEvent->fPartType[0]->GetPDGMass() )
	-fDefaultEvent->fPartType[0]->GetPDGMass();

	fParticleGun->SetParticleDefinition(fDefaultEvent->fPartType[0]);
	fParticleGun->SetParticleMomentumDirection(fDefaultEvent->fPartMom[0].unit());
	fParticleGun->SetParticleEnergy( kinE  );
	fParticleGun->SetParticlePosition( fDefaultEvent->fPartPos[0] );


	fIO->SetEventData(fDefaultEvent);
	fParticleGun->GeneratePrimaryVertex(anEvent);

}

G4ParticleGun* g4rcPrimaryGeneratorAction::GetParticleGun() {
  return fParticleGun;
} 

