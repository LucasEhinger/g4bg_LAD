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
#include "G4RunManager.hh"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"

#include "G4SystemOfUnits.hh"

g4rcPrimaryGeneratorAction::g4rcPrimaryGeneratorAction()
{

	G4int n_particle = 1;

	fXmin = -1. * mm;
	fXmax = 1. * mm;
	fYmin = -1. * mm;
	fYmax = 1. * mm;

	fZ = -35. * cm;

	fEbeam = 10.9 * GeV;

	fThetaMin = 95. * deg;
	fThetaMax = 160. * deg;
	fPhiMin = -17. * deg;
	fPhiMax = 17. * deg;

	radius_start = 0. * cm;

	fPprotonMax = 1 * GeV;

	fParticleGun = new G4ParticleGun(n_particle);
	fDefaultEvent = new g4rcEvent();
	numberOfEvents = 0;
}

g4rcPrimaryGeneratorAction::~g4rcPrimaryGeneratorAction()
{

	delete fParticleGun;
	delete fDefaultEvent;
}

void g4rcPrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent)
{

	/*  Generate event, set IO data */

	// Use default, static single generator
	// Update this just in case things changed
	// from the command user interface
	fDefaultEvent->Reset();

	// Set data //////////////////////////////////
	// Magic happens here

	double xPos, yPos, zPos;
	double px, py, pz;

	if (fProtBool == 0)
	{
		xPos = CLHEP::RandFlat::shoot(fXmin, fXmax);
		yPos = CLHEP::RandFlat::shoot(fYmin, fYmax);

		zPos = fZ;

		double mass = 0.511 * MeV;
		double p = sqrt(fEbeam * fEbeam - mass * mass);

		fDefaultEvent->ProduceNewParticle(G4ThreeVector(xPos, yPos, zPos), G4ThreeVector(0., 0., p), "e-");
	}
	else
	{
		int eventID = anEvent->GetEventID();
		// double theta_p = CLHEP::RandFlat::shoot(fThetaMin, fThetaMax);
		// double phi_p = CLHEP::RandFlat::shoot(fPhiMin, fPhiMax);
		// theta_p=127.*deg;
		numberOfEvents = G4RunManager::GetRunManager()->GetNumberOfEventsToBeProcessed();
		int n_floor = floor(sqrt(numberOfEvents));
		int n_ceil = ceil(sqrt(numberOfEvents));
		int n_mod = eventID % n_floor;
		int n_div = eventID / n_floor;
		double theta_p = fThetaMin + ((double)n_mod) / n_floor * (fThetaMax - fThetaMin);
		double phi_p = 0. * deg;
		// double fPproton = CLHEP::RandFlat::shoot(0, fPprotonMax);
		double fPproton = ((double)n_div) / n_ceil * fPprotonMax;

		px = fPproton * sin(theta_p) * cos(phi_p);
		py = fPproton * sin(theta_p) * sin(phi_p);
		pz = fPproton * cos(theta_p);

		xPos = radius_start * sin(theta_p) * cos(phi_p);
		yPos = radius_start * sin(theta_p) * sin(phi_p);
		zPos = radius_start * cos(theta_p);

		fDefaultEvent->ProduceNewParticle(G4ThreeVector(xPos, yPos, zPos), G4ThreeVector(px, py, pz), "proton");
	}

	// Register and create event
	double kinE = sqrt(fDefaultEvent->fPartMom[0].mag() * fDefaultEvent->fPartMom[0].mag() + fDefaultEvent->fPartType[0]->GetPDGMass() * fDefaultEvent->fPartType[0]->GetPDGMass()) - fDefaultEvent->fPartType[0]->GetPDGMass();

	fParticleGun->SetParticleDefinition(fDefaultEvent->fPartType[0]);
	fParticleGun->SetParticleMomentumDirection(fDefaultEvent->fPartMom[0].unit());
	fParticleGun->SetParticleEnergy(kinE);
	fParticleGun->SetParticlePosition(fDefaultEvent->fPartPos[0]);

	fIO->SetEventData(fDefaultEvent);
	fParticleGun->GeneratePrimaryVertex(anEvent);
}

G4ParticleGun *g4rcPrimaryGeneratorAction::GetParticleGun()
{
	return fParticleGun;
}
