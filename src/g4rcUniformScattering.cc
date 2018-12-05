#include "g4rcUniformScattering.hh"

#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"

#include "CLHEP/Random/RandFlat.h"

#include <math.h>

g4rcUniformScattering::g4rcUniformScattering(const G4String& processName)
 : G4VDiscreteProcess(processName) {

	fVertexZ = 0.;
	fThetaCentral = 17.5*deg;
	
}



g4rcUniformScattering::~g4rcUniformScattering() {
}



G4bool g4rcUniformScattering::IsApplicable(const G4ParticleDefinition& aParticle) {

//	return (aParticle.GetParticleName() == "electron");
	return true;

}



G4double g4rcUniformScattering::GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*) {

	return DBL_MIN;

}



G4double g4rcUniformScattering::PostStepGetPhysicalInteractionLength(const G4Track& aTrack) {

	G4double length = DBL_MAX;

	G4double z = aTrack.GetPosition().z();

	if (z < fVertexZ) {

		length = fVertexZ - z;

	}

	return length;

}



G4ParticleChange* g4rcUniformScattering::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep) {

	G4double z = aTrack.GetPosition().z();

	G4ParticleChange* theParticleChange;
	theParticleChange->Initialize(aTrack);

	if (fabs(fVertexZ - z) < 1.e-2) {
		theParticleChange = Scatter(aTrack, aStep);
	}

	return theParticleChange;

}



G4ParticleChange* g4rcUniformScattering::Scatter(const G4Track& aTrack, const G4Step& aStep) {

	G4ParticleChange fParticleChange;

	fParticleChange.Initialize(aTrack);
	
	G4double fThetaMin = fThetaCentral - 10.*deg;
	G4double fThetaMax = fThetaCentral + 10.*deg;

	G4double fPhiMin = -20.*deg;
	G4double fPhiMax = +20.*deg;

	G4double fCosThMin = cos(fThetaMax);
	G4double fCosThMax = cos(fThetaMin);

	G4double fEpMin = 2.6*GeV;
	G4double fEpMax = 3.6*GeV;

	G4double theta = acos(CLHEP::RandFlat::shoot(fCosThMin, fCosThMax));
	G4double phi = CLHEP::RandFlat::shoot(fPhiMin, fPhiMax);
	G4double Eprime = CLHEP::RandFlat::shoot(fEpMin, fEpMax);

	fParticleChange.ProposeEnergy(Eprime);

	G4ThreeVector p = G4ThreeVector(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta));

	fParticleChange.ProposeMomentumDirection(p);
	
	return &fParticleChange;			

}



void g4rcUniformScattering::SetVertexZ(G4double z) {
	fVertexZ = z;
}



void g4rcUniformScattering::SetCentralScatteringAngle(G4double th) {
	fThetaCentral = th;
}



















