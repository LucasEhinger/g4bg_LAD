#include "g4rcUniformScattering.hh"

#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"

#include "CLHEP/Random/RandFlat.h"

#include <math.h>

g4rcUniformScattering::g4rcUniformScattering(const G4String& processName)
 : G4VDiscreteProcess(processName, fNotDefined) {

	fVertexZ = 7.5*cm;
	fThetaCentral = 17.5*deg;

	fEcut = 1.*MeV;

	fHasScattered = false;
	
}



g4rcUniformScattering::~g4rcUniformScattering() {
}



G4double g4rcUniformScattering::PostStepGetPhysicalInteractionLength(const G4Track& aTrack, G4double previousStepSuze, G4ForceCondition* condition) {

	*condition = NotForced;

	G4double length = DBL_MAX;

	// Only limit step if track is primary

	if(aTrack.GetTrackID()==1) {

		G4double z = aTrack.GetPosition().z();
	
		if(z < fVertexZ) {
			length = fVertexZ - z;
		}
	}

	return length;
}



G4VParticleChange* g4rcUniformScattering::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep) {

	aParticleChange.Initialize(aTrack);

	if(!fHasScattered) {
		G4double Mp = 938.272*MeV;

		G4double Epre = aTrack.GetTotalEnergy();
		G4double Ekin = aTrack.GetKineticEnergy();

		// Choose theta
		G4double fThetaMin = fThetaCentral - 10.*deg;
		G4double fThetaMax = fThetaCentral + 10.*deg;	
		G4double fCosThMin = cos(fThetaMax);
		G4double fCosThMax = cos(fThetaMin);
		G4double theta = acos(CLHEP::RandFlat::shoot(fCosThMin, fCosThMax));

		// Set upper Q2 limit based on angle, and choose Q2
		G4double fQ2Min = 0.;
		G4double fQ2Max = 2.*Epre*Epre*(1. - cos(theta));
		G4double Q2_true = CLHEP::RandFlat::shoot(fQ2Min, fQ2Max);

		// Choose phi;
		G4double fPhiMin = -20.*deg;
		G4double fPhiMax = +20.*deg;
		G4double phi = CLHEP::RandFlat::shoot(fPhiMin, fPhiMax);

		G4double internal_loss1;
		G4double E0, Ef, nu_true;

		nu_true = -1.;

		while(nu_true < 0.) {
			internal_loss1 = RadiateInternal(Q2_true, Ekin);
			E0 = Epre - internal_loss1;
			Ef = Q2_true/(2.*E0*(1. - cos(theta)));
			nu_true = E0 - Ef;
		}

		G4double x_true = Q2_true/(2.*Mp*nu_true);

		Ekin = Ekin - internal_loss1 - nu_true;
		G4double internal_loss2 = RadiateInternal(Q2_true, Ekin);

		G4double Epost = Ef - internal_loss2;
		G4ThreeVector p = G4ThreeVector(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta));

		fEpre = Epre;
		fE0 = E0;
		fEp = Ef;
		fEpost = Epost;
		fTheta = theta/deg;
		fQ2true = Q2_true;
		fxBtrue = x_true;

		aParticleChange.ProposeEnergy(Epost);
		aParticleChange.ProposeMomentumDirection(p);

		fHasScattered = true;
	}
	return &aParticleChange;			

}


G4double g4rcUniformScattering::RadiateInternal(G4double Q2, G4double E) {

	G4double Euler = 0.5772157;
	G4double alpha = 1./137.;

	// Equivalent radiator thickness from Mo and Tsai
	G4double bt = (alpha/pi)*(log(Q2/(electron_mass_c2*electron_mass_c2)) - 1.);

	G4double eps = fEcut/E;

	G4double prob, prob_sample, sample, Eloss, value;

	value = 1.;
	
	prob = 1. - pow(eps, bt) - (bt/(bt+1.))*(1. - pow(eps, bt+1.)) + (3./4.)*(bt/(bt+2.))*(1. - pow(eps, bt+2.));
	prob/=(1. - bt*Euler + ((bt*bt)/12.)*(6.*Euler*Euler + pi*pi));

	prob_sample = G4UniformRand();

	G4double env, ref;

	Eloss = 0.;

	if (prob_sample <= prob) {
		do {
			sample = G4UniformRand();
			Eloss = fEcut*pow(E/fEcut, sample);
			env = 1./Eloss;
			value = (1./Eloss)*pow(Eloss/E, bt)*(1. - (Eloss/E) + (3./4.)*pow(Eloss/E, 2.));

			sample = G4UniformRand();
			ref = value/env;
		} while (sample > ref);
	}

	return Eloss;

}



void g4rcUniformScattering::SetVertexZ(G4double z) {
	fVertexZ = z;
}



void g4rcUniformScattering::SetCentralScatteringAngle(G4double th) {
	fThetaCentral = th;
}

