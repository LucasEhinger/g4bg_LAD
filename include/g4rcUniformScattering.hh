#ifndef G4RCUNIFORMSCATTERING_H
#define G4RCUNIFORMSCATTERING_H

#include "G4VDiscreteProcess.hh"

class g4rcUniformScattering : public G4VDiscreteProcess {

public:

	g4rcUniformScattering(const G4String& processName = "UniformScattering");
	virtual ~g4rcUniformScattering();

private:

	g4rcUniformScattering(const g4rcUniformScattering &right);
	g4rcUniformScattering& operator=(const g4rcUniformScattering &right);

	G4double fVertexZ;
	G4bool fEcut;
	G4double fThetaCentral;

	G4double RadiateInternal(G4double, G4double);

public:

	virtual G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*) {return -1.;}
	virtual G4double PostStepGetPhysicalInteractionLength(const G4Track&, G4double, G4ForceCondition*);
	virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

	void SetVertexZ(G4double);
	void SetCentralScatteringAngle(G4double);

	G4bool fHasScattered;
	G4double fEpre;
	G4double fE0;
	G4double fEp;
	G4double fEpost;
	G4double fTheta;
	G4double fPhi;
	G4double fQ2true;
	G4double fxBtrue;
	G4double fThTarg;
	G4double fPhTarg;

};

#endif
