#ifndef G4RCUNIFORMSCATTERING_H
#define G4RCUNIFORMSCATTERING_H

#include "G4VDiscreteProcess.hh"
#include "G4ParticleChange.hh"

class g4rcUniformScattering : public G4VDiscreteProcess {

public:

	g4rcUniformScattering(const G4String& processName = "UniformScattering");
	~g4rcUniformScattering();

private:

	g4rcUniformScattering(const g4rcUniformScattering &right);
	g4rcUniformScattering& operator=(const g4rcUniformScattering &right);

	G4double fVertexZ;
	G4double fThetaCentral;

public:

	G4bool IsApplicable(const G4ParticleDefinition&);

	G4double PostStepGetPhysicalInteractionLength(const G4Track&);

	G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*);

	G4ParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

	G4ParticleChange* Scatter(const G4Track&, const G4Step&);

	void SetVertexZ(G4double);
	void SetCentralScatteringAngle(G4double);

};

#endif
