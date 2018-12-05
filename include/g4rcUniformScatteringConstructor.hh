#ifndef G4RCUNIFORMSCATTERINGCONSTRUCTOR_HH
#define G4RCUNIFORMSCATTERINGCONSTRUCTOR_HH

#include "g4rcUniformScattering.hh"

#include "G4VPhysicsConstructor.hh"


class g4rcUniformScatteringConstructor : public G4VPhysicsConstructor {

public: 
	
	g4rcUniformScatteringConstructor();
	virtual ~g4rcUniformScatteringConstructor();

	virtual void ConstructParticle();
	virtual void ConstructProcess();
	
	g4rcUniformScattering* GetUniformScatteringProcess() {return theUniformScatteringProcess;}


private:
	
	g4rcUniformScattering* theUniformScatteringProcess;

};


#endif























