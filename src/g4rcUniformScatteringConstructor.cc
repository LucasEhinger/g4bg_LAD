#include "g4rcUniformScatteringConstructor.hh"

#include "G4ProcessManager.hh"
#include "G4Electron.hh"

g4rcUniformScatteringConstructor::g4rcUniformScatteringConstructor() 
 : G4VPhysicsConstructor("UniformScattering") {

	theUniformScatteringProcess = NULL;

}

g4rcUniformScatteringConstructor::~g4rcUniformScatteringConstructor() {
}

void g4rcUniformScatteringConstructor::ConstructParticle() {
	G4Electron::ElectronDefinition();	
}

void g4rcUniformScatteringConstructor::ConstructProcess() {

	theUniformScatteringProcess = new g4rcUniformScattering();

	G4ProcessManager* pManager = G4Electron::Electron()->GetProcessManager();

	pManager->AddDiscreteProcess(theUniformScatteringProcess);
	
}
