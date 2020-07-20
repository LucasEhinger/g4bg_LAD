#ifndef g4rcDetectorConstruction_h
#define g4rcDetectorConstruction_h 

#include "G4VUserDetectorConstruction.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"

class g4rcMaterial;
class G4SDManager;

class g4rcDetectorConstruction : public G4VUserDetectorConstruction {
	
public:
	g4rcDetectorConstruction();
	~g4rcDetectorConstruction();
	
	G4VPhysicalVolume* Construct();

	void SetTarget(G4String);

public:
	G4String fHRS;
	G4double fHRSAngle;	

private:
	G4int fTargIndex;

	G4String fTarg;

	g4rcMaterial* fMaterial;

	void AddGEM(G4LogicalVolume*, int, bool, double, double, G4RotationMatrix*, G4ThreeVector);
	G4double fGEMCenter[1];

	G4SDManager* SDman;

	void GetTargetIndex(G4String);

};

#endif
