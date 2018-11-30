#ifndef g4rcDetectorConstruction_h
#define g4rcDetectorConstruction_h 

#include "G4VUserDetectorConstruction.hh"

class g4rcMaterial;

class g4rcDetectorConstruction : public G4VUserDetectorConstruction {
	
public:
	g4rcDetectorConstruction();
	~g4rcDetectorConstruction();
	
	G4VPhysicalVolume* Construct();


public:
	G4String fTarg, fHRS;	

private:
	G4double berillium_window;
	G4double chamber_window;
	G4double q1_window;
	// cell parameters
	G4double length;
	G4double radius;
	G4double entrance_window;
	G4double entrance_wall;
	G4double mid_wall;
	G4double exit_wall;
	G4double exit_window;

	g4rcMaterial* fMaterial;

};

#endif
