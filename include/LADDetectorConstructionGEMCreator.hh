#ifndef LADDETECTORCONSTRUCTIONGEMCREATOR_HH
#define LADDETECTORCONSTRUCTIONGEMCREATOR_HH

// Include any necessary headers here
// #include "G4VUserDetectorConstruction.hh"
#include "G4LogicalVolume.hh"
#include "G4VUserDetectorConstruction.hh"
class G4SDManager;
class g4rcMaterial;

class LADDetectorConstructionGEMCreator {
public:
    // Constructor
    LADDetectorConstructionGEMCreator();

    // Destructor
    ~LADDetectorConstructionGEMCreator();

    // Declare any member functions here
    void BuildGEM(G4LogicalVolume *worldLV, G4SDManager* SDman);

private:
    // Declare any member variables here
    g4rcMaterial* fMaterial;

};

#endif // LADDETECTORCONSTRUCTIONGEMCREATOR_HH