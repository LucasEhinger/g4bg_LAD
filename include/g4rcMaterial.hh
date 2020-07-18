#ifndef g4rcMaterial_h
#define g4rcMaterial_h  

#include <G4SystemOfUnits.hh>
#include <G4PhysicalConstants.hh>

#include "G4NistManager.hh"
#include "G4Isotope.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4MaterialTable.hh"

class g4rcMaterial {

public:
	g4rcMaterial();
	virtual ~g4rcMaterial();

	//Static method which returns the singleton pointer of this class.
	static g4rcMaterial* GetMaterialManager();

private:
	static g4rcMaterial* fMaterialManager;

	G4NistManager* fNistMan;

private:
	void ConstructMaterials();

public:

	G4Material* vacuum;
	G4Material* H2_gas;
	G4Material* D2_gas;
	G4Material* T2_gas;
	G4Material* He_gas;
	G4Material* air;
	G4Material* aluminum;
	G4Material* beryllium;
	G4Material* kapton;
	G4Material* H2_liquid;
	G4Material* D2_liquid;


	G4Material* Galaxy;
	G4Material* NemaG10;
	G4Material* CO2;
	G4Material* ArCO2;
	G4Material* Kapton;
	G4Material* Kapton0d2;
	G4Material* Kapton0d8;
	G4Material* Copper;
	G4Material* Copper0d2;
	G4Material* Copper0d75;
	G4Material* Copper0d8;



};
typedef g4rcMaterial g4rcMaterialManager;

#endif

