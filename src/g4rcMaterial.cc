#include <stdio.h>
#include <math.h>
#include <string>
using namespace std;

#include "g4rcMaterial.hh"


g4rcMaterial* g4rcMaterial::fMaterialManager=0;
	
g4rcMaterial* g4rcMaterial::GetMaterialManager() {
 
	if(!fMaterialManager) 
	{
            fMaterialManager = new g4rcMaterial();
	}
	return fMaterialManager; 
}

/////////////////////////////////////////////////////////////////////

g4rcMaterial::g4rcMaterial() {

	fNistMan = G4NistManager::Instance();
	ConstructMaterials();
	fMaterialManager=this;

}

g4rcMaterial::~g4rcMaterial() {
}


/////////////////////////////////////////////////////////////////////
void g4rcMaterial::ConstructMaterials() {

	// Define elements

	G4int z, n;
	G4double a;
	G4double density;
	G4int ncomponents, natoms;
	G4double fractionmass;


	G4double fExtDensityRatio = 1.0;

	// Define elements from NIST material table
	G4Element *H  = fNistMan->FindOrBuildElement(z = 1);
	G4Element *He = fNistMan->FindOrBuildElement(z = 2);
	G4Element *Be = fNistMan->FindOrBuildElement(z = 4);
	G4Element *C  = fNistMan->FindOrBuildElement(z = 6);
	G4Element *N  = fNistMan->FindOrBuildElement(z = 7);
	G4Element *O  = fNistMan->FindOrBuildElement(z = 8);
	G4Element *F  = fNistMan->FindOrBuildElement(z = 9);
	G4Element *Na = fNistMan->FindOrBuildElement(z = 11);
	G4Element *Al = fNistMan->FindOrBuildElement(z = 13);
	G4Element *Si = fNistMan->FindOrBuildElement(z = 14);
	G4Element *P  = fNistMan->FindOrBuildElement(z = 15);
	G4Element *S  = fNistMan->FindOrBuildElement(z = 16);
	G4Element *Ar = fNistMan->FindOrBuildElement(z = 18);
	G4Element *K  = fNistMan->FindOrBuildElement(z = 19);
	G4Element *Cr = fNistMan->FindOrBuildElement(z = 24);
	G4Element *Mn = fNistMan->FindOrBuildElement(z = 25);
	G4Element *Fe = fNistMan->FindOrBuildElement(z = 26);
	G4Element *Ni = fNistMan->FindOrBuildElement(z = 28);
	G4Element *Cu = fNistMan->FindOrBuildElement(z = 29);
	G4Element *Zn = fNistMan->FindOrBuildElement(z = 30);
	G4Element *As = fNistMan->FindOrBuildElement(z = 33);
	G4Element *W  = fNistMan->FindOrBuildElement(z = 74);
	G4Element *Pb = fNistMan->FindOrBuildElement(z = 82);


	// Space Vacuum
	Galaxy = new G4Material("Galaxy", density = universe_mean_density, ncomponents = 1, kStateGas, 0.1 * kelvin, 1.0e-19 * pascal);
	Galaxy->AddElement(H, fractionmass = 1.0);
//	fVisAtts[Galaxy->GetName()] = new G4VisAttributes(G4VisAttributes::Invisible);

	// GEM Frame G10
	NemaG10 = new G4Material("NemaG10", density = fExtDensityRatio * 1.700 * g / cm3, ncomponents = 4);
	NemaG10->AddElement(Si, natoms = 1);
	NemaG10->AddElement(O, natoms = 2);
	NemaG10->AddElement(C, natoms = 3);
	NemaG10->AddElement(H, natoms = 3);
//	fVisAtts[NemaG10->GetName()] = new G4VisAttributes(G4Colour::Brown());

	// Ar/CO2 Gas
	CO2 = new G4Material("CO2", density = fExtDensityRatio * 1.842e-3 * g / cm3, ncomponents = 2);
	CO2->AddElement(C, natoms = 1);
	CO2->AddElement(O, natoms = 2);
	ArCO2 = new G4Material("ArCO2", density = fExtDensityRatio * 1.715e-3 * g / cm3, ncomponents = 2);
	ArCO2->AddElement(Ar, fractionmass = 0.7);
	ArCO2->AddMaterial(CO2, fractionmass = 0.3);
//	fVisAtts[ArCO2->GetName()] = new G4VisAttributes(G4Colour::Yellow());

	// Kapton
	Kapton = new G4Material("Kapton", density = fExtDensityRatio * 1.42 * g / cm3, ncomponents = 4);
	Kapton->AddElement(H, fractionmass = 0.0273);
	Kapton->AddElement(C, fractionmass = 0.7213);
	Kapton->AddElement(N, fractionmass = 0.0765);
	Kapton->AddElement(O, fractionmass = 0.1749);
//	fVisAtts[Kapton->GetName()] = new G4VisAttributes(G4Colour::Brown());
	Kapton0d2 = new G4Material("Kapton0.2", Kapton->GetDensity() * 0.2, Kapton);
//	fVisAtts[Kapton0d2->GetName()] = new G4VisAttributes(G4Colour::Brown());
	Kapton0d8 = new G4Material("Kapton0.8", Kapton->GetDensity() * 0.8, Kapton);
//	fVisAtts[Kapton0d8->GetName()] = new G4VisAttributes(G4Colour::Brown());

	// Copper C101
	Copper = new G4Material("Copper", density = fExtDensityRatio * 8.92 * g / cm3, ncomponents = 1);
	Copper->AddElement(Cu, natoms = 1);
//	fVisAtts[Copper->GetName()] = new G4VisAttributes(G4Colour::Brown());
	Copper0d2 = new G4Material("Copper0.2", Copper->GetDensity() * 0.2, Copper);
//	fVisAtts[Copper0d2->GetName()] = new G4VisAttributes(G4Colour::Brown());
	Copper0d75 = new G4Material("Copper0.75", Copper->GetDensity() * 0.75, Copper);
//	fVisAtts[Copper0d75->GetName()] = new G4VisAttributes(G4Colour::Brown());
	Copper0d8 = new G4Material("Copper0.8", Copper->GetDensity() * 0.8, Copper);
//	fVisAtts[Copper0d8->GetName()] = new G4VisAttributes(G4Colour::Brown());


	H2_liquid = fNistMan->FindOrBuildMaterial("G4_lH2");

	poly = fNistMan->FindOrBuildMaterial("G4_POLYETHYLENE");
	
	// ----  gas and liquid deuterium

	// Deuteron isotope
	G4Isotope* deuteron  = new G4Isotope("deuteron", 1, 2, 2.0141018*g/mole);
	
	// Deuterium element
	G4Element* deuterium = new G4Element("deuterium", "deuterium", 1);
	deuterium->AddIsotope(deuteron, 1);
	
	// Liquid Deuterium
	D2_liquid = new G4Material("LD2", 0.174*g/cm3, 1, kStateLiquid, 19.0*kelvin);
	D2_liquid->AddElement(deuterium, 2);


	//A pressure that causes the Hg column to rise 1 millimeter is called ONE torr (1mm Hg)
	// Vacuum of 1.e-6 torr at room temperature,  1 atmosphere = 101325*pascal = 760 *torr
	double rho_vacuum = 1.e-6/760.0 * 1.29*mg/cm3; //0.001 of air density
	double P_vacuum = 1.e-6/760.0 *atmosphere;
	vacuum = new G4Material("vacuum", 1, 28.97*g/mole, rho_vacuum, kStateGas, 293.15*kelvin, P_vacuum);

	// Air
	air = new G4Material("air", density = 1.225*mg/cm3, ncomponents = 2, kStateGas, CLHEP::STP_Temperature, CLHEP::STP_Pressure);
	air->AddElement(N, 80.*perCent);
	air->AddElement(O, 20.*perCent);

	// Aluminum
	aluminum = new G4Material("aluminum", density = 2.7*g/cm3, ncomponents = 1, kStateSolid, CLHEP::STP_Temperature, CLHEP::STP_Pressure);
	aluminum->AddElement(Al, 1);
	
	// Beryllium
	beryllium = new G4Material("beryllium", density = 1.85*g/cm3, ncomponents = 1, kStateSolid, CLHEP::STP_Temperature, CLHEP::STP_Pressure);
	beryllium->AddElement(Be,1); 

	// Kapton
	kapton = new G4Material("kapton", density = 1.42*g/cm3, ncomponents = 4, kStateSolid, CLHEP::STP_Temperature, CLHEP::STP_Pressure);
	kapton->AddElement(H, 10);
	kapton->AddElement(C, 22);
	kapton->AddElement(N, 2);
	kapton->AddElement(O, 5);

}


