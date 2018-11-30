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

	ConstructMaterials();
	fMaterialManager=this;

}

g4rcMaterial::~g4rcMaterial() {
}


/////////////////////////////////////////////////////////////////////
void g4rcMaterial::ConstructMaterials() {

	// Define elements

	// Molar mass for target gases 
	// Assigned to variables for use in pressure calculation 
	double M_H = 1.01*g/mole;
	double M_D = 2.01*g/mole;
	double M_T = 3.02*g/mole;
	double M_He = 3.02*g/mole; 

	G4Element* H = new G4Element("Hydrogen", "H", 1, M_H);
	G4Element* D = new G4Element("Deuterium", "D", 1, M_D);
	G4Element* T = new G4Element("Tritium", "T", 1, M_T);	
	G4Element* He = new G4Element("Helium", "He", 2, M_He);
	G4Element* Be = new G4Element("Beryllium", "Be", 4, 9.01*g/mole);
	G4Element* C = new G4Element("Carbon", "C", 6, 12.01*g/mole);
	G4Element* N = new G4Element("Nitrogen", "N", 7 , 14.01*g/mole);
	G4Element* O = new G4Element("Oxygen"  , "O", 8 , 16.00*g/mole);
	G4Element* Al = new G4Element("Aluminum", "Al", 13, 26.98*g/mole);

	// Define materials

	// Target gases
	double R = 8.314;

	double rho_H2 = 2.832*mg/cm3;
	double rho_D2 = 5.686*mg/cm3;
	double rho_T2 = 3.404*mg/cm3; 
	double rho_He = 2.135*mg/cm3;
	
	double temp = 45.0*kelvin;

	double P_H2 = rho_H2*(R/M_H)*temp;
	double P_D2 = rho_D2*(R/M_D)*temp;
	double P_T2 = rho_T2*(R/M_T)*temp;
	double P_He = rho_He*(R/M_He)*temp;

	H2_gas = new G4Material("H2_gas", rho_H2, 2, kStateGas, temp, P_H2);
	H2_gas->AddElement(H,1); 
	H2_gas->AddElement(H,1); 

	D2_gas = new G4Material("D2_gas", rho_D2, 2, kStateGas, temp, P_D2);
	D2_gas->AddElement(D,1); 
	D2_gas->AddElement(D,1); 

	T2_gas = new G4Material("T2_gas", rho_T2, 2, kStateGas, temp, P_T2);
	T2_gas->AddElement(T,1); 
	T2_gas->AddElement(T,1); 

	He_gas = new G4Material("He_gas", rho_He, 1, kStateGas, temp, P_He);
	He_gas->AddElement(He,1);

	// Air

	air = new G4Material("Air", 1.29*mg/cm3, 2);
	air->AddElement(N, 70.*perCent);
	air->AddElement(O, 30.*perCent);

	// Aluminum
	aluminum = new G4Material("aluminum", 2.7*g/cm3, 1);
	aluminum->AddElement(Al, 1);
	
	// Beryllium
	beryllium = new G4Material("beryllium", 1.85*g/cm3, 1);
	beryllium->AddElement(Be,1); 

	// Kapton
	kapton = new G4Material("Kapton", 1.42*g/cm3, 4);
	kapton->AddElement(H, 10);
	kapton->AddElement(C, 22);
	kapton->AddElement(O, 5);
	kapton->AddElement(N, 2);

}


