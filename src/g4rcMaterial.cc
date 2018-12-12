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

	double Z, A;
	G4Element* H = new G4Element("Hydrogen","H", 	Z = 1, 	A = M_H);
	G4Element* D = new G4Element("Deuterium","D", 	Z = 1, 	A = M_D);
	G4Element* T = new G4Element("Tritium","T", 	Z = 1, 	A = M_T);	
	G4Element* He = new G4Element("Helium","He", 	Z = 2, 	A = M_He);
	G4Element* Be = new G4Element("Beryllium","Be",	Z = 4, 	A = 9.01*g/mole);
	G4Element* C = new G4Element("Carbon", "C", 	Z = 6, 	A = 12.01*g/mole);
	G4Element* N = new G4Element("Nitrogen", "N", 	Z = 7, 	A = 14.01*g/mole);
	G4Element* O = new G4Element("Oxygen", "O", 	Z = 8, 	A = 16.00*g/mole);
	G4Element* Al = new G4Element("Aluminum", "Al",	Z = 13,	A = 26.98*g/mole);

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

	double density, nComp, pressure; 
	H2_gas = new G4Material("H2_gas", density = rho_H2, nComp = 1, kStateGas, temp, pressure = P_H2);
	H2_gas->AddElement(H,2); 

	D2_gas = new G4Material("D2_gas", density = rho_D2, nComp = 1, kStateGas, temp, pressure = P_D2);
	D2_gas->AddElement(D,2); 

	T2_gas = new G4Material("T2_gas", density = rho_T2, nComp = 1, kStateGas, temp, pressure = P_T2);
	T2_gas->AddElement(T,2); 

	He_gas = new G4Material("He_gas", density = rho_He, nComp = 1, kStateGas, temp, pressure = P_He);
	He_gas->AddElement(He,1);


	//A pressure that causes the Hg column to rise 1 millimeter is called ONE torr (1mm Hg)
	// Vacuum of 1.e-6 torr at room temperature,  1 atmosphere = 101325*pascal = 760 *torr
	double rho_vacuum = 1.e-6/760.0 * 1.29*mg/cm3; //0.001 of air density
	double P_vacuum = 1.e-6/760.0 *atmosphere;
	vacuum = new G4Material("vacuum", 1, 28.97*g/mole, rho_vacuum, kStateGas, 293.15*kelvin, P_vacuum);

	// Air
	air = new G4Material("air", density = 1.225*mg/cm3, nComp = 2, kStateGas, CLHEP::STP_Temperature, CLHEP::STP_Pressure);
	air->AddElement(N, 80.*perCent);
	air->AddElement(O, 20.*perCent);

	// Aluminum
	aluminum = new G4Material("aluminum", density = 2.7*g/cm3, nComp = 1, kStateSolid, CLHEP::STP_Temperature, CLHEP::STP_Pressure);
	aluminum->AddElement(Al, 1);
	
	// Beryllium
	beryllium = new G4Material("beryllium", density = 1.85*g/cm3, nComp = 1, kStateSolid, CLHEP::STP_Temperature, CLHEP::STP_Pressure);
	beryllium->AddElement(Be,1); 

	// Kapton
	kapton = new G4Material("kapton", density = 1.42*g/cm3, nComp = 4, kStateSolid, CLHEP::STP_Temperature, CLHEP::STP_Pressure);
	kapton->AddElement(H, 10);
	kapton->AddElement(C, 22);
	kapton->AddElement(N, 2);
	kapton->AddElement(O, 5);

}


