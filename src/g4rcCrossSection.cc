#include "g4rcCrossSection.hh"

#include <iostream>

#include "G4SystemOfUnits.hh"

extern "C" {
	
	float sigmodel_calc_(float* E0, float* Ep, float* theta, int* Z, int* A, float* M, int* xflag, int* model);

}

g4rcCrossSection::g4rcCrossSection(G4String model) {

	fModel = model;

}

g4rcCrossSection::~g4rcCrossSection() {

}

G4double g4rcCrossSection::CalculateCrossSection(G4double theta, G4double Eprime) {


	float E0 = 10.6;
	float Ep = Eprime/GeV;
	float th = theta/deg;
	int Z = 1;
	int A = 2;
	float M = 1.809;
	int xf = 3;
	int model = 1;

	float mod_xs = sigmodel_calc_(&E0, &Ep, &th, &Z, &A, &M, &xf, &model); 
	std::cout << "GEANT4 Ep = " << Ep << "\n\n";	
	return mod_xs;

}

