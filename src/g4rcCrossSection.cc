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


	float E0 = 10.589;
	float Ep = Eprime/GeV;
	float th = theta/deg;
	int Z = 1;
	int A = 1;
	float M = 0.93828;
	int xf = 1;
	int model = 1;

	float mod_xs = sigmodel_calc_(&E0, &Ep, &th, &Z, &A, &M, &xf, &model); 

	return mod_xs;

}

