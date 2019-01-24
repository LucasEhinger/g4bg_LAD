#include "g4rcCrossSection.hh"

#include <iostream>

#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

extern "C" {
	
	float sigmodel_calc_(float* E0, float* Ep, float* theta, int* Z, int* A, float* M, int* xflag, int* model);

}

g4rcCrossSection::g4rcCrossSection() {


	// Default to hydrogen
	fTarget = "H1";
	Z = 1;
	A = 1;
	M = 0.93828;

	// 1 = QE + DIS
	// 2 = QE only
	// 3 = DIS only
	xf = 1;

}

g4rcCrossSection::~g4rcCrossSection() {

}

G4double g4rcCrossSection::CalculateCrossSection(double pE0, double pEp, double pTh, G4String pModel) {

	int model;

	if(pModel == "bodek") {

		model = 1;

	} else if(pModel == "christy") {

		if(fTarget == "H2") {
			model = 3;
		} else {
			model = 2;
		}

	} else {
		G4cout << "Invalid model specified!" << G4endl;
		exit(1);
	}

	float E0 = pE0/GeV;
	float Ep = pEp/GeV;
	float th = pTh/deg;

	G4cout << "E0 = " << E0 << "\tEp = " << Ep << "\tth = " << th << G4endl;
	
	float mod_xs = sigmodel_calc_(&E0, &Ep, &th, &Z, &A, &M, &xf, &model); 

	return mod_xs;

}

void g4rcCrossSection::SetTarget(G4String pTarget) {
	
	fTarget = pTarget;

	const int nTarg = 4;

	G4String targList[nTarg] = {"H1", "H2", "H3", "He3"};
	int ZList[nTarg] = {1, 1, 1, 2};
	int AList[nTarg] = {1, 2, 3, 3};
	float MList[nTarg] = {0.93828, 1.87537, 2.8095, 2.8094};

	for(int i = 0; i<nTarg; i++) {
		if(fTarget == targList[i]) {
			Z = ZList[i];
			A = AList[i];
			M = MList[i];
			break;
		}
	}

}
