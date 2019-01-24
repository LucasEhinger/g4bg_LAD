#ifndef G4RCCROSSSECTION_HH
#define G4RCCROSSSECTION_HH

#include "G4String.hh"

class g4rcCrossSection {

public:

	g4rcCrossSection();
	~g4rcCrossSection();	


private:

	int Z;
	int A;
	float M;
	int xf;
	

	G4String fTarget;


public:

	G4double CalculateCrossSection(double, double, double, G4String);
	void SetTarget(G4String);

};

#endif
