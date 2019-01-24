#ifndef G4RCCROSSSECTION_HH
#define G4RCCROSSSECTION_HH

#include "G4String.hh"

class g4rcCrossSection {

public:

	g4rcCrossSection(G4String);
	~g4rcCrossSection();	


private:

	G4String fModel;


public:

	G4double CalculateCrossSection(G4double, G4double);
		

};

#endif
