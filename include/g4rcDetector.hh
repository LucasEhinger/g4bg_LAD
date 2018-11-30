#ifndef __QSIMDETECTOR_HH
#define __QSIMDETECTOR_HH

#include "G4VSensitiveDetector.hh"
#include "g4rcDetectorHit.hh"

#include <map>

/*! 
      Default detector class.  This will record information on:

      - Primary generated hit information

*/

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

class g4rcDetector : public G4VSensitiveDetector {
    public:
	g4rcDetector( G4String name, G4int detnum );
	virtual ~g4rcDetector();

	virtual void Initialize(G4HCofThisEvent*);
	virtual G4bool ProcessHits(G4Step*,G4TouchableHistory*);
	virtual void EndOfEvent(G4HCofThisEvent*);

    private:
	g4rcDetectorHitsCollection *fHitColl;
	G4int fHCID;

	G4bool fTrackSecondaries;
	G4int fDetNo;

	G4double fDontRecordThresh;
};

#endif//__QSIMDETECTOR_HH
