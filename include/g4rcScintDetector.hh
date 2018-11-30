#ifndef __QSIMSCINTDETECTOR_HH
#define __QSIMSCINTDETECTOR_HH

#include "G4VSensitiveDetector.hh"
#include "g4rcScintDetectorHit.hh"

#include <map>

/*! 
      Default detector class.  This will record information on:

      - Primary generated hit information

*/

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

class g4rcScintDetector : public G4VSensitiveDetector {
    public:
	g4rcScintDetector( G4String name, G4int detnum );
	virtual ~g4rcScintDetector();

	virtual void Initialize(G4HCofThisEvent*);
	virtual G4bool ProcessHits(G4Step*,G4TouchableHistory*);
	virtual void EndOfEvent(G4HCofThisEvent*);

    private:
	g4rcScintDetectorHitsCollection *fHitColl;
	G4int fHCID;

	G4bool fTrackSecondaries;
	G4int fDetNo;

	std::map<int, g4rcScintDetectorHit *> fSumMap;

	G4double fDontRecordThresh;
};

#endif//__QSIMSCINTDETECTOR_HH
