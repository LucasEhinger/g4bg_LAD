#ifndef __QSIMDETECTORHIT_HH
#define __QSIMDETECTORHIT_HH

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

class g4rcDetectorHit : public G4VHit {
    public:
	g4rcDetectorHit(G4int, G4int);
	~g4rcDetectorHit();

	g4rcDetectorHit(const g4rcDetectorHit &right);
	const g4rcDetectorHit& operator=(const g4rcDetectorHit &right);
	G4int operator==(const g4rcDetectorHit &right) const;

	inline void *operator new(size_t);
	inline void operator delete(void *aHit);
	void *operator new(size_t,void*p){return p;}

    private:

    public:
	G4int fDetID;
	G4int fCopyID;

	// Position and momentum in lab coordinates
	G4ThreeVector f3X;
	G4ThreeVector f3P;
	// Total momentum, energy, mass
	G4double fP, fE, fM;
	// Origin
	G4ThreeVector f3V;
	G4ThreeVector f3D;
	// Geant4 track ID, particle type, and mother ID
	G4int    fTrID, fPID, fmTrID;
	// Process generator type
	G4int    fGen;
};


typedef G4THitsCollection<g4rcDetectorHit> g4rcDetectorHitsCollection;

extern G4Allocator<g4rcDetectorHit> g4rcDetectorHitAllocator;

inline void* g4rcDetectorHit::operator new(size_t){
    void *aHit;
    aHit = (void *) g4rcDetectorHitAllocator.MallocSingle();
    return aHit;
}

inline void g4rcDetectorHit::operator delete(void *aHit){
    g4rcDetectorHitAllocator.FreeSingle( (g4rcDetectorHit*) aHit);
}

#endif//__QSIMDETECTORHIT_HH
