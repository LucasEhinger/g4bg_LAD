#ifndef __QSIMSCINTDETECTORHIT_HH
#define __QSIMSCINTDETECTORHIT_HH

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

class g4rcScintDetectorHit : public G4VHit {
    public:
	g4rcScintDetectorHit(G4int, G4int);
	~g4rcScintDetectorHit();

	g4rcScintDetectorHit(const g4rcScintDetectorHit &right);
	const g4rcScintDetectorHit& operator=(const g4rcScintDetectorHit &right);
	G4int operator==(const g4rcScintDetectorHit &right) const;

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

	G4double fEdep;
};


typedef G4THitsCollection<g4rcScintDetectorHit> g4rcScintDetectorHitsCollection;

extern G4Allocator<g4rcScintDetectorHit> g4rcScintDetectorHitAllocator;

inline void* g4rcScintDetectorHit::operator new(size_t){
    void *aHit;
    aHit = (void *) g4rcScintDetectorHitAllocator.MallocSingle();
    return aHit;
}

inline void g4rcScintDetectorHit::operator delete(void *aHit){
    g4rcScintDetectorHitAllocator.FreeSingle( (g4rcScintDetectorHit*) aHit);
}

#endif//__QSIMSCINTDETECTORHIT_HH
