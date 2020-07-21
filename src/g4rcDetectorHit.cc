#include "g4rcDetectorHit.hh"

G4Allocator<g4rcDetectorHit> g4rcDetectorHitAllocator;

g4rcDetectorHit::g4rcDetectorHit(G4int det, G4int copy){
    fDetID  = det;
    fCopyID = copy;

    f3X = G4ThreeVector(-1e9, -1e9, -1e9);
    f3P = G4ThreeVector(-1e9, -1e9, -1e9);
    f3V = G4ThreeVector(-1e9, -1e9, -1e9);
    f3D = G4ThreeVector(-1e9, -1e9, -1e9);

    fEdep = 0.0;	

    fP  = -1.0;
    fE  = -1.0;
    fM  = -1.0;

    fTrID  = -1;
    fPID   = (G4int) 1e9;
    fmTrID = -1;

    fGen   = 1;
}

g4rcDetectorHit::~g4rcDetectorHit(){
}

g4rcDetectorHit::g4rcDetectorHit(const g4rcDetectorHit &right) : G4VHit(){
    // copy constructor

    fDetID  = right.fDetID;
    fCopyID = right.fCopyID;
    f3X     = right.f3X;
    f3P     = right.f3P;
    f3D     = right.f3D;
    f3V     = right.f3V;

    fEdep   = right.fEdep;

    fP      = right.fP;
    fE      = right.fE;
    fM      = right.fM;

    fTrID   = right.fTrID;
    fPID    = right.fPID;
    fmTrID  = right.fmTrID;
    fGen    = right.fGen;
}

const g4rcDetectorHit& g4rcDetectorHit::operator =(const g4rcDetectorHit &right){
    (*this) = right;
    return *this;
}

G4int g4rcDetectorHit::operator==(const g4rcDetectorHit &right ) const {
    return (this==&right) ? 1 : 0;
}

void g4rcDetectorHit::AddEdep(G4double thisEdep) {
	fEdep += thisEdep;
}

G4double g4rcDetectorHit::GetEdep() {
	return fEdep;
}

void g4rcDetectorHit::SetTrackID(G4int thisTrackID) {
	fTrackID = thisTrackID;
}

G4int g4rcDetectorHit::GetTrackID() {
	return fTrackID;
}
