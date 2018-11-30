#include "g4rcSteppingAction.hh"

#include "G4VVisManager.hh"
#include "G4Polyline.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4SteppingManager.hh"

g4rcSteppingAction::g4rcSteppingAction()
:drawFlag(false)
{

}

void g4rcSteppingAction::UserSteppingAction(const G4Step *aStep) {
}


