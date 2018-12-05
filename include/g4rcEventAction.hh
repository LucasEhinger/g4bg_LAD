
#ifndef g4rcEventAction_h
#define g4rcEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class G4Event;
class g4rcIO;
class g4rcSteppingAction;

class g4rcEventAction : public G4UserEventAction
{
  public:
    g4rcEventAction();
    virtual ~g4rcEventAction();

  public:
    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);

	void SetIO( g4rcIO *io ){ fIO = io; }
	void SetSteppingAction(g4rcSteppingAction* step) { fStep = step; }

  private:

	g4rcIO *fIO;
 	g4rcSteppingAction* fStep;

  public:
};

#endif

    
