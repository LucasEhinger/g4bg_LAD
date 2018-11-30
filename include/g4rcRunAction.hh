
#ifndef g4rcRunAction_h
#define g4rcRunAction_h 1

#include "globals.hh"
#include "G4UserRunAction.hh"

class G4Timer;
class G4Run;
class g4rcIO;

class g4rcRunAction : public G4UserRunAction
{
  public:
    g4rcRunAction();
    ~g4rcRunAction();

  public:
    void BeginOfRunAction(const G4Run* aRun);
    void EndOfRunAction(const G4Run* aRun);

    void SetIO( g4rcIO *io ){ fIO = io; }

  private:
    G4Timer* timer;

    g4rcIO *fIO;
};

#endif

