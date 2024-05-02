// Make this appear first!
#include "G4Timer.hh"

#include "g4rcRunAction.hh"
#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "g4rcIO.hh"

g4rcRunAction::g4rcRunAction()
{
  timer = new G4Timer;
}

g4rcRunAction::~g4rcRunAction()
{
  delete timer;
}

void g4rcRunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  //  timer->Start();
  fIO->InitializeTree();

}

void g4rcRunAction::EndOfRunAction(const G4Run* aRun)
{
  timer->Stop();
  G4cout << "number of event = " << aRun->GetNumberOfEvent() << G4endl;
  //       << " " << *timer << G4endl;

  fIO->WriteTree();
}

