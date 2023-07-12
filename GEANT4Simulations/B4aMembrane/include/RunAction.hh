#ifndef B4RunAction_h
#define B4RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include <vector>

class EventAction;
class G4Run;


class RunAction : public G4UserRunAction
{
  public:
    RunAction(EventAction* eventAction);
    ~RunAction() override;

    EventAction *fEventAction;

    void BeginOfRunAction(const G4Run*) override;
    void EndOfRunAction(const G4Run*) override;
};



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

