#ifndef B4RunAction_h
#define B4RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include <vector>

#include "RunMessenger.hh"

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

    RunMessenger* fRunMessenger;

    G4String fOutputFileName;
    G4int fSeedN;

    void SetOutputFileName(G4String name);
    void SetSeed(G4int seedN);

};



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

