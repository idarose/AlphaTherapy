#ifndef B4aSteppingAction_h
#define B4aSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"


class DetectorConstruction;
class EventAction;
class PrimaryGeneratorAction;

class SteppingAction : public G4UserSteppingAction
{
public:
  SteppingAction(DetectorConstruction* detConstruction,
                 EventAction* eventAction, PrimaryGeneratorAction* primaryGeneratorAction);
  ~SteppingAction() override;

  void UserSteppingAction(const G4Step* step) override;

private:
  DetectorConstruction* fDetConstruction;
  EventAction* fEventAction = nullptr;

  // Volume type numbers defined
  G4int volumeTypeMembrane = 1;
  G4int volumeTypeCytoplasm = 2;
  G4int volumeTypeNucleus = 3;


  PrimaryGeneratorAction* fPrimaryGeneratorAction;

  G4int numberCells;
};



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
