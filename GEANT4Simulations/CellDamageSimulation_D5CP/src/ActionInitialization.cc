#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"
#include "DetectorConstruction.hh"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ActionInitialization::ActionInitialization
                            (DetectorConstruction* detConstruction)
 : fDetConstruction(detConstruction)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ActionInitialization::~ActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ActionInitialization::BuildForMaster() const
{
    auto primaryGeneratorAction = new PrimaryGeneratorAction(fDetConstruction);
    auto eventAction = new EventAction();
    SetUserAction(new RunAction(eventAction));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ActionInitialization::Build() const
{
    auto primaryGeneratorAction = new PrimaryGeneratorAction(fDetConstruction);
    SetUserAction(primaryGeneratorAction);
    auto eventAction = new EventAction();
    SetUserAction(new RunAction(eventAction));
    SetUserAction(eventAction);
    SetUserAction(new SteppingAction(fDetConstruction,eventAction,primaryGeneratorAction));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


