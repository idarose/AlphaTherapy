#include "EventAction.hh"
#include "RunAction.hh"

#include "G4AnalysisManager.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* /*event*/)
{
    // Clearing vectors for new event
    energyDepVec.clear();
    positionZVec.clear();
    stepLengthVec.clear();

    // ResetStep();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* event)
{
    auto analysisManager = G4AnalysisManager::Instance();

    // Adding row
    analysisManager->AddNtupleRow();

}


void EventAction::StoreInteractionInformation(G4double energyDep, G4double positionZ, G4double stepLength)
{
    //Adding values to vectors
    energyDepVec.push_back(energyDep);
    positionZVec.push_back(positionZ);
    stepLengthVec.push_back(stepLength);
}

