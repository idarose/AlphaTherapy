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
    cellIDVec.clear();
    volumeTypeVec.clear();
    kineticEnergyVec.clear();
    particleTypeVec.clear();
    interactionTimeVec.clear();

    // Default value of boolean
    booleanFirstInteractionInSolution = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* event)
{
    auto analysisManager = G4AnalysisManager::Instance();


    //Adding a row in TTree file only if first interaction was in solution
    if(booleanFirstInteractionInSolution)
    {
        if(energyDepVec.size()>0)
        {
            analysisManager->AddNtupleRow();
        }
    }


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::StoreInteractionInformation(G4double energyDep, G4int cellID, G4int volumeType, G4double kineticEnergy, G4int particleType, G4double interactionTime)
{
    //Adding values to vectors
    energyDepVec.push_back(energyDep);
    cellIDVec.push_back(cellID);
    volumeTypeVec.push_back(volumeType);
    kineticEnergyVec.push_back(kineticEnergy);
    particleTypeVec.push_back(particleType);
    interactionTimeVec.push_back(interactionTime);
}

