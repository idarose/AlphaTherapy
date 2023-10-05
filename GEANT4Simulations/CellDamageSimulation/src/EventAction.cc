#include "EventAction.hh"
#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4AnalysisManager.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(PrimaryGeneratorAction* primaryGeneratorAction) : fPrimaryGeneratorAction(primaryGeneratorAction)
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
    firstInteractionTimeVec.clear();
    firstInteractionVolumeVec.clear();

    // Resetting step number
    ResetStepNumber();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* event)
{
    auto analysisManager = G4AnalysisManager::Instance();


    // If simulating decays in solution, only store decays happening within first two hours
    if(fPrimaryGeneratorAction->GetInitialRadionuclide_location()==0)
    {
        if(firstInteractionTimeVec[0]/3600. < 2.0)
        {
            // Adding row
            analysisManager->AddNtupleRow();
        }
    }
    else if(fPrimaryGeneratorAction->GetInitialRadionuclide_location()==1)
    {
        // Adding row
        analysisManager->AddNtupleRow();
    }
    else if(fPrimaryGeneratorAction->GetInitialRadionuclide_location()==2)
    {
        // Adding row
        analysisManager->AddNtupleRow();
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::SetFirstInteractionInfo(G4double firstInteractionTime, G4int firstInteractionVolume)
{
    firstInteractionTimeVec.push_back(firstInteractionTime);
    firstInteractionVolumeVec.push_back(firstInteractionVolume);
}


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

