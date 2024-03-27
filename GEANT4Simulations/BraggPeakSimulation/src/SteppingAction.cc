#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"
#include "G4SystemOfUnits.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4Step.hh"
#include "G4RunManager.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


SteppingAction::SteppingAction(DetectorConstruction* detConstruction,
                               EventAction* eventAction, PrimaryGeneratorAction* primaryGeneratorAction)
  : fDetConstruction(detConstruction),
    fEventAction(eventAction), fPrimaryGeneratorAction(primaryGeneratorAction)
{
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{
    //Getting the volume copy number for the interaction

    G4int particleType = step->GetTrack()->GetParticleDefinition()->GetPDGEncoding();
    if(particleType==1000020040)
    {
        // G4StepPoint* preStepPoint = step->GetPreStepPoint();
        G4StepPoint* postStepPoint = step->GetPostStepPoint();

        // G4double preZPosition = preStepPoint->GetPosition().z();
        G4double postZPosition = postStepPoint->GetPosition().z();

        G4double preStepEnergy = step->GetPreStepPoint()->GetTotalEnergy();
        G4double postStepEnergy = step->GetPostStepPoint()->GetTotalEnergy();
        G4double energyLoss = preStepEnergy - postStepEnergy;

        G4double stepLength = step->GetStepLength();

        G4double energyDep = step->GetTotalEnergyDeposit();

        fEventAction->StoreInteractionInformation(energyLoss/MeV,postZPosition/um, stepLength/um);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


