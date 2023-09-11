#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"
#include "G4SystemOfUnits.hh"

#include "G4Step.hh"
#include "G4RunManager.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction* detConstruction,
                               EventAction* eventAction)
  : fDetConstruction(detConstruction),
    fEventAction(eventAction)
{
    numberCells = fDetConstruction->GetNumberCells();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{
    //Getting the volume copy number for the interaction
    G4int volumeCopyNumber = step->GetPreStepPoint()->GetTouchableHandle()
      ->GetVolume()->GetCopyNo();


    /*
    Cell ID numbers work like this:

    Cell IDs are numbered from 1 to number of cells.


    The cell membrane has a volume copy number equal to the cell ID number.
    The cell cytoplasm of a cell has a volume copy number equal to the cell ID number plus the number of cells.
    The cell nuleus has a volume copy number equal to the cell ID plus two times the number of cells.
    */


    // Checking is the first interaction was in the solution
    if(fEventAction->GetStepNumber() == 0)
    {
        if(volumeCopyNumber == -3)
        {
            fEventAction->BooleanFirstInteractionInSolutionMakeTrue();
            // G4cout << "First interaction in solution, bool: " << fEventAction->GetBooleanFirstInteractionInSolution() << G4endl;
            // G4cout << "First interaction in solution : " << volumeCopyNumber << " bool: " << fEventAction->GetBooleanFirstInteractionInSolution() << G4endl;
        }
        else
        {
            fEventAction->BooleanFirstInteractionInSolutionMakeFalse();
            // G4cout << "First interaction not in solution " << volumeCopyNumber <<  " bool : " << fEventAction->GetBooleanFirstInteractionInSolution() << G4endl;
        }
    }

    // G4cout << numberCells << G4endl;
    // Update step number
    fEventAction->AddStepNumber();

    // Storing interaction information only if first interaction was in solution
    if(fEventAction->GetBooleanFirstInteractionInSolution())
    {
        //----------------------------------
        //Checking if the interaction took place in the membrane of the cell
        if(volumeCopyNumber <= numberCells && volumeCopyNumber>= 1)
        {
            //Getting the energy deposited for the interaction, in MeV
            G4double energyDeposition = step->GetTotalEnergyDeposit()/MeV;

            //Getting the particle type number for the interaction
            G4int particleType = step->GetTrack()->GetParticleDefinition()->GetPDGEncoding();

            //Getting the kinetic energy of the particle for the interaction, in MeV
            G4double kineticEnergy = step->GetPreStepPoint()->GetKineticEnergy()/MeV;

            //Getting the interaction time, in seconds
            G4double interactionTime = step->GetPreStepPoint()->GetGlobalTime()/s;

            // Calculating cell ID number
            G4int cellID = volumeCopyNumber;


            // Storing information
            fEventAction->StoreInteractionInformation(energyDeposition, cellID, volumeTypeMembrane, kineticEnergy, particleType, interactionTime);
        }

        //----------------------------------
        //Checking if the interaction took place in the cytoplasm of the cell
        if(volumeCopyNumber >= (numberCells+1) && volumeCopyNumber <= 2*numberCells)
        {
            //Getting the energy deposited for the interaction, in MeV
            G4double energyDeposition = step->GetTotalEnergyDeposit()/MeV;

            //Getting the particle type number for the interaction
            G4int particleType = step->GetTrack()->GetParticleDefinition()->GetPDGEncoding();

            //Getting the kinetic energy of the particle for the interaction, in MeV
            G4double kineticEnergy = step->GetPreStepPoint()->GetKineticEnergy()/MeV;

            //Getting the interaction time
            G4double interactionTime = step->GetPreStepPoint()->GetGlobalTime()/s;

            // Calculating cell ID number
            G4int cellID = volumeCopyNumber - numberCells;

            // Storing information
            fEventAction->StoreInteractionInformation(energyDeposition, cellID, volumeTypeCytoplasm, kineticEnergy, particleType, interactionTime);
        }

        //----------------------------------
        //Checking if the interaction took place in the nucleus
        else if(volumeCopyNumber >= (2*numberCells+1) && volumeCopyNumber <= 3*numberCells)
        {
            //Getting the energy deposited for the interaction, in MeV
            G4double energyDeposition = step->GetTotalEnergyDeposit()/MeV;

            //Getting the particle type number for the interaction
            G4int particleType = step->GetTrack()->GetParticleDefinition()->GetPDGEncoding();

            //Getting the kinetic energy of the particle for the interaction, in MeV
            G4double kineticEnergy = step->GetPreStepPoint()->GetKineticEnergy()/MeV;

            //Getting the interaction time
            G4double interactionTime = step->GetPreStepPoint()->GetGlobalTime()/s;

            // Calculating cell ID number
            G4int cellID = volumeCopyNumber - 2*numberCells;

            // Storing information
            fEventAction->StoreInteractionInformation(energyDeposition, cellID, volumeTypeNucleus, kineticEnergy, particleType, interactionTime);
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


