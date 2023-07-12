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
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{
    //Getting number of cells
    G4int numberCells = fDetConstruction->GetNumberCells();

    //Getting the volume copy number for the interaction
    G4int volumeCopyNumber = step->GetPreStepPoint()->GetTouchableHandle()
      ->GetVolume()->GetCopyNo();

    //Getting the energy deposited for the interaction, in MeV
    G4double energyDeposition = step->GetTotalEnergyDeposit()/MeV;

    //Getting the particle type number for the interaction
    G4int particleType = step->GetTrack()->GetParticleDefinition()->GetPDGEncoding();

    //Getting the kinetic energy of the particle for the interaction, in MeV
    G4double kineticEnergy = step->GetPreStepPoint()->GetKineticEnergy()/MeV;

    //Getting the interaction time
    G4double interactionTime = step->GetPreStepPoint()->GetGlobalTime()/s;

    G4int cellID;


    /*
    Cell ID numbers work like this:

    Cell IDs are numbered from 1 to number of cells.

    The cell cytoplasm of a cell has a volume copy number equal to the cell ID number.
    The cell membrane has a volume copy number equal to the cell ID plus the number of cells.
    The cell nuleus has a volume copy number equal to the cell ID plus two times the number of cells.
    */



    //------------------------
    // Checking if first interaction took place inside a cell nucleus
    if(fEventAction->GetVolumeTypeVec().size() == 0)
    {
        if(volumeCopyNumber >= (2*numberCells+1) && volumeCopyNumber <= 3*numberCells)
        {
            //First interaction was in cell nucleus
            fEventAction->ChangeBooleanFirstInteractionNotInCellNucleus(true);
        }
    }



    if(fEventAction->GetBooleanFirstInteractionNotInCellNucleus())
    {
        //----------------------------------
        //Checking if the interaction took place in the membrane of the cell
        if(volumeCopyNumber >= (numberCells+1) && volumeCopyNumber <= 2*numberCells)
        {
            cellID = volumeCopyNumber - numberCells;

            fEventAction->StoreInteractionInformation(energyDeposition, cellID, volumeTypeMembrane, kineticEnergy, particleType, interactionTime);
        }

        //----------------------------------
        //Checking if the interaction took place in the cytoplasm
        else if(volumeCopyNumber <= numberCells && volumeCopyNumber >= 1)
        {
            cellID = volumeCopyNumber;

            fEventAction->StoreInteractionInformation(energyDeposition, cellID, volumeTypeCytoplasm, kineticEnergy, particleType, interactionTime);
        }

        //----------------------------------
        //Checking if the interaction took place in the nucleus
        else if(volumeCopyNumber >= (2*numberCells+1) && volumeCopyNumber <= 3*numberCells)
        {
            cellID = volumeCopyNumber - 2*numberCells;

            fEventAction->StoreInteractionInformation(energyDeposition, cellID, volumeTypeNucleus, kineticEnergy, particleType, interactionTime);
        }
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


