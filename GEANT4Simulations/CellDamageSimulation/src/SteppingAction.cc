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


    //Getting the interaction time
    G4double interactionTime = step->GetPreStepPoint()->GetGlobalTime()/s;



    //------------------------
    // Storing interaction time and volume type of first interaction
    if(fEventAction->GetStepNumber() == 0)
    {
        //Checking if the interaction took place in the membrane of the cell
        if(volumeCopyNumber <= numberCells && volumeCopyNumber>= 1)
        {
            fEventAction->SetFirstInteractionInfo(step->GetPreStepPoint()->GetGlobalTime()/s, 1);
        }

        //Checking if the interaction took place in the cytoplasm of the cell
        else if(volumeCopyNumber >= (numberCells+1) && volumeCopyNumber <= 2*numberCells)
        {
            fEventAction->SetFirstInteractionInfo(step->GetPreStepPoint()->GetGlobalTime()/s, 2);
        }

        //Checking if the interaction took place in the nucleus of the cell
        else if(volumeCopyNumber >= (2*numberCells+1) && volumeCopyNumber <= 3*numberCells)
        {
            fEventAction->SetFirstInteractionInfo(step->GetPreStepPoint()->GetGlobalTime()/s, 3);
        }

        // Checking if interaction took place in solution
        else if(volumeCopyNumber == -3)
        {
            fEventAction->SetFirstInteractionInfo(step->GetPreStepPoint()->GetGlobalTime()/s, 0);
        }
    }

    // Add step number
    fEventAction->AddStepNumber();




    //----------------------------------
    /*
    Cell ID numbers work like this:

    Cell IDs are numbered from 1 to number of cells.

    The cell membrane has a volume copy number equal to the cell ID number.
    The cell cytoplasm of a cell has a volume copy number equal to the cell ID number plus the number of cells.
    The cell nuleus has a volume copy number equal to the cell ID plus two times the number of cells.
    */



    //----------------------------------
    // Lambda function for filling vectors in EventAction with interaction information
    auto FillVectors = [&]()
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

            // //Getting the interaction time
            // G4double interactionTime = step->GetPreStepPoint()->GetGlobalTime()/s;

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

            // //Getting the interaction time
            // G4double interactionTime = step->GetPreStepPoint()->GetGlobalTime()/s;

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

            // //Getting the interaction time
            // G4double interactionTime = step->GetPreStepPoint()->GetGlobalTime()/s;

            // Calculating cell ID number
            G4int cellID = volumeCopyNumber - 2*numberCells;

            // Storing information
            fEventAction->StoreInteractionInformation(energyDeposition, cellID, volumeTypeNucleus, kineticEnergy, particleType, interactionTime);
        }
    };



    //----------------------------------
    // Filling vectors

    // If initial radionuclide position is in solution, only store interactions taking place
    // within first hour
    if(fPrimaryGeneratorAction->GetInitialRadionuclide_location()==0)
    {
        if(interactionTime/3600. <= 1.0)
        {
            FillVectors();
        }
    }
    // Else, store every interaction
    else
    {
        FillVectors();
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


