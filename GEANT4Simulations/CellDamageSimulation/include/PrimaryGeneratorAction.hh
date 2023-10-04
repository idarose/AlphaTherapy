#ifndef B4PrimaryGeneratorAction_h
#define B4PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include <vector>


class G4ParticleGun;
class G4Event;
class G4ParticleDefinition;

class DetectorConstruction;
class PrimaryGeneratorMessenger;


class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
    public:
        PrimaryGeneratorAction(DetectorConstruction* detConstruction);
        ~PrimaryGeneratorAction() override;

        void GeneratePrimaries(G4Event* event) override;

        // set methods
        void SetRandomFlag(G4bool value);

        void SetInitialRadionuclide_Z(G4int value);
        void SetInitialRadionuclide_A(G4int value);
        void SetInitialRadionuclide_excitationEnergy(G4double value);
        void SetInitialRadionuclide_location(G4int value);
        void DefineInitialRadionuclide();

        int GetInitialRadionuclide_location(){return initialRadionuclide_location;};

        G4int       initialRadionuclide_Z;
        G4int       initialRadionuclide_A;
        G4double    initialRadionuclide_excitationEnergy;
        G4int       initialRadionuclide_location;

    private:
        G4ParticleGun* fParticleGun = nullptr; // G4 particle gun
        DetectorConstruction* fDetConstruction;
        PrimaryGeneratorMessenger* fPrimaryGeneratorMessenger;

        G4int numberCells;
        std::vector<G4ThreeVector> cellPositions;

        G4double cellTubeRMin;
        G4double cellTubeHeight;
        G4double cellCytoplasmRMax;
        G4double cellRMax;
        G4double thickness_cellMembrane;

        G4ParticleDefinition* ion;
};


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

