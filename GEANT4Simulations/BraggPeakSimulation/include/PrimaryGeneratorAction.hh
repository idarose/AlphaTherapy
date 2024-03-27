#ifndef B4PrimaryGeneratorAction_h
#define B4PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include <vector>
#include <TF1.h>
#include "G4Types.hh"

#include "G4AutoLock.hh"

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

        void SetInitialAlphaEnergy(G4double value);

        G4double    initialAlphaEnergy;

    private:
        G4ParticleGun* fParticleGun = nullptr; // G4 particle gun
        DetectorConstruction* fDetConstruction;
        PrimaryGeneratorMessenger* fPrimaryGeneratorMessenger;

        G4ParticleDefinition* ion;
};


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

