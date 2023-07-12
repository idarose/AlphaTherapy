#ifndef B4PrimaryGeneratorAction_h
#define B4PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"


class G4ParticleGun;
class G4Event;
class DetectorConstruction;


class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
    public:
      PrimaryGeneratorAction(DetectorConstruction* detConstruction);
      ~PrimaryGeneratorAction() override;

      void GeneratePrimaries(G4Event* event) override;

      // set methods
      void SetRandomFlag(G4bool value);

    private:
      G4ParticleGun* fParticleGun = nullptr; // G4 particle gun
      DetectorConstruction* fDetConstruction;
};


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
