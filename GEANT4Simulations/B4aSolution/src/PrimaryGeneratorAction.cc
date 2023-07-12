#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4IonTable.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* detConstruction): fDetConstruction(detConstruction)
{
    G4int n_particle = 1;
    fParticleGun  = new G4ParticleGun(n_particle);

    // default particle kinematic
    //
    auto particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("alpha");
    fParticleGun->SetParticleDefinition(particleDefinition);
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
    fParticleGun->SetParticleEnergy(0.0*keV);
    // fParticleGun->SetParticleEnergy(1.0*eV);
    fParticleGun->SetParticleCharge(0);

    // // default particle kinematic
    // G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    // G4String particleName;


    // fParticleGun->SetParticleMomentumDirection(G4ThreeVector(1.0,0.0,0.0));
    // fParticleGun->SetParticleEnergy(0.0*MeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    //Generating new particle for new event

    //Making particle a 212Pb nucleus
    G4ParticleDefinition* ion_212Pb = G4IonTable::GetIonTable()->GetIon(83,212,0.0);

    fParticleGun->SetParticleDefinition(ion_212Pb);


    //-------------------------------
    // Distributing radionuclides randomly within the cell tube. There are methods in eventAction and steppingAction that sorts the decays happening in solution from decays happening inside cell



    //Cell tube parameters
    G4double cellTubeRMin_ = fDetConstruction->GetCellTubeRMin(); //Inner radius
    G4double cellTubeHeight_ = fDetConstruction->GetCellTubeHeight(); //Height in z-direction


    //------------------------
    //Generating random postion for particle inside cell Tube
    G4double r_particle = std::pow(G4UniformRand(), 1.0/2.0) * cellTubeRMin_;
    G4double theta_particle = G4UniformRand()*2*CLHEP::pi;
    G4double z_particle = (G4UniformRand() - 0.5)*(cellTubeHeight_);

    G4double x_particle = r_particle*std::cos(theta_particle);
    G4double y_particle = r_particle*std::sin(theta_particle);

    G4ThreeVector randPosition = G4ThreeVector(x_particle,y_particle,z_particle);


    //------------------------
    fParticleGun->SetParticlePosition(randPosition);
    fParticleGun->GeneratePrimaryVertex(anEvent);




    // //-----------------------
    // //Generating isotropic direction for momentum
    // G4double theta_momentum = std::acos(1 - 2*G4UniformRand());
    // G4double phi_momentum = 2*CLHEP::pi*G4UniformRand();

    // G4double x_momentum = std::sin(theta_momentum)*std::cos(phi_momentum);
    // G4double y_momentum = std::sin(theta_momentum)*std::sin(phi_momentum);
    // G4double z_momentum = std::cos(theta_momentum);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


