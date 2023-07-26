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

// https://stackoverflow.com/questions/6942273/how-to-get-a-random-element-from-a-c-container

#include  <random>
#include  <iterator>


// Methods to extract random element from a vector
template<typename Iter, typename RandomGenerator>
Iter select_randomly(Iter start, Iter end, RandomGenerator& gen) {
    std::uniform_int_distribution<> dis(0, std::distance(start, end) - 1);
    std::advance(start, dis(gen));
    return start;
}

template<typename Iter>
Iter select_randomly(Iter start, Iter end) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    return select_randomly(start, end, gen);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* detConstruction): fDetConstruction(detConstruction)
{
    //-------------------------------
    G4int n_particle = 1;
    fParticleGun  = new G4ParticleGun(n_particle);

    //-------------------------------
    // default particle kinematics
    auto particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("alpha");
    fParticleGun->SetParticleDefinition(particleDefinition);
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
    fParticleGun->SetParticleEnergy(0.0*keV);
    // fParticleGun->SetParticleEnergy(1.0*eV);
    fParticleGun->SetParticleCharge(0);


    //-------------------------------
    numberCells = fDetConstruction->GetNumberCells();
    cellPositions = fDetConstruction->GetCellPositions();

    cellTubeRMin = fDetConstruction->GetCellTubeRMin();
    cellTubeHeight = fDetConstruction->GetCellTubeHeight();
    cellCytoplasmRMax = fDetConstruction->GetCellCytoplasmRMax();
    thickness_cellMembrane = fDetConstruction->GetThickness_cellMembrane();
    cellRMax = thickness_cellMembrane + cellCytoplasmRMax;
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

    //------------------------
    //Generating random postion for particle inside cell Tube
    G4double r_particle = std::pow(G4UniformRand(), 1.0/2.0) * cellTubeRMin;
    G4double theta_particle = G4UniformRand()*2*CLHEP::pi;
    G4double z_particle = (G4UniformRand() - 0.5)*(cellTubeHeight);

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


