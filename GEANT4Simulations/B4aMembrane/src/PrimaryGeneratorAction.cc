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
    //Generating new pericle for new event

    //Making particle a 212Pb nucleus
    G4ParticleDefinition* ion_212Pb = G4IonTable::GetIonTable()->GetIon(83,212,0.0);

    fParticleGun->SetParticleDefinition(ion_212Pb);



    //-------------------------------
    // Distributing radionuclides randomly within the cell memebrane


    //cell parameters
    G4double cellRMax_ = fDetConstruction->GetCellRMax(); //Outer radius
    G4double thickness_membrane_ = fDetConstruction->GetThickness_cellMembrane();
    G4double cytoplasmRMax_ = cellRMax_ - thickness_membrane_;


    // Vector with cell positions
    std::vector<G4ThreeVector> cellPositions_ = fDetConstruction->GetCellPositions();


    // Transformed distribution to generate random radii
    G4double T = (std::pow(cellRMax_, 3.0) - std::pow(cytoplasmRMax_, 3.0))*G4UniformRand() + std::pow(cytoplasmRMax_, 3.0);

    // Generating radius relative to cell center
    G4double radius = std::pow(T, 1.0/3.0);


    // Generating angles
    G4double theta = std::acos(1 - (2.0*G4UniformRand()));
    G4double phi = 2*CLHEP::pi*G4UniformRand();




    // Transforming from spherical coord, to cartesian
    G4double x_particle = radius*std::sin(theta)*std::cos(phi);
    G4double y_particle = radius*std::sin(theta)*sin(phi);
    G4double z_particle = radius*std::cos(theta);



    // Selecting radom cell
    G4ThreeVector centerOfCell = *select_randomly(cellPositions_.begin(), cellPositions_.end());



    // Generating position inside memebrane of this cell
    G4ThreeVector particlePosRelativeToCenterCell = G4ThreeVector(x_particle,y_particle,z_particle);

    // Finding the absolute position
    G4ThreeVector particlePosRelativeWorld = particlePosRelativeToCenterCell + centerOfCell;



    fParticleGun->SetParticlePosition(particlePosRelativeWorld);
    fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


