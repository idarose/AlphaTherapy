#include "PrimaryGeneratorAction.hh"
#include "PrimaryGeneratorMessenger.hh"
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

PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* detConstruction):
fDetConstruction(detConstruction),
initialRadionuclide_excitationEnergy(0)
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

    fPrimaryGeneratorMessenger = new PrimaryGeneratorMessenger(this);
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
    // G4ParticleDefinition* ion_212Pb = G4IonTable::GetIonTable()->GetIon(82,212,0.0);
    // fParticleGun->SetParticleDefinition(ion_212Pb);


    //-------------------------------
    // Distributing radionuclides randomly within the cell tube. There are methods in eventAction and steppingAction that sorts the decays happening in solution from decays happening inside cell
    if(initialRadionuclide_location==0)
    {
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
    }

    //-------------------------------
    // Distributing radionuclides randomly within the cell memebrane
    if(initialRadionuclide_location==1)
    {
        //Generating random radius
        G4double T = (std::pow(cellRMax, 3.0) - std::pow(cellCytoplasmRMax, 3.0))*G4UniformRand() + std::pow(cellCytoplasmRMax, 3.0);
        G4double radius = std::pow(T, 1.0/3.0);

        // Generating random angles
        G4double theta = std::acos(1 - (2.0*G4UniformRand()));
        G4double phi = 2*CLHEP::pi*G4UniformRand();

        // Transforming from spherical coord, to cartesian
        G4double x_particle = radius*std::sin(theta)*std::cos(phi);
        G4double y_particle = radius*std::sin(theta)*sin(phi);
        G4double z_particle = radius*std::cos(theta);

        // Selecting radom cell
        G4ThreeVector centerOfCell = *select_randomly(cellPositions.begin(), cellPositions.end());


        // Generating position inside memebrane of this cell
        G4ThreeVector particlePosRelativeToCenterCell = G4ThreeVector(x_particle,y_particle,z_particle);

        // Finding position relative to world volume
        G4ThreeVector particlePosRelativeWorld = particlePosRelativeToCenterCell + centerOfCell;


        fParticleGun->SetParticlePosition(particlePosRelativeWorld);
        fParticleGun->GeneratePrimaryVertex(anEvent);
    }

    //-------------------------------
    // Generating unifmroly distibuted particles inside the cell cytoplasm. Conditions in eventAction and steppingAction discard events where first interaction is in cell nuclei
    if(initialRadionuclide_location==2)
    {
        G4double radius = std::pow(G4UniformRand(), 1.0/3.0)*(cellCytoplasmRMax);
        G4double theta = std::acos(1 - 2.0*G4UniformRand());
        G4double phi = 2*G4UniformRand()*CLHEP::pi;

        G4double x_particle = radius*std::sin(theta)*std::cos(phi);
        G4double y_particle = radius*std::sin(theta)*std::sin(phi);
        G4double z_particle = radius*std::cos(theta);


        // Selecting radom cell
        G4ThreeVector centerOfCell = *select_randomly(cellPositions.begin(), cellPositions.end());


        // Generating position inside memebrane of this cell
        G4ThreeVector particlePosRelativeToCenterCell = G4ThreeVector(x_particle,y_particle,z_particle);

        // Finding the absolute position
        G4ThreeVector particlePosRelativeWorld = particlePosRelativeToCenterCell + centerOfCell;


        fParticleGun->SetParticlePosition(particlePosRelativeWorld);
        fParticleGun->GeneratePrimaryVertex(anEvent);
    }
}

//----------------------------------------------------------------------------

void PrimaryGeneratorAction::SetInitialRadionuclide_Z(G4int value)
{
    initialRadionuclide_Z = value;
}

//----------------------------------------------------------------------------

void PrimaryGeneratorAction::SetInitialRadionuclide_A(G4int value)
{
    initialRadionuclide_A = value;
}

//----------------------------------------------------------------------------

void PrimaryGeneratorAction::SetInitialRadionuclide_excitationEnergy(G4double value)
{
    initialRadionuclide_excitationEnergy = value;
}

//----------------------------------------------------------------------------

void PrimaryGeneratorAction::SetInitialRadionuclide_location(G4int value)
{
    initialRadionuclide_location = value;
}

//----------------------------------------------------------------------------

void PrimaryGeneratorAction::DefineInitialRadionuclide()
{
    ion = G4IonTable::GetIonTable()->GetIon(initialRadionuclide_Z, initialRadionuclide_A, initialRadionuclide_excitationEnergy);
    fParticleGun->SetParticleDefinition(ion);
}

//----------------------------------------------------------------------------

