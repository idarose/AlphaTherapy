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
#include  <TFile.h>
#include  <TROOT.h>


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

void PrimaryGeneratorAction::LoadDecayCurveRadionuclide(G4int initialRadionuclide_location, G4int sampleActivity, std::string cellLineName)
{
    // Method to load file containing number of decays in solution and inside cells
    // as a function of time. This curve is used to draw random decay times in
    // GeneratingPrimaries().

    std::string filePath;
    // TF1* fDecayCurveRadionuclide = nullptr;

    if(initialRadionuclide_location == 0)
    {
        filePath = "../DecaysPb212Bi212/Decays212Pb212Bi_" + cellLineName + "_Solution_Activity_" + std::to_string(sampleActivity) + "kBq.root";

        TFile* inputFile = new TFile(filePath.c_str(), "READ");
        inputFile->GetObject("fDecaysSolution", fDecayCurve);
        inputFile->Close();

        minTime = 1.;
        maxTime = 2.;
    }
    if(initialRadionuclide_location == 1 || initialRadionuclide_location == 2)
    {
        filePath = "../DecaysPb212Bi212/Decays212Pb212Bi_"  + cellLineName + "_Cells_Activity_" + std::to_string(sampleActivity) + "kBq.root";

        TFile* inputFile = new TFile(filePath.c_str(), "READ");
        inputFile->GetObject("fDecaysCells", fDecayCurve);
        inputFile->Close();

        minTime = 1.;
        maxTime = 26.;
    }

    if(fDecayCurve)
    {
        G4cout << " DecayCurvePointer was initialized!" << G4endl;
        // fDecayCurve = std::shared_ptr<TF1>(fDecayCurveRadionuclide);
        maxValueDecayCurve = 1.01*fDecayCurve->GetMaximum();
    }
    else
    {
        G4cout << " DecayCurvePointer was NOT initialized!" << G4endl;
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* detConstruction):
fDetConstruction(detConstruction),
initialRadionuclide_excitationEnergy(0)
{
    TF1::DefaultAddToGlobalList(false);
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
    // gROOT->GetListOfFunctions()->Remove(fDecayCurve.get());
    // gROOT->GetListOfFunctions()->Remove(fDecayCurve);

    delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    //Generating new particle for new event

    bool foundDecayTime = false;

    auto GenerateDecayTime = [&]()
    {
        double r = G4UniformRand()*maxValueDecayCurve;
        double x = minTime + (maxTime - minTime)*G4UniformRand();
        double y = fDecayCurve->Eval(x);
        if(r<=y)
        {
            fParticleGun->SetParticleTime(x*3600.*s);
            foundDecayTime = true;
        }
    };


    //-------------------------------
    // Distributing radionuclides randomly within the cell tube (aka. solution)
    if(initialRadionuclide_location==0)
    {
        //------------------------
        // Generating decay time

        while(!foundDecayTime)
        {
            GenerateDecayTime();
        }


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
    // Distributing radionuclides randomly within the cell membrane
    if(initialRadionuclide_location==1)
    {
        //------------------------
        // Generating decay time

        while(!foundDecayTime)
        {
            GenerateDecayTime();
        }

        //-----------------------
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
    // Generating uniformly distibuted particles inside the cell cytoplasm
    if(initialRadionuclide_location==2)
    {
        //------------------------
        // Generating decay time

        while(!foundDecayTime)
        {
            GenerateDecayTime();
        }

        //-----------------------
        //Generating random radius
        G4double radius = std::pow(G4UniformRand(), 1.0/3.0)*(cellCytoplasmRMax);
        G4double theta = std::acos(1 - 2.0*G4UniformRand());
        G4double phi = 2*G4UniformRand()*CLHEP::pi;

        G4double x_particle = radius*std::sin(theta)*std::cos(phi);
        G4double y_particle = radius*std::sin(theta)*std::sin(phi);
        G4double z_particle = radius*std::cos(theta);


        // Selecting radom cell
        G4ThreeVector centerOfCell = *select_randomly(cellPositions.begin(), cellPositions.end());


        // Generating position inside membrane of this cell
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
    ion->SetPDGLifeTime(0.*ns);
}

//----------------------------------------------------------------------------

void PrimaryGeneratorAction::SetSampleActivity(G4int value)
{
    sampleActivity = value;
}


//----------------------------------------------------------------------------

void PrimaryGeneratorAction::SetCellLineName(std::string value)
{
    cellLineName = value;
}

//----------------------------------------------------------------------------

void PrimaryGeneratorAction::SetDecayCurveRadionuclide()
{
    LoadDecayCurveRadionuclide(initialRadionuclide_location, sampleActivity, cellLineName);
}
