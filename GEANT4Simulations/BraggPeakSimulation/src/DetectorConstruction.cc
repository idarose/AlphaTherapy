#include "DetectorConstruction.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"
#include "Randomize.hh"
#include "G4ProductionCuts.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"





//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


DetectorConstruction::DetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    // Define volumes
    return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{
    auto nist = G4NistManager::Instance();


    //For setting a constant seed for the sample cell positions:

    G4int seed = 13732;
    CLHEP::HepRandom::setTheEngine(new CLHEP::MTwistEngine);
    // CLHEP::HepRandom::setTheSeed(time(NULL));
    CLHEP::HepRandom::setTheSeed(seed);


    //-------------------------------
    //Making world filled with water

    //World dimensions
    G4double worldSizeXY = 50.*cm;
    G4double worldSizeZ  = 50.*cm;


    //World material
    G4Material* world_mat = nist->FindOrBuildMaterial("G4_WATER");

    //Making world solid
    G4Box* worldS
    = new G4Box("World",           // its name
                 worldSizeXY/2, worldSizeXY/2, worldSizeZ/2); // its size

    //Making world logical
    G4LogicalVolume* worldLV
    = new G4LogicalVolume(
                 worldS,           // its solid
                 world_mat,  // its material
                 "World");         // its name

    //Making world physical
    G4PVPlacement* worldPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 worldLV,          // its logical volume
                 "World",          // its name
                 0,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps

    return worldPV;
}



