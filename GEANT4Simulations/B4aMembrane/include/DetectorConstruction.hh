#ifndef B4DetectorConstruction_h
#define B4DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
# include "G4ThreeVector.hh"

class G4VPhysicalVolume;
class G4GlobalMagFieldMessenger;


class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction();
    ~DetectorConstruction() override;

  public:
    G4VPhysicalVolume* Construct() override;


    // Functions for extracting geometry constants from detectorConstruction
    G4int& GetNumberCells() {return numberCells;};
    std::vector<G4ThreeVector>& GetCellPositions() {return cellPositions;};

    G4double& GetThickness_cellMembrane() {return thickness_cellMembrane;};
    G4double& GetCellTubeRMin() {return cellTubeRMin;};
    G4double& GetCellTubeHeight() {return cellTubeHeight;};
    G4double& GetCellRMax() {return cellRMax;};

  private:
    G4VPhysicalVolume* DefineVolumes();

    G4bool fCheckOverlaps = false; // option to activate checking of volumes overlaps

    // Constants of geometry
    G4int numberCells;
    std::vector<G4ThreeVector> cellPositions;

    G4double cellTubeRMin;
    G4double cellTubeHeight;
    G4double cellRMax;
    G4double thickness_cellMembrane;
};




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

