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
    //Making world filled with air

    //World dimensions
    G4double worldSizeXY = 10.*cm;
    G4double worldSizeZ  = 10.*cm;


    //World material
    G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");

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

    //-------------------------------
    //Creating cell tube cylinder

    //Thickness of tube
    G4double thickness_cellTube = 0.1*mm;

    //Tube material
    G4Material* cellTube_mat = nist->FindOrBuildMaterial("G4_POLYSTYRENE");

    //Tube parameters
    cellTubeRMin = 0.5*mm; //Inner radius
    cellTubeHeight = 1.*mm; //Height in z-diraction
    G4double cellTubeRMax = cellTubeRMin + thickness_cellTube; //Outer radius
    G4double cellTubeSPhi = 0.*deg; //Start angle for Tube section
    G4double cellTubeEPhi = 360.*deg; //End angle for Tube section

    //Making Cell Tube solid
    G4Tubs* cellTubeSV =
    new G4Tubs("CellTube",
            cellTubeRMin,
            cellTubeRMax,
            0.5*cellTubeHeight,
            cellTubeSPhi,
            cellTubeEPhi);

    //Making Cell Tube logical
    G4LogicalVolume* cellTubeLV =
    new G4LogicalVolume(cellTubeSV,
                    cellTube_mat,
                    "CellTube");

    //Making Cell Tube physical
    new G4PVPlacement(0,
                G4ThreeVector(),
                cellTubeLV,
                "CellTube",
                worldLV,
                false,
                -2,
                fCheckOverlaps);

    //-------------------------------
    //Making water solution filling the cell tube

    //Solution material
    G4Material* solution_mat = nist->FindOrBuildMaterial("G4_WATER");

    //Solution parameters
    G4double solutionRMin = 0.*mm;
    G4double solutionRMax = cellTubeRMin;
    G4double solutionHeight = cellTubeHeight;
    G4double solutionSPhi = 0.*deg;
    G4double solutionEPhi = 360.*deg;

    //Making solution solid
    G4Tubs* solidSolution =
    new G4Tubs("Solution",
        solutionRMin,
        solutionRMax,
        0.5*solutionHeight,
        solutionSPhi,
        solutionEPhi);

    //Making solution logical
    G4LogicalVolume* logicSolution =
    new G4LogicalVolume(solidSolution,
                    solution_mat,
                    "Solution");

    //Making solution physical
    new G4PVPlacement(0,
                G4ThreeVector(),
                logicSolution,
                "Solution",
                worldLV,
                false,
                -3,
                fCheckOverlaps);


    //-------------------------------
    //Making Cell

    //cell parameters
    G4double cellRMin = 0.*um; //Inner radius
    cellRMax = 9.*um; //Outer radius
    G4double cellSPhi = 0.*deg; //Phi Start
    G4double cellEPhi = 360.*deg; //Phi end
    G4double cellSTheta = 0.*deg; //Start Theta
    G4double cellETheta = 180.*deg; //End Theta

    //Cell Material
    G4Material* cell_mat = nist->FindOrBuildMaterial("G4_WATER");

    //-------------------------------
    //Making Cell Membrane

    //Thickness of cell membrane
    thickness_membrane = 4.*nm;

    //Cell membrane parameters
    G4double membraneRMax = cellRMax;
    G4double membraneRMin = membraneRMax - thickness_membrane;
    G4double membraneSPhi = 0.*deg;
    G4double membraneEPhi = 360.*deg;
    G4double membraneSTheta = 0.*deg;
    G4double membraneETheta = 180.*deg;

    //Cell membrane material
    G4Material* membrane_mat = nist->FindOrBuildMaterial("G4_WATER");

    //Making membrane solid
    G4Sphere* membraneSV =
    new G4Sphere("Membrane",
            membraneRMin,
            membraneRMax,
            membraneSPhi,
            membraneEPhi,
            membraneSTheta,
            membraneETheta);

    //Making membrane logical
    G4LogicalVolume* membraneLV =
    new G4LogicalVolume(membraneSV,
                    membrane_mat,
                    "Membrane");

    auto cellMembraneVisAtt = new G4VisAttributes(G4Colour(0.0,0.5,1.0));
    // cellMembraneVisAtt->SetForceSolid(true);
    membraneLV->SetVisAttributes(cellMembraneVisAtt);

    //-------------------------------
    //Making Cell Nuclues

    //Cell Nucleus parameters
    G4double cellNucleusRMin = 0.*um; //Inner radius
    G4double cellNucleusRMax = 2.5*um; //Outer radius
    G4double cellNucleusSPhi = 0.*deg; //Phi Start
    G4double cellNucleusEPhi = 360.*deg; //Phi end
    G4double cellNucleusSTheta = 0.*deg; //Start Theta
    G4double cellNucleusETheta = 180.*deg; //End Theta

    //Cell Nucleus Material
    G4Material* cellNucleus_mat = nist->FindOrBuildMaterial("G4_WATER");

    //Makeing cell nucleus solid
    G4Sphere* cellNucleusSV =
    new G4Sphere("CellNucleus",
            cellNucleusRMin,
            cellNucleusRMax,
            cellNucleusSPhi,
            cellNucleusEPhi,
            cellNucleusSTheta,
            cellNucleusETheta);

    //Making cell nucleus logical
    G4LogicalVolume* cellNucleusLV =
    new G4LogicalVolume(cellNucleusSV,
                    cellNucleus_mat,
                    "CellNucleus");


    //-------------------------------
    //Placing 1000 cells uniformly within the cell tube

    //Number of cells to be placed
    G4double numberCells_sample = 1000000;
    G4double milliLiter_sample = 0.2;
    G4double cellDensity_sample = numberCells_sample/(milliLiter_sample*1000*mm*mm*mm);

    G4double volumeCellTube = CLHEP::pi*std::pow(cellTubeRMin,2.0)*cellTubeHeight;

    numberCells = volumeCellTube*cellDensity_sample;

    //Counter
    G4int cellCounter = 0;

    while(cellCounter < numberCells)
    {
        //Generating uniformly distributed position inside cell tube
        G4double r_cell = std::pow(CLHEP::RandFlat::shoot(), 1.0/2) * (cellTubeRMin - cellRMax);
        G4double theta_cell = CLHEP::RandFlat::shoot()*2*CLHEP::pi;
        G4double z_cell = (CLHEP::RandFlat::shoot() - 0.5)*(cellTubeHeight - 2*cellRMax - 0.1*mm);

        G4double x_cell = r_cell*std::cos(theta_cell);
        G4double y_cell = r_cell*std::sin(theta_cell);

        G4ThreeVector newVec = G4ThreeVector(x_cell,y_cell,z_cell);

        //Boolean, true if new position is found, or else false
        bool foundNewPosition = true;

        //Check for overlaps with other cells
        for(G4int i=0; i<cellPositions.size(); i++)
        {
            //Calculate distance between cell-centres
            G4ThreeVector storedVec = cellPositions[i];
            G4ThreeVector diffVec = storedVec - newVec;
            G4double diffDistance = diffVec.mag();

            //Make bolean false if overlap
            if(diffDistance <= 2*cellRMax)
            {
                foundNewPosition = false;
            }
        }

        //Place cell in new generated position
        if(foundNewPosition)
        {
            //-----------------------
            //Store new position
            cellPositions.push_back(newVec);

            //Update number of cells
            cellCounter += 1;

            //Copy number for cells
            G4int cellCopyNumber = cellCounter;

            //Copy number for membranes
            G4int membraneCopyNumber = cellCounter + numberCells;

            //Copy number for nucleus
            G4int nucleusCopynumber = cellCounter + 2*numberCells;

            //Copy number for cell nucleus
            G4int cellNucleusCopyNumber = numberCells + cellCounter;

            //-----------------------
            //Making cell solid
            G4Sphere* solidCell =
            new G4Sphere("Cell",
                        cellRMin,
                        cellRMax,
                        cellSPhi,
                        cellEPhi,
                        cellSTheta,
                        cellETheta);

            //Making cell logical shape
            G4LogicalVolume* logicCell =
            new G4LogicalVolume(solidCell,
                                cell_mat,
                                "Cell");
            //Placeing cell
            new G4PVPlacement(0,
                        newVec,
                        logicCell,
                        "Cell",
                        logicSolution,
                        false,
                        cellCopyNumber,
                        fCheckOverlaps);


            //----------------------------
            //Placing memebrane
            //Making membrane physical
            new G4PVPlacement(0,
                            G4ThreeVector(),
                            membraneLV,
                            "Membrane",
                            logicCell,
                            false,
                            membraneCopyNumber,
                            fCheckOverlaps);

            //------------------------
            //Place cell nuclei

            //Generate uniformly distributes cell nuclei
            G4double r_nucleus = std::pow(CLHEP::RandFlat::shoot(),1.0/3) * (cellRMax - thickness_membrane - cellNucleusRMax);
            G4double theta_nucleus = std::acos(1 - 2*CLHEP::RandFlat::shoot());
            G4double phi_nucleus = 2*CLHEP::pi*CLHEP::RandFlat::shoot();

            G4double x_nucleus = r_nucleus*std::sin(theta_nucleus)*std::cos(phi_nucleus);
            G4double y_nucleus = r_nucleus*std::sin(theta_nucleus)*std::sin(phi_nucleus);
            G4double z_nucleus = r_nucleus*std::cos(theta_nucleus);

            G4ThreeVector nucluesPosition = G4ThreeVector(x_nucleus,y_nucleus,z_nucleus);

            //Making cell nucleus physical
            new G4PVPlacement(0,
                            nucluesPosition,
                            cellNucleusLV,
                            "CellNucleus",
                            logicCell,
                            false,
                            nucleusCopynumber,
                            fCheckOverlaps);
        }
    }
    return worldPV;
}



