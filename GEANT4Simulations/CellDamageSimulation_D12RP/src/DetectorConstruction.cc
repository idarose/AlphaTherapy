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

    auto worldVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
    worldLV->SetVisAttributes(worldVisAtt);

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

    auto cellCellTubeVisAtt = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
    cellTubeLV->SetVisAttributes(cellCellTubeVisAtt);

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
    G4Tubs* cellTubeSolutionSV =
    new G4Tubs("Solution",
        solutionRMin,
        solutionRMax,
        0.5*solutionHeight,
        solutionSPhi,
        solutionEPhi);

    //Making solution logical
    G4LogicalVolume* cellTubeSolutionLV =
    new G4LogicalVolume(cellTubeSolutionSV,
                    solution_mat,
                    "Solution");

    //Making solution physical
    new G4PVPlacement(0,
                G4ThreeVector(),
                cellTubeSolutionLV,
                "Solution",
                worldLV,
                false,
                -3,
                fCheckOverlaps);

    //-------------------------------
    //Making Cell Cytoplasm

    //Thickness of cell membrane
    thickness_cellMembrane = 4.*nm;

    //cell cytoplasm parameters
    G4double cellCytoplasmRMin = 0.*um; //Inner radius
    cellCytoplasmRMax = 9.*um - thickness_cellMembrane; //Outer radius
    G4double cellCytoplasmSPhi = 0.*deg; //Phi Start
    G4double cellCytoplasmEPhi = 360.*deg; //Phi end
    G4double cellCytoplasmSTheta = 0.*deg; //Start Theta
    G4double cellCytoplasmETheta = 180.*deg; //End Theta

    //Cell cytoplasm Material
    G4Material* cell_mat = nist->FindOrBuildMaterial("G4_WATER");

    // Vector to store cytoplasm logical volumes
    std::vector<G4LogicalVolume *> cellCytoplasmLVVec;


    //-------------------------------
    //Number of cells to be placed
    G4double numberCells_sample = 0.5*std::pow(10.0,6.0);
    G4double milliLiter_sample = 0.2;
    G4double volume_sample = milliLiter_sample*1000*mm*mm*mm;
    G4double cellDensity_sample = numberCells_sample/volume_sample;

    G4double volumeCellTube = CLHEP::pi*std::pow(cellTubeRMin,2.0)*cellTubeHeight;

    numberCells = volumeCellTube*cellDensity_sample;

    //Defining the G4Region for the cell cytoplasm
    G4Region* region_cellCytoplasm = new G4Region("region_cellCytoplasm");


    //Defining the G4ProductionCut for the cell cytoplasm
    G4ProductionCuts* cut_cellCytoplasm= new G4ProductionCuts;
    cut_cellCytoplasm->SetProductionCut(0.1*nm);

    //Linking cut_cellCytoplasm with region_cellCytoplasm
    region_cellCytoplasm->SetProductionCuts(cut_cellCytoplasm);

    //-------------------------------
    // Making and storing cell cytoplasm logical volumes
    for(int i=0; i<numberCells; i++)
    {
        //-----------------------
        //Making cell cytoplasm solid
        G4Sphere* cellCytoplasmSV =
        new G4Sphere("Cell",
                    cellCytoplasmRMin,
                    cellCytoplasmRMax,
                    cellCytoplasmSPhi,
                    cellCytoplasmEPhi,
                    cellCytoplasmSTheta,
                    cellCytoplasmETheta);

        //Making cell cytoplasm logical shape
        G4LogicalVolume* cellCytoplasmLV =
        new G4LogicalVolume(cellCytoplasmSV,
                            cell_mat,
                            "Cell");

        auto cellCytoplasmVisAtt = new G4VisAttributes(false);
        cellCytoplasmLV->SetVisAttributes(cellCytoplasmVisAtt);

        cellCytoplasmLVVec.push_back(cellCytoplasmLV);

        region_cellCytoplasm->AddRootLogicalVolume(cellCytoplasmLVVec.back());
    }


    //-------------------------------
    //Making Cell Membrane

    //Cell membrane parameters
    G4double cellMembraneRMax = cellCytoplasmRMax + thickness_cellMembrane;
    G4double cellMembraneRMin = cellCytoplasmRMax;
    G4double cellMembraneSPhi = 0.*deg;
    G4double cellMembraneEPhi = 360.*deg;
    G4double cellMembraneSTheta = 0.*deg;
    G4double cellMembraneETheta = 180.*deg;

    //Cell membrane material
    G4Material* cellMembrane_mat = nist->FindOrBuildMaterial("G4_WATER");

    //Making membrane solid
    G4Sphere* cellMembraneSV =
    new G4Sphere("cellMembrane",
            cellMembraneRMin,
            cellMembraneRMax,
            cellMembraneSPhi,
            cellMembraneEPhi,
            cellMembraneSTheta,
            cellMembraneETheta);

    //Making membrane logical
    G4LogicalVolume* cellMembraneLV =
    new G4LogicalVolume(cellMembraneSV,
                    cellMembrane_mat,
                    "cellMembrane");

    //Defining the G4Region for the cell membrane
    G4Region* region_cellMembrane = new G4Region("region_cellMembrane");
    region_cellMembrane->AddRootLogicalVolume(cellMembraneLV);

    //Defining the G4ProductionCut for the cell membrane
    G4ProductionCuts* cut_cellMembrane= new G4ProductionCuts;
    cut_cellMembrane->SetProductionCut(0.1*nm);

    //Linking cut_cellMembrane with region_cellMembrane
    region_cellMembrane->SetProductionCuts(cut_cellMembrane);



    auto cellMembraneVisAtt = new G4VisAttributes(G4Colour(0.0,0.5,1.0));
    cellMembraneLV->SetVisAttributes(cellMembraneVisAtt);

    //-------------------------------
    //Making Cell Nuclues

    //Cell Nucleus parameters
    G4double cellNucleusRMin = 0.*um; //Inner radius
    G4double cellNucleusRMax = 6.*um; //Outer radius
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

     //Defining the G4Region for the cell nucleus
    G4Region* region_cellNucleus = new G4Region("region_cellNucleus");
    region_cellNucleus->AddRootLogicalVolume(cellNucleusLV);

    //Defining the G4ProductionCut for the cell nucleus
    G4ProductionCuts* cut_cellNucleus= new G4ProductionCuts;
    cut_cellNucleus->SetProductionCut(0.1*nm);

    //Linking cut_cellNucleus with region_cellNucleus
    region_cellNucleus->SetProductionCuts(cut_cellNucleus);

    //Counter
    G4int numberCells_counter = 0;

    //-------------------------
    G4double cellRmax = cellCytoplasmRMax + thickness_cellMembrane;
    while(numberCells_counter < numberCells)
    {
        //Generating uniformly distributed position inside cell tube
        G4double r_cell = std::pow(CLHEP::RandFlat::shoot(), 1./2.) * (cellTubeRMin - cellRmax);
        G4double theta_cell = CLHEP::RandFlat::shoot()*2*CLHEP::pi;
        G4double z_cell = (CLHEP::RandFlat::shoot() - 0.5)*(cellTubeHeight - 2*cellRmax);

        G4double x_cell = r_cell*std::cos(theta_cell);
        G4double y_cell = r_cell*std::sin(theta_cell);

        G4ThreeVector newVec = G4ThreeVector(x_cell,y_cell,z_cell);

        //Boolean, true if new position is found, false if no new position is found
        bool foundNewPosition = true;

        //Check for overlaps with other cells
        for(G4int i=0; i<cellPositions.size(); i++)
        {
            //Calculate distance between cell-centres
            G4ThreeVector storedVec = cellPositions[i];
            G4ThreeVector diffVec = storedVec - newVec;
            G4double diffDistance = diffVec.mag();

            //Make bolean false if overlap
            if(diffDistance <= (2*cellRmax))
            {
                foundNewPosition = false;
            }
        }

        //Place cell in new generated position
        if(foundNewPosition)
        {
            G4LogicalVolume* cellCytoplasmLV = cellCytoplasmLVVec[numberCells_counter];

            //-----------------------
            //Store new position
            cellPositions.push_back(newVec);

            //Update number of cells
            numberCells_counter += 1;


            /*
            Cell ID numbers work like this:

            Cell IDs are numbered from 1 to number of cells.

            The cell membrane has a volume copy number equal to the cell ID number.
            The cell cytoplasm of a cell has a volume copy number equal to the cell ID number plus the number of cells.
            The cell nuleus has a volume copy number equal to the cell ID plus two times the number of cells.
            */

            //Copy number for cell membrane
            G4int cellMembraneCopyNumber = numberCells_counter;

            //Copy number for cell cytoplasm
            G4int cellCytoplasmCopyNumber = numberCells_counter + numberCells;

            //Copy number for nucleus
            G4int cellNucleusCopynumber = numberCells_counter + 2*numberCells;


            //----------------------------
            //Placing cell membrane
            new G4PVPlacement(0,
                            newVec,
                            cellMembraneLV,
                            "cellMembrane",
                            cellTubeSolutionLV,
                            false,
                            cellMembraneCopyNumber,
                            fCheckOverlaps);

            //----------------------------
            //Placeing cell cytoplasm
            new G4PVPlacement(0,
                        newVec,
                        cellCytoplasmLV,
                        "Cell",
                        cellTubeSolutionLV,
                        false,
                        cellCytoplasmCopyNumber,
                        fCheckOverlaps);


            //------------------------
            // Generate uniformly distributed cell nucleus
            // Will generate a random position within the cell cytoplasm

            G4double r_nucleus = std::pow(CLHEP::RandFlat::shoot(),1./3.) * (cellCytoplasmRMax - cellNucleusRMax);
            G4double theta_nucleus = std::acos(1. - 2.*CLHEP::RandFlat::shoot());
            G4double phi_nucleus = 2*CLHEP::pi*CLHEP::RandFlat::shoot();

            G4double x_nucleus = r_nucleus*std::sin(theta_nucleus)*std::cos(phi_nucleus);
            G4double y_nucleus = r_nucleus*std::sin(theta_nucleus)*std::sin(phi_nucleus);
            G4double z_nucleus = r_nucleus*std::cos(theta_nucleus);

            G4ThreeVector nucleusPosition = G4ThreeVector(x_nucleus,y_nucleus,z_nucleus);

            //Placing cell nucleus
            new G4PVPlacement(0,
                            nucleusPosition,
                            cellNucleusLV,
                            "CellNucleus",
                            cellCytoplasmLV,
                            false,
                            cellNucleusCopynumber,
                            fCheckOverlaps);
        }
    }
    return worldPV;
}



