#include "../include/DecayDynamics.hpp"

//------------------–----------
DecayDynamics::DecayDynamics(int activitySample_in, std::string cellLine_in, std::string cellGeometry_in)
{
    //--------------------
    activitySample = activitySample_in;
    cellLine = cellLine_in;
    cellGeometry = cellGeometry_in;


    //-----------------------
    double VolumeSample = 0.2*1000; // mm^3
    double volumeCellTube = TMath::Pi()*std::pow(0.5,2.0)*1.0; // mm^3
    volumeRatio = volumeCellTube/VolumeSample;

    numberCells = 500000.*volumeRatio;

    //--------------------------
    double densityWater = 1000. ; // kg/m^3
    double radiusCell = 9.0e-6; // m
    double radiusCytoplasm = radiusCell - 4.0e-9; // m

    double radiusNucleus;

    if(cellGeometry=="D12RP"||cellGeometry=="D12CP")
    {
        radiusNucleus = 6.0e-6; // m
    }

    if(cellGeometry=="D5RP"||cellGeometry=="D5CP")
    {
        radiusNucleus = 2.5e-6; // m
    }

    massNucleus = (4./3.)*TMath::Pi()*std::pow(radiusNucleus,3.)*densityWater; // kg
    massCytoplasm = (4./3.)*TMath::Pi()*std::pow(radiusCytoplasm,3.)*densityWater - massNucleus; // kg
    massCell = (4./3.)*TMath::Pi()*std::pow(radiusCell,3.)*densityWater; // kg
    massMembrane = massCell - (4./3.)*TMath::Pi()*std::pow(radiusCytoplasm,3.)*densityWater; // kg
    // std::cout << massNucleus << std::endl;
}


//------------------–----------
std::vector<double> DecayDynamics::ReadFileFromMathematica(std::string filename)
{
    //----------------------
    //  Imports data from file where the first column contains data in string form
    //  and second column contains data in double form. Returns the columns of double
    //  values in the form of a vector

    std::vector<double> input_data;

    std::fstream myfile(filename, std::ios_base::in);


    if (myfile.is_open())
    {
        std::string line;
        std::string x;
        double y;

        while(std::getline(myfile,line))
        {
            std::stringstream mystream(line);
            mystream >> x >> y;
            input_data.push_back(y);
        }
    }
    else
    {
        std::cout << "Unable to open file " << filename << std::endl;
    }
    myfile.close();

    return input_data;

}


//------------------–----------
void DecayDynamics::LoadDataFromMathematicaCalculations(std::string filepathToMathematicaOutput)
{

    //----------------------
    // Loads data from calculations, assuming output files are structured as shown in the file "212PbDecayDynamics.nb"

    std::string filepathSolutionData = filepathToMathematicaOutput + "/" + cellLine + "/Solution/Activity_" + std::to_string(activitySample) + "kBq/NumberDecays.dat";

    // Importing data for decays occuring in solution
    std::vector<double> decayDataSolution = ReadFileFromMathematica(filepathSolutionData);
    numberDecays212PbInSolution1hTo2h = decayDataSolution[0];


    // If there are no radionuclides internalized there is no output file for "Cells"
    // So only read "Cells" files if there is uptake
    if(cellLine!="PC3_Flu")
    {
        // Importing data for decays occuring in cells
        std::string filepathCellData = filepathToMathematicaOutput + "/" + cellLine + "/Cells/Activity_" + std::to_string(activitySample) + "kBq/NumberDecays.dat";
        std::vector<double> decayDataCells = ReadFileFromMathematica(filepathCellData);

        numberDecays212PbInMembrane1h2To26h = decayDataCells[7];
        numberDecays212PbInCytoplasm1hTo26h = decayDataCells[14];
    }
    else
    {
        numberDecays212PbInMembrane1h2To26h = 0.0;
        numberDecays212PbInCytoplasm1hTo26h = 0.0;
    }
}
