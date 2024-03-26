#include "include/CellSurvival.hpp"
#include "include/SurvivalFit.hpp"
#include "include/HitAnalysis.hpp"
#include "include/DoseAnalysis.hpp"
#include <vector>
#include <tuple>

int main(int argc, char *argv[])
{


    // std::vector<std::tuple<double,double,double>> data_cellSurvival_C4_2;
    // data_cellSurvival_C4_2.push_back(std::make_tuple(5.0,    0.86,   std::sqrt(std::pow(0.117,2.) + std::pow(0.1*0.86,2.))));
    // data_cellSurvival_C4_2.push_back(std::make_tuple(10.0,    0.69,    std::sqrt(std::pow(0.014,2.0) + std::pow(0.1*0.69,2.))));
    // data_cellSurvival_C4_2.push_back(std::make_tuple(25.0,    0.51,    std::sqrt(std::pow(0.096,2.) + std::pow(0.1*0.51,2.))));
    // data_cellSurvival_C4_2.push_back(std::make_tuple(50.0,    0.26,   std::sqrt(std::pow(0.035,2.) + std::pow(0.1*0.26,2.))));
    // data_cellSurvival_C4_2.push_back(std::make_tuple(75.0,    0.13,    std::sqrt(std::pow(0.031,2.) + std::pow(0.1*0.13,2.))));
    // data_cellSurvival_C4_2.push_back(std::make_tuple(100.0,    0.08,    std::sqrt(std::pow(0.035,2.) + std::pow(0.1*0.08,2.))));
    // data_cellSurvival_C4_2.push_back(std::make_tuple(150.0,    0.04,    std::sqrt(std::pow(0.024,2.) + std::pow(0.1*0.04,2.))));

    std::vector<std::tuple<double,double,double>> data_cellSurvival_C4_2;
    data_cellSurvival_C4_2.push_back(std::make_tuple(5.0,    0.86,   std::sqrt(std::pow(0.117,2.) + std::pow(0.05*0.86,2.))));
    data_cellSurvival_C4_2.push_back(std::make_tuple(10.0,    0.69,    std::sqrt(std::pow(0.014,2.0) + std::pow(0.05*0.69,2.))));
    data_cellSurvival_C4_2.push_back(std::make_tuple(25.0,    0.51,    std::sqrt(std::pow(0.096,2.) + std::pow(0.05*0.51,2.))));
    data_cellSurvival_C4_2.push_back(std::make_tuple(50.0,    0.26,   std::sqrt(std::pow(0.035,2.) + std::pow(0.05*0.26,2.))));
    data_cellSurvival_C4_2.push_back(std::make_tuple(75.0,    0.13,    std::sqrt(std::pow(0.031,2.) + std::pow(0.05*0.13,2.))));
    data_cellSurvival_C4_2.push_back(std::make_tuple(100.0,    0.08,    std::sqrt(std::pow(0.035,2.) + std::pow(0.05*0.08,2.))));
    data_cellSurvival_C4_2.push_back(std::make_tuple(150.0,    0.04,    std::sqrt(std::pow(0.024,2.) + std::pow(0.05*0.04,2.))));


    std::vector<std::tuple<double,double,double>> data_cellSurvival_PC3_Flu;
    data_cellSurvival_PC3_Flu.push_back(std::make_tuple(10.0,    0.955,    std::sqrt(std::pow(0.0955,2.) + std::pow(0.1*0.955,2.))));
    data_cellSurvival_PC3_Flu.push_back(std::make_tuple(25.0,    0.724,    std::sqrt(std::pow(0.0724,2.) + std::pow(0.1*0.724,2.))));
    data_cellSurvival_PC3_Flu.push_back(std::make_tuple(50.0,    0.733,    std::sqrt(std::pow(0.0733,2.) + std::pow(0.1*0.733,2.))));
    data_cellSurvival_PC3_Flu.push_back(std::make_tuple(75.0,    0.798,    std::sqrt(std::pow(0.0798,2.) + std::pow(0.1*0.798,2.))));
    data_cellSurvival_PC3_Flu.push_back(std::make_tuple(100.0,    0.729,    std::sqrt(std::pow(0.0729,2.) + std::pow(0.1*0.729,2.))));
    data_cellSurvival_PC3_Flu.push_back(std::make_tuple(150.0,    0.690,    std::sqrt(std::pow(0.0690,2.) + std::pow(0.1*0.690,2.))));

    std::vector<std::tuple<double,double,double>> data_cellSurvival_PC3_PIP;
    data_cellSurvival_PC3_PIP.push_back(std::make_tuple(10.0,    0.630,    std::sqrt(std::pow(0.063,2.) + std::pow(0.1*0.630,2.))));
    data_cellSurvival_PC3_PIP.push_back(std::make_tuple(25.0,    0.317,    std::sqrt(std::pow(0.05,2.) + std::pow(0.1*0.317,2.))));
    data_cellSurvival_PC3_PIP.push_back(std::make_tuple(50.0,    0.071,    std::sqrt(std::pow(0.05,2.) + std::pow(0.1*0.071,2.))));
    data_cellSurvival_PC3_PIP.push_back(std::make_tuple(75.0,    0.032,    std::sqrt(std::pow(0.05,2.) + std::pow(0.1*0.032,2.))));
    data_cellSurvival_PC3_PIP.push_back(std::make_tuple(100.0,    0.014,    std::sqrt(std::pow(0.05,2.) + std::pow(0.1*0.014,2.))));
    data_cellSurvival_PC3_PIP.push_back(std::make_tuple(150.0,    0.006,    std::sqrt(std::pow(0.05,2.) + std::pow(0.1*0.006,2.))));

    //------------------–----------
    if (argc != 5) {
        std::cerr << "Usage: " << argv[0] << " string_argument integer_argument" << std::endl;
        return 1;
    }

    //------------------–----------
    std::string cellLine = argv[1];
    std::string cellGeometry = argv[2];
    std::string modelType = argv[3];
    int volumeType;

    //--------------------------
    try {
        volumeType = std::stoi(argv[4]);
        if(volumeType<=0||volumeType>4){
            throw std::invalid_argument("Invalid colume number, has to be within {1,2,3}.");
        }
    } catch (const std::invalid_argument& e) {
        std::cerr << "Error: third argument is not a valid integer! " << e.what() << std::endl;
        return 4;
    } catch (const std::out_of_range& e) {
        std::cerr << "Error: third argument is out of range for an integer! " << e.what() << std::endl;
        return 5;
    }


    //--------------------------
    std::vector<std::tuple<double,double,double>> data_CellSurvival;

    if(cellLine=="C4_2")
    {
        data_CellSurvival = data_cellSurvival_C4_2;
    }
    else if(cellLine=="PC3_PIP")
    {
        data_CellSurvival = data_cellSurvival_PC3_PIP;
    }
    else if(cellLine=="PC3_Flu")
    {
        data_CellSurvival = data_cellSurvival_PC3_Flu;
    }

    //-----------------------------
    CellSurvival cellSurvival = CellSurvival(cellLine, cellGeometry);
    cellSurvival.AddCellSurvivalData(data_CellSurvival);

    SurvivalFit survivalFit;
    survivalFit.FitCellSurvival(cellSurvival, modelType, volumeType);

    HitAnalysis hitAnalysis(survivalFit);
    hitAnalysis.MakeHitAnalysis(41);

    DoseAnalysis doseAnalysis(survivalFit);
    doseAnalysis.MakeDoseAnalysis();


    //-------------------
    std::string volumeName;
    if(volumeType==1)
    {
        volumeName = "Membrane";
    }
    if(volumeType==2)
    {
        volumeName = "Cytoplasm";
    }
    if(volumeType==3)
    {
        volumeName = "Nucleus";
    }
    if(volumeType==4)
    {
        volumeName = "TotalCell";
    }

    std::cout << volumeName << std::endl;
    std::string outputDir = "Output_" + cellGeometry + "/";
    std::string outputName = outputDir + "Output_" + cellLine + "_" + volumeName + "_" + modelType + ".root";
    TFile *outputFile = new TFile(outputName.c_str(), "RECREATE");

    survivalFit.WriteToFile();
    hitAnalysis.WriteToFile();
    doseAnalysis.WriteToFile();

    outputFile->Write();
    outputFile->Close();

    return 0;
}