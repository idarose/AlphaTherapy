#include "../include/CellSurvival.hpp"
#include "../include/SurvivalFit.hpp"
#include "../include/HitAnalysis.hpp"
#include <vector>
#include <tuple>


int main()
{


    std::vector<std::tuple<double,double,double>> data_cellSurvival_C4_2;

    data_cellSurvival_C4_2.push_back(std::make_tuple(5.0,    0.86,   std::sqrt(std::pow(0.117,2.) + std::pow(0.1*0.86,2.))));
    data_cellSurvival_C4_2.push_back(std::make_tuple(10.0,    0.69,    std::sqrt(std::pow(0.014,2.0) + std::pow(0.1*0.69,2.))));
    data_cellSurvival_C4_2.push_back(std::make_tuple(25.0,    0.51,    std::sqrt(std::pow(0.096,2.) + std::pow(0.1*0.51,2.))));
    data_cellSurvival_C4_2.push_back(std::make_tuple(50.0,    0.26,   std::sqrt(std::pow(0.035,2.) + std::pow(0.1*0.26,2.))));
    data_cellSurvival_C4_2.push_back(std::make_tuple(75.0,    0.13,    std::sqrt(std::pow(0.031,2.) + std::pow(0.1*0.13,2.))));
    data_cellSurvival_C4_2.push_back(std::make_tuple(100.0,    0.08,    std::sqrt(std::pow(0.035,2.) + std::pow(0.1*0.08,2.))));
    data_cellSurvival_C4_2.push_back(std::make_tuple(150.0,    0.04,    std::sqrt(std::pow(0.024,2.) + std::pow(0.1*0.04,2.))));


    CellSurvival cellSurvival_C4_2 = CellSurvival("C4_2", "D12CP");
    cellSurvival_C4_2.AddCellSurvivalData(data_cellSurvival_C4_2);


    SurvivalFit Fit_C4_2;
    Fit_C4_2.FitCellSurvival(cellSurvival_C4_2, "LM", 3);

    HitAnalysis HitAnalysis_C4_2 = HitAnalysis(Fit_C4_2);
    HitAnalysis_C4_2.MakeHitAnalysis(41);

    std::string outputName = "Output_AnalyseCellSurvival_Test.root";
    TFile *outputFile = new TFile(outputName.c_str(), "RECREATE");

    Fit_C4_2.WriteToFile();
    HitAnalysis_C4_2.WriteToFile();

    outputFile->Write();
    outputFile->Close();


    // std::vector<std::tuple<double,double,double>> data_cellSurvival_PC3_PIP;

    // // data_cellSurvival_PC3_PIP.push_back(make_tuple(0.0,    1.0,    0.05));
    // data_cellSurvival_PC3_PIP.push_back(std::make_tuple(10.0,    0.630,    std::sqrt(std::pow(0.063,2.) + std::pow(0.1*0.630,2.))));
    // data_cellSurvival_PC3_PIP.push_back(std::make_tuple(25.0,    0.317,    std::sqrt(std::pow(0.05,2.) + std::pow(0.1*0.317,2.))));
    // data_cellSurvival_PC3_PIP.push_back(std::make_tuple(50.0,    0.071,    std::sqrt(std::pow(0.05,2.) + std::pow(0.1*0.071,2.))));
    // data_cellSurvival_PC3_PIP.push_back(std::make_tuple(75.0,    0.032,    std::sqrt(std::pow(0.05,2.) + std::pow(0.1*0.032,2.))));
    // data_cellSurvival_PC3_PIP.push_back(std::make_tuple(100.0,    0.014,    std::sqrt(std::pow(0.05,2.) + std::pow(0.1*0.014,2.))));
    // data_cellSurvival_PC3_PIP.push_back(std::make_tuple(150.0,    0.006,    std::sqrt(std::pow(0.05,2.) + std::pow(0.1*0.006,2.))));


    // CellSurvival cellSurvival_PC3_PIP = CellSurvival("PC3_PIP");
    // cellSurvival_PC3_PIP.AddCellSurvivalData(data_cellSurvival_PC3_PIP);


    // SurvivalFit Fit_PC3_PIP;
    // Fit_PC3_PIP.FitCellSurvival(cellSurvival_PC3_PIP, "LM", 3);

    // HitAnalysis HitAnalysis_PC3_PIP = HitAnalysis(Fit_PC3_PIP);
    // HitAnalysis_PC3_PIP.MakeHitAnalysis(41);

    // std::string outputName = "Output_AnalyseCellSurvival_Test.root";
    // TFile *outputFile = new TFile(outputName.c_str(), "RECREATE");

    // Fit_PC3_PIP.WriteToFile();
    // HitAnalysis_PC3_PIP.WriteToFile();

    // outputFile->Write();
    // outputFile->Close();

    return 0;
}