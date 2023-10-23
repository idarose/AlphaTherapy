#include "TFile.h"
#include "TCanvas.h"
#include "TStopwatch.h"
#include "TGraphAsymmErrors.h"
#include "TCutG.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH2Poly.h"
#include "TGraphErrors.h"
#include "TGraph2D.h"
#include "TArrow.h"
#include "TLatex.h"
#include "TF1Convolution.h"
#include "TFitResult.h"
#include "TParameter.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TMathText.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TAxis.h"
#include "TGaxis.h"
#include <iomanip>
#include <sstream>



//----------------------------
class CellSurvival
{
    public:
        CellSurvival(std::vector<std::tuple<double,double,double>>  dataCellSurvival_in);
        void GraphCellSurvivalData();
        TGraphErrors* GetGraphObject(){return graphClonogenicSurvival;};
        std::vector<std::tuple<double,double,double>> GetCellSurvivalData(){return dataCellSurvival;};


    private:
        std::vector<std::tuple<double,double,double>> dataCellSurvival;
        TGraphErrors* graphClonogenicSurvival;
};



CellSurvival::CellSurvival(std::vector<std::tuple<double,double,double>>  dataCellSurvival_in)
{
    dataCellSurvival = dataCellSurvival_in;
}



void CellSurvival::GraphCellSurvivalData()
{
    int graphPointN_data = 0;

    for(auto & entry : dataCellSurvival)
    {
        double activity = std::get<0>(entry);
        double cellSurvival = std::get<1>(entry);
        double cellSurvivalUncertainty = std::get<2>(entry);

        graphClonogenicSurvival->SetPoint(graphPointN_data, activity, cellSurvival);
        graphClonogenicSurvival->SetPointError(graphPointN_data, 0.0, cellSurvivalUncertainty);

        graphPointN_data++;
    }
}


//----------------------------
class EnergyDepositionHistograms
{
    public:
        EnergyDepositionHistograms(std::string cellLine_in);
        void LoadHistogramsFromAnalysis(CellSurvival cellSurvivalInstance);

    private:
        std::vector<std::tuple<double,TH1D*>> activity_histograms_vec;
        std::string cellLine;
};



EnergyDepositionHistograms::EnergyDepositionHistograms(std::string cellLine_in)
{
    cellLine = cellLine_in;
}



void EnergyDepositionHistograms::LoadHistogramsFromAnalysis(CellSurvival cellSurvivalInstance)
{
    std::string generalFilePath = "../OutputAnalysisCode/Output_" + cellLine + "_";
    std::vector<std::tuple<double,double,double>> cellSurvivalData = cellSurvivalInstance.GetCellSurvivalData();


    std::string filePath;
    for(int i=0; i<cellSurvivalData.size(); i++)
    {
        // Extracting activity
        double activity = std::get<0>(cellSurvivalData[i]);

        // Making filename
        filePath = generalFilePath + std::to_string(((int)activity)) + "kBq.root";
        TFile* inputFile = new TFile(filePath.c_str(), "READ");

        // Extracting energy deposition histogram for nucleus
        TH1D* hEnergyDeps_Nucleus = nullptr;
        std::string histogramName = "hEnergyDeps_212Pb_" + cellLine + "_" + std::to_string(((int)activity)) +"kBq_Nucleus";
        inputFile->GetObject(histogramName.c_str(), hEnergyDeps_Nucleus);
        hEnergyDeps_Nucleus->SetDirectory(0);
        inputFile->Close();

        // Storing histogram
        activity_histograms_vec.push_back(std::make_tuple(activity, hEnergyDeps_Nucleus));
    }

    for(auto & entry : activity_histograms_vec)
    {
        double activity = std::get<0>(entry);
        auto hist = std::get<1>(entry);

        double integral = hist->Integral();

        std::cout << "activity: " << activity << ",\t integral: " << integral << std::endl;
    }


}




//----------------------------
void mainFittingCode()
{
    std::vector<std::tuple<double,double,double>> data_cellSurvival_PC3_Flu;
    data_cellSurvival_PC3_Flu.push_back(make_tuple(10.0,    0.955,    0.0955));
    data_cellSurvival_PC3_Flu.push_back(make_tuple(25.0,    0.724,    0.0724));
    data_cellSurvival_PC3_Flu.push_back(make_tuple(50.0,    0.733,    0.0733));
    data_cellSurvival_PC3_Flu.push_back(make_tuple(75.0,    0.798,    0.0798));
    data_cellSurvival_PC3_Flu.push_back(make_tuple(100.0,    0.729,    0.0720));
    data_cellSurvival_PC3_Flu.push_back(make_tuple(150.0,    0.690,    0.0698));

    CellSurvival cellSurvival_PC3_Flu = CellSurvival(data_cellSurvival_PC3_Flu);

    EnergyDepositionHistograms energyDepHistograms_PC3_Flu =    EnergyDepositionHistograms("PC3_Flu");
    energyDepHistograms_PC3_Flu.LoadHistogramsFromAnalysis(cellSurvival_PC3_Flu);

};