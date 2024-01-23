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
        CellSurvival(std::string cellLine_in);

        void AddCellSurvivalData(std::vector<std::tuple<double,double,double>>  dataCellSurvival_in);

        std::vector<std::tuple<double,double,double>> GetCellSurvivalData(){return dataCellSurvival;};
        std::string GetCellLine(){return cellLine;};


    private:
        std::vector<std::tuple<double,double,double>> dataCellSurvival;
        std::string cellLine;
};

//----------------------------
CellSurvival::CellSurvival(std::string cellLine_in)
{
    cellLine = cellLine_in;
}

//----------------------------
void CellSurvival::AddCellSurvivalData(std::vector<std::tuple<double,double,double>>  dataCellSurvival_in)
{
    dataCellSurvival = dataCellSurvival_in;
}



//----------------------------
class EnergyDepositionHistograms
{
    public:
        EnergyDepositionHistograms(double t);
        void LoadHistogramsFromAnalysis(CellSurvival cellSurvivalInstance);
        std::vector<std::tuple<double,TH1D*,TH1D*>> GetEnergyDepositionHistograms(){return activity_histograms_vec;};

    private:
        double n;
        std::vector<std::tuple<double,TH1D*,TH1D*>> activity_histograms_vec;
};


EnergyDepositionHistograms::EnergyDepositionHistograms(double t)
{
    n = t;
}


//----------------------------
void EnergyDepositionHistograms::LoadHistogramsFromAnalysis(CellSurvival cellSurvivalInstance)
{
    // General filepath to output
    std::string generalFilePath = "../AnalysisCode/Output_" + cellSurvivalInstance.GetCellLine() + "_";

    std::vector<std::tuple<double,double,double>> cellSurvivalData = cellSurvivalInstance.GetCellSurvivalData();

    // Specific filepath
    std::string filePath;

    // Looping through all activities
    for(int i=0; i<cellSurvivalData.size(); i++)
    {
        // Extracting activity
        double activity = std::get<0>(cellSurvivalData[i]);

        if(activity<=0.)
        {
            std::cout << "Called on zero" << std::endl;

            // Extracting activity
            double activity = std::get<0>(cellSurvivalData[i+1]);

            // Making filename
            filePath = generalFilePath + std::to_string(((int)activity)) + "kBq.root";

            // Read file containing histogram
            TFile* inputFile = new TFile(filePath.c_str(), "READ");

            // Extracting energy deposition histogram for nucleus with keV binning
            TH1D* hEnergyDeps_Nucleus_keV = nullptr;
            std::string histogramName_keV = "hEnergyDeps_212Pb_" + cellSurvivalInstance.GetCellLine() + "_" + std::to_string(((int)activity)) +"kBq_Nucleus_keVBinning";

            inputFile->GetObject(histogramName_keV.c_str(), hEnergyDeps_Nucleus_keV);
            hEnergyDeps_Nucleus_keV->SetDirectory(0);
            hEnergyDeps_Nucleus_keV->Reset();

            // Extracting energy deposition histogram for nucleus with eVbinning
            TH1D* hEnergyDeps_Nucleus_eV = nullptr;
            std::string histogramName_eV = "hEnergyDeps_212Pb_" + cellSurvivalInstance.GetCellLine() + "_" + std::to_string(((int)activity)) +"kBq_Nucleus_eVBinning";

            inputFile->GetObject(histogramName_eV.c_str(), hEnergyDeps_Nucleus_eV);
            hEnergyDeps_Nucleus_eV->SetDirectory(0);
            inputFile->Close();

            hEnergyDeps_Nucleus_eV->Reset();

            activity_histograms_vec.push_back(std::make_tuple(0.,hEnergyDeps_Nucleus_eV,hEnergyDeps_Nucleus_keV));
        }
        else
        {
            // Making filename
            filePath = generalFilePath + std::to_string(((int)activity)) + "kBq.root";

            // Read file containing histogram
            TFile* inputFile = new TFile(filePath.c_str(), "READ");

            // Extracting energy deposition histogram for nucleus with keV binning
            TH1D* hEnergyDeps_Nucleus_keV = nullptr;
            std::string histogramName_keV = "hEnergyDeps_212Pb_" + cellSurvivalInstance.GetCellLine() + "_" + std::to_string(((int)activity)) +"kBq_Nucleus_keVBinning";

            inputFile->GetObject(histogramName_keV.c_str(), hEnergyDeps_Nucleus_keV);
            hEnergyDeps_Nucleus_keV->SetDirectory(0);

            // Extracting energy deposition histogram for nucleus with eV binning
            TH1D* hEnergyDeps_Nucleus_eV = nullptr;
            std::string histogramName_eV = "hEnergyDeps_212Pb_" + cellSurvivalInstance.GetCellLine() + "_" + std::to_string(((int)activity)) +"kBq_Nucleus_eVBinning";

            inputFile->GetObject(histogramName_eV.c_str(), hEnergyDeps_Nucleus_eV);
            hEnergyDeps_Nucleus_eV->SetDirectory(0);
            inputFile->Close();

            // Storing histogram
            activity_histograms_vec.push_back(std::make_tuple(activity, hEnergyDeps_Nucleus_eV, hEnergyDeps_Nucleus_keV));
        }
    }

    for(auto & entry : activity_histograms_vec)
    {
        double activity = std::get<0>(entry);
        auto hist = std::get<2>(entry);

        double integral = hist->Integral();

        std::cout << "activity: " << activity << ",\t integral: " << integral << std::endl;
    }
}


void FitCellSurvival(CellSurvival cellSurvivalInstance, std::string modelName)
{

    //--------------------------
    // Graphing cell survival data
    TGraphErrors *grClonogenicSurvival = new TGraphErrors();
    grClonogenicSurvival->SetName("grClonogenicSurvival");

    std::vector<std::tuple<double, double, double>> dataCellSurvival = cellSurvivalInstance.GetCellSurvivalData();

    auto GraphData = [&](TGraphErrors* grClonogenicSurvival, std::vector<std::tuple<double, double, double>> dataCellSurvival)
    {
        int graphPointN_data = 0;

        for(auto & entry : dataCellSurvival)
        {
            double activity = std::get<0>(entry);
            double cellSurvival = std::get<1>(entry);
            double cellSurvivalUncertainty = std::get<2>(entry);

            grClonogenicSurvival->SetPoint(graphPointN_data, activity, cellSurvival);
            grClonogenicSurvival->SetPointError(graphPointN_data, 0.0, cellSurvivalUncertainty);

            graphPointN_data++;
        }
    };

    GraphData(grClonogenicSurvival, dataCellSurvival);


    //------------------------------
    // Loading energy deposition histograms
    EnergyDepositionHistograms energyDepHistograms = EnergyDepositionHistograms(1.0);
    energyDepHistograms.LoadHistogramsFromAnalysis(cellSurvivalInstance);
    std::vector<std::tuple<double, TH1D*, TH1D*>> energyDepHistogramsVec = energyDepHistograms.GetEnergyDepositionHistograms();

    // ------------------------------
    // Lambda function to calculate the cell survival fraction for an energydepostion histogram
    auto CalculateCellSurvivalFraction = [&](TH1D* hEnergyDeposition_eV, TH1D* hEnergyDeposition_keV, double* par)
    {
        double alpha = par[0];
        double beta = par[1];

        // double fractionOfComponentsHit_eVHistogram = hEnergyDeposition_eV->Integral();
        double fractionOfComponentsHit_keVHistogram = hEnergyDeposition_keV->Integral();

        // std::cout << "eV : " << fractionOfComponentsHit_eVHistogram << " keV : " << fractionOfComponentsHit_keVHistogram << std::endl;
        double fractionOfComponentsHit = fractionOfComponentsHit_keVHistogram;

        double fractionOfComponentsMissed = (1.0 - fractionOfComponentsHit);

        // All cells that are missed are immidiately added as survived
        double fractionOfTotalCellsSurviving = fractionOfComponentsMissed;

        //------------------------------
        // Looping over eV binned histogram
        for(int i=0; i<hEnergyDeposition_eV->GetNbinsX(); i++)
        {
            double energyDeposition = hEnergyDeposition_eV->GetBinCenter(i+1);
            double fractionOfTotalCells = hEnergyDeposition_eV->GetBinContent(i+1);

            double cellSurvivalFraction = TMath::Exp(-(alpha*energyDeposition + beta*TMath::Power(energyDeposition, 2.0)));

            double fractionOfTotalCells_survivingFraction = cellSurvivalFraction*fractionOfTotalCells;

            fractionOfTotalCellsSurviving += fractionOfTotalCells_survivingFraction;
        }

        // Looping over keV binned histogram
        for(int i=1000; i<hEnergyDeposition_keV->GetNbinsX(); i++)
        {
            double energyDeposition = hEnergyDeposition_keV->GetBinCenter(i+1);
            double fractionOfTotalCells = hEnergyDeposition_keV->GetBinContent(i+1);

            double cellSurvivalFraction = TMath::Exp(-(alpha*energyDeposition + beta*TMath::Power(energyDeposition, 2.0)));

            double fractionOfTotalCells_survivingFraction = cellSurvivalFraction*fractionOfTotalCells;

            fractionOfTotalCellsSurviving += fractionOfTotalCells_survivingFraction;
        }

        // ------------------------------
        return fractionOfTotalCellsSurviving;

    };


    auto GenerateGraph_CellSurvivalFraction = [&](std::vector<std::tuple<double,TH1D*,TH1D*>> vec, double *par)
    {
        TGraphErrors gr;
        int graphPointN = 0;

        for(auto & entry : vec)
        {
            double activity = std::get<0>(entry);
            auto histogram_eV = std::get<1>(entry);
            auto histogram_keV = std::get<2>(entry);

            double cellSurvival = CalculateCellSurvivalFraction(histogram_eV, histogram_keV, par);

            gr.SetPoint(graphPointN, activity, cellSurvival);

            graphPointN++;
        }

        return gr;
    };

    int nParameters = 2;
    double savedParameters[nParameters];
    TGraphErrors grFittedCellSurvival;

    auto fFittedCellSurvival = new TF1("fFittedCellSurvival",
        [&](double *x, double *p)
        {
            double activity = x[0];

            bool foundChangeInParameters = false;

            for(int i=0; i<nParameters; i++)
            {
                if(savedParameters[i]!=p[i])
                {
                    foundChangeInParameters = true;
                    savedParameters[i] = p[i];
                }
            }

            if(foundChangeInParameters)
            {
                grFittedCellSurvival = GenerateGraph_CellSurvivalFraction(energyDepHistogramsVec, p);
            }

            double cellSurvival = grFittedCellSurvival.Eval(activity, 0, "");

            return cellSurvival;
        }, 0.0, 150.0, nParameters);

    fFittedCellSurvival->SetParLimits(0,0.0,1.0e+5);
    fFittedCellSurvival->SetParLimits(1,0.0,1.0e+5);
    fFittedCellSurvival->SetNpx(10000);
    fFittedCellSurvival->SetParameter(0, 1.0e+00);


    if(modelName=="LQM")
    {
        fFittedCellSurvival->SetParameter(1, 1.0e+00);
    }
    else if(modelName=="LM")
    {
        std::cout << "YEPP" << std::endl;
        fFittedCellSurvival->FixParameter(1, 0.0e+00);
    }

    std::cout << "--------------------------------" << std::endl;
    std::cout << "Fitting  cell line " << cellSurvivalInstance.GetCellLine() << " using " << modelName << " : " << std::endl;

    grClonogenicSurvival->Fit(fFittedCellSurvival, "", "", 0.0, 150.0);

    double chi_sq = fFittedCellSurvival->GetChisquare();
    double deg_freedom = ((double)fFittedCellSurvival->GetNDF());

    std::cout << "Chi squared: " << chi_sq << std::endl;
    std::cout << "Deg Freedom: " << deg_freedom << std::endl;
    std::cout << "Reduced Chi squared: " << chi_sq/deg_freedom << std::endl;

    std::string outputFileName = "Output_Fitting_" + cellSurvivalInstance.GetCellLine() + "_" + modelName + ".root";
    TFile *outputFile = new TFile(outputFileName.c_str(), "RECREATE");
    grClonogenicSurvival->Write();
    fFittedCellSurvival->Write();
    grFittedCellSurvival.Write();

    outputFile->Write();
    outputFile->Close();
}



//----------------------------
void mainFittingCode()
{
    // std::vector<std::tuple<double,double,double>> data_cellSurvival_PC3_Flu;
    // data_cellSurvival_PC3_Flu.push_back(make_tuple(10.0,    0.955,    0.0955));
    // data_cellSurvival_PC3_Flu.push_back(make_tuple(25.0,    0.724,    0.0724));
    // data_cellSurvival_PC3_Flu.push_back(make_tuple(50.0,    0.733,    0.0733));
    // data_cellSurvival_PC3_Flu.push_back(make_tuple(75.0,    0.798,    0.0798));
    // data_cellSurvival_PC3_Flu.push_back(make_tuple(100.0,    0.729,    0.0720));
    // data_cellSurvival_PC3_Flu.push_back(make_tuple(150.0,    0.690,    0.0698));

    // CellSurvival cellSurvival_PC3_Flu = CellSurvival(data_cellSurvival_PC3_Flu);

    std::vector<std::tuple<double,double,double>> data_cellSurvival_C4_2;
    data_cellSurvival_C4_2.push_back(make_tuple(0.0,    1.0,    0.0));
    data_cellSurvival_C4_2.push_back(make_tuple(5.0,    0.86,   0.117));
    data_cellSurvival_C4_2.push_back(make_tuple(10.0,    0.69,    0.014));
    data_cellSurvival_C4_2.push_back(make_tuple(25.0,    0.51,    0.096));
    data_cellSurvival_C4_2.push_back(make_tuple(50.0,    0.26,    0.035));
    data_cellSurvival_C4_2.push_back(make_tuple(75.0,    0.13,    0.031));
    data_cellSurvival_C4_2.push_back(make_tuple(100.0,    0.08,    0.035));
    data_cellSurvival_C4_2.push_back(make_tuple(150.0,    0.04,    0.024));

    CellSurvival cellSurvival_C4_2 = CellSurvival("C4_2");
    cellSurvival_C4_2.AddCellSurvivalData(data_cellSurvival_C4_2);
    // cellSurvival_C4_2.GraphCellSurvivalData();

    FitCellSurvival(cellSurvival_C4_2, "LM");

};