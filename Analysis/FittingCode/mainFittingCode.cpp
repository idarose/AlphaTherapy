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
#include "TGraph.h"
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
#include <vector>
#include <tuple>



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
        void LoadHistogramsFromAnalysis(CellSurvival cellSurvivalInstance, std::string regionName);
        std::vector<std::tuple<double, TH1D*, TH1D*>> Get_activity_energyDeps_histograms_vec(){return activity_energyDeps_histograms_vec;};
        std::vector<std::tuple<double, TH1D*>> Get_hitMultiplicity_histograms_vec(){return hitMultiplicity_histograms_vec;};
        std::vector<std::tuple<double, TH2D*>> Get_dose_hits_histograms_vec(){return dose_hits_histograms_vec;};

    private:
        double n;
        std::vector<std::tuple<double, TH1D*, TH1D*>> activity_energyDeps_histograms_vec;

        std::vector<std::tuple<double, TH1D*>> hitMultiplicity_histograms_vec;

        std::vector<std::tuple<double, TH2D*>> dose_hits_histograms_vec;
};


EnergyDepositionHistograms::EnergyDepositionHistograms(double t)
{
    n = t;
}


//----------------------------
void EnergyDepositionHistograms::LoadHistogramsFromAnalysis(CellSurvival cellSurvivalInstance, std::string regionName)
{
    std::vector<std::tuple<double,double,double>> cellSurvivalData = cellSurvivalInstance.GetCellSurvivalData();
    // General filepath to output
    std::string generalFilePath = "/Volumes/SamsungT7/OutputFromAnalysis/Output_" + cellSurvivalInstance.GetCellLine() + "_";

    // Specific filepath
    std::string filePath;

    // Looping through all activities
    for(int i=0; i<cellSurvivalData.size(); i++)
    {

        // Extracting activity
        double activity = std::get<0>(cellSurvivalData[i]);

        if(activity<=0.)
        {
            // Extracting activity
            double activity = std::get<0>(cellSurvivalData[i+1]);

            // Making filename
            filePath = generalFilePath + std::to_string(((int)activity)) + "kBq.root";

            // Read file containing histogram
            TFile* inputFile = new TFile(filePath.c_str(), "READ");

            //-------------------------------------
            // Extracting energy deposition histogram for nucleus with mGY binning
            TH1D* hEnergyDeps_Nucleus_keV = nullptr;
            std::string histogramName_keV = "i0_hDose_212Pb_" + cellSurvivalInstance.GetCellLine() + "_" + std::to_string(((int)activity)) +"kBq_" + regionName + "_mGyBinning";
            std::string histogramNewName = "i0_hDose_212Pb_" + cellSurvivalInstance.GetCellLine() + "_" + std::to_string(0) +"kBq_" + regionName + "_mGyBinning";


            inputFile->GetObject(histogramName_keV.c_str(), hEnergyDeps_Nucleus_keV);
            hEnergyDeps_Nucleus_keV->SetDirectory(0);
            hEnergyDeps_Nucleus_keV->SetName(histogramNewName.c_str());
            hEnergyDeps_Nucleus_keV->Reset();


            //-------------------------------------
            // Extracting energy deposition histogram for nucleus with uGyBinning
            TH1D* hEnergyDeps_Nucleus_eV = nullptr;
            std::string histogramName_eV = "i0_hDose_212Pb_" + cellSurvivalInstance.GetCellLine() + "_" + std::to_string(((int)activity)) +"kBq_" + regionName + "_uGyBinning";
            histogramNewName = "i0_hDose_212Pb_" + cellSurvivalInstance.GetCellLine() + "_" + std::to_string(0) +"kBq_" + regionName + "_uGyBinning";

            inputFile->GetObject(histogramName_eV.c_str(), hEnergyDeps_Nucleus_eV);
            hEnergyDeps_Nucleus_eV->SetDirectory(0);
            hEnergyDeps_Nucleus_eV->SetName(histogramNewName.c_str());
            hEnergyDeps_Nucleus_eV->Reset();


            //-------------------------------------
            // Extracting 2D energydeps vs number of hits histogram
            TH2D* hEnergyDeps_hitsAlpha_Nucleus = nullptr;
            std::string histogramName_energyDeps_hitsAlpha_Nucleus = "i0_hDose_212Pb_" + cellSurvivalInstance.GetCellLine() + "_" + std::to_string(((int)activity)) +"kBq_HitsAlpha_" + regionName;
            histogramNewName = "i0_hDose_212Pb_" + cellSurvivalInstance.GetCellLine() + "_" + std::to_string(0) +"kBq_HitsAlpha_" + regionName;

            // std::cout << histogramNewName << std::endl;

            inputFile->GetObject(histogramName_energyDeps_hitsAlpha_Nucleus.c_str(), hEnergyDeps_hitsAlpha_Nucleus);
            hEnergyDeps_hitsAlpha_Nucleus->SetDirectory(0);
            hEnergyDeps_hitsAlpha_Nucleus->SetName(histogramNewName.c_str());
            hEnergyDeps_hitsAlpha_Nucleus->Reset();


            //-------------------------------------
            // Extracting hit multiplicity histogram for nucleus
            TH1D* hHitMultiplicity_Nucleus = nullptr;
            std::string histogramName_hitMultiplicity = "i0_hFractionHitsAlpha_212Pb_" + cellSurvivalInstance.GetCellLine() + "_" + std::to_string(((int)activity)) +"kBq_" + regionName;

            inputFile->GetObject(histogramName_hitMultiplicity.c_str(), hHitMultiplicity_Nucleus);
            hHitMultiplicity_Nucleus->SetDirectory(0);
            hHitMultiplicity_Nucleus->SetName(histogramNewName.c_str());
            hHitMultiplicity_Nucleus->Reset();


            //-------------------------------------
            inputFile->Close();


            activity_energyDeps_histograms_vec.push_back(std::make_tuple(0.,hEnergyDeps_Nucleus_eV,hEnergyDeps_Nucleus_keV));
            hitMultiplicity_histograms_vec.push_back(std::make_tuple(0,hHitMultiplicity_Nucleus));
            dose_hits_histograms_vec.push_back(std::make_tuple(0.,hEnergyDeps_hitsAlpha_Nucleus));
        }

        else
        {
            // Making filename
            filePath = generalFilePath + std::to_string(((int)activity)) + "kBq.root";

            // Read file containing histogram
            TFile* inputFile = new TFile(filePath.c_str(), "READ");


            //-------------------------------------
            // Extracting energy deposition histogram for nucleus with keV binning
            TH1D* hEnergyDeps_Nucleus_keV = nullptr;
            std::string histogramName_keV = "i0_hDose_212Pb_" + cellSurvivalInstance.GetCellLine() + "_" + std::to_string(((int)activity)) +"kBq_" + regionName + "_mGyBinning";

            inputFile->GetObject(histogramName_keV.c_str(), hEnergyDeps_Nucleus_keV);
            hEnergyDeps_Nucleus_keV->SetDirectory(0);


            //-------------------------------------
            // Extracting energy deposition histogram for nucleus with eV binning
            TH1D* hEnergyDeps_Nucleus_eV = nullptr;
            std::string histogramName_eV = "i0_hDose_212Pb_" + cellSurvivalInstance.GetCellLine() + "_" + std::to_string(((int)activity)) +"kBq_" + regionName + "_uGyBinning";

            inputFile->GetObject(histogramName_eV.c_str(), hEnergyDeps_Nucleus_eV);
            hEnergyDeps_Nucleus_eV->SetDirectory(0);

            //-------------------------------------
            // Extracting 2D energydeps vs number of hits histogram
            TH2D* hEnergyDeps_hitsAlpha_Nucleus = nullptr;
            std::string histogramName_energyDeps_hitsAlpha_Nucleus = "i0_hDose_212Pb_" + cellSurvivalInstance.GetCellLine() + "_" + std::to_string(((int)activity)) +"kBq_HitsAlpha_" + regionName;

            inputFile->GetObject(histogramName_energyDeps_hitsAlpha_Nucleus.c_str(), hEnergyDeps_hitsAlpha_Nucleus);
            hEnergyDeps_hitsAlpha_Nucleus->SetDirectory(0);

            //-------------------------------------
            // Extracting hit multiplicity histogram for nucleus
            TH1D* hHitMultiplicity_Nucleus = nullptr;
            std::string histogramName_hitMultiplicity = "i0_hFractionHitsAlpha_212Pb_" + cellSurvivalInstance.GetCellLine() + "_" + std::to_string(((int)activity)) +"kBq_" + regionName;

            inputFile->GetObject(histogramName_hitMultiplicity.c_str(), hHitMultiplicity_Nucleus);
            hHitMultiplicity_Nucleus->SetDirectory(0);


            //-------------------------------------
            inputFile->Close();


            activity_energyDeps_histograms_vec.push_back(std::make_tuple(activity,hEnergyDeps_Nucleus_eV,hEnergyDeps_Nucleus_keV));
            hitMultiplicity_histograms_vec.push_back(std::make_tuple(activity,hHitMultiplicity_Nucleus));
            dose_hits_histograms_vec.push_back(std::make_tuple(activity,hEnergyDeps_hitsAlpha_Nucleus));
        }
    }

    // for(auto & entry : activity_energyDeps_histograms_vec)
    // {
    //     double activity = std::get<0>(entry);
    //     auto hist = std::get<2>(entry);

    //     double integral = hist->Integral();

    //     std::cout << "activity: " << activity << ",\t integral: " << integral << std::endl;
    // }
}


void FitCellSurvival(CellSurvival cellSurvivalInstance, std::string modelName, int region)
{
    //------------------------
    // Defininng volume used for fit
    std::string regionName;
    if(region==1)
    {
        // Fitting for membrane
        regionName = "Membrane";
    }
    else if(region==2)
    {
        //Fitting for cytoplasm
        regionName = "Cytoplasm";
    }
    else if(region==3)
    {
        // Fitting for nucleus
        regionName = "Nucleus";
    }
    else if(region==4)
    {
        regionName = "TotalCell";
    }


    //--------------------------
    // Getting and graphing cell survival data
    std::vector<std::tuple<double,double,double>> data_cellSurvival;
    data_cellSurvival = cellSurvivalInstance.GetCellSurvivalData();


    TGraphErrors *gr_clonogenicSurvival = new TGraphErrors();
    gr_clonogenicSurvival->SetName("gr_clonogenicSurvival");

    auto GraphData = [&](TGraphErrors* gr_clonogenicSurvival_CellLine, std::vector<std::tuple<double, double, double>> data_cellSurvival_CellLine)
    {
        int graphPointN_data = 0;
        for(auto & entry : data_cellSurvival_CellLine)
        {
            double activity = std::get<0>(entry);
            double cellSurvival = std::get<1>(entry);
            double cellSurvivalUncertainty = std::get<2>(entry);

            gr_clonogenicSurvival_CellLine->SetPoint(graphPointN_data, activity, cellSurvival);
            gr_clonogenicSurvival_CellLine->SetPointError(graphPointN_data, 0.0, cellSurvivalUncertainty);

            graphPointN_data++;
        }
    };

    GraphData(gr_clonogenicSurvival, data_cellSurvival);


    //------------------------------
    // Loading energy deposition histograms
    EnergyDepositionHistograms energyDepHistograms = EnergyDepositionHistograms(1.0);
    energyDepHistograms.LoadHistogramsFromAnalysis(cellSurvivalInstance, regionName);


    std::vector<std::tuple<double,TH1D*,TH1D*>> vec_activities_histograms = energyDepHistograms.Get_activity_energyDeps_histograms_vec();

    std::vector<std::tuple<double, TH1D*>> hitMultiplicity_histograms_vec = energyDepHistograms.Get_hitMultiplicity_histograms_vec();

    std::vector<std::tuple<double, TH2D*>> dose_hits_histograms_vec = energyDepHistograms.Get_dose_hits_histograms_vec();

    //------------------------------
    // Defining decay dynamics to calculate scaling factor for dose histograms
    double VolumeSample = 0.2*1000; // mm^3
    double volumeCellTube = TMath::Pi()*std::pow(0.5,2.0)*1.0; // mm^3
    double volumeRatio = volumeCellTube/VolumeSample;

    double numberIterations =1.;
    double numberCells = 500000.*volumeRatio;
    double scalingFactor = numberIterations*numberCells;


    //------------------------------
    auto CalculateCellSurvivalFraction = [&](TH1D *h_energyDeposition_eV, TH1D *h_energyDeposition_keV, double *par)
    {
        double alpha = par[0];
        double beta = par[1];

        // std::cout << "alpha: " << alpha << std::endl;
        // std::cout << "beta: " << beta << std::endl;

        //------------------------------------------------
        double fractionOfComponentsHit = h_energyDeposition_keV->Integral();
        double fractionOfComponentsMissed = (1.0 - fractionOfComponentsHit);

        //------------------------------------------------
        //      The fraction of missed cells obviously all survive so that is immediately added.
        double fractionOfTotalCellsHitSurviving = fractionOfComponentsMissed;

        // Uncertainty in cell survival
        std::vector<double> deltaFractionOfTotalCellsSurviving_Vec;

        double deltaAlpha = 0.05*alpha;


        // Looping over uGy binned histogram
        for(int i=0; i<h_energyDeposition_eV->GetNbinsX();i++)
        {
            double energyDeposition = h_energyDeposition_eV->GetBinCenter(i+1);
            double fractionOfTotalCellsHit = h_energyDeposition_eV->GetBinContent(i+1);

            // Uncertainty in fraction of cells hit at this dose
            double deltaFractionOfTotalCellsHit = (1./scalingFactor)*std::sqrt(scalingFactor*fractionOfTotalCellsHit);

            //----------------------------
            double cellSurvivalfraction = TMath::Exp(-(alpha*energyDeposition + beta*TMath::Power(energyDeposition, 2.0)));

            double fractionOfTotalCellsHit_survivingFraction = cellSurvivalfraction*fractionOfTotalCellsHit;

            // Uncertainty in surviving fraction of fraction of cells hit
            double deltaFractionOfTotalCellsHit_survivingFraction = fractionOfTotalCellsHit_survivingFraction*std::sqrt( std::pow(energyDeposition*deltaAlpha,2.) + std::pow((1./fractionOfTotalCellsHit)*deltaFractionOfTotalCellsHit,2.));

            // Adding uncertainty to storage vector
            deltaFractionOfTotalCellsSurviving_Vec.push_back(deltaFractionOfTotalCellsHit_survivingFraction);

            fractionOfTotalCellsHitSurviving += fractionOfTotalCellsHit_survivingFraction;
        }

        // Looping over mGy binned histogram
        for(int i=1000; i<h_energyDeposition_keV->GetNbinsX(); i++)
        {
            double energyDeposition = h_energyDeposition_keV->GetBinCenter(i+1);
            double fractionOfTotalCellsHit = h_energyDeposition_keV->GetBinContent(i+1);

            // Uncertainty in fraction of cells hit at this dose
            double deltaFractionOfTotalCellsHit = (1./scalingFactor)*std::sqrt(scalingFactor*fractionOfTotalCellsHit);

            //----------------------------
            double cellSurvivalfraction = TMath::Exp(-(alpha*energyDeposition + beta*TMath::Power(energyDeposition, 2.0)));

            double fractionOfTotalCellsHit_survivingFraction = cellSurvivalfraction*fractionOfTotalCellsHit;

            // Uncertainty in surviving fraction of fraction of cells hit
            double deltaFractionOfTotalCellsHit_survivingFraction = fractionOfTotalCellsHit_survivingFraction*std::sqrt( std::pow(energyDeposition*deltaAlpha,2.) + std::pow((1./fractionOfTotalCellsHit)*deltaFractionOfTotalCellsHit,2.));

            // Adding uncertainty to storage vector
            deltaFractionOfTotalCellsSurviving_Vec.push_back(deltaFractionOfTotalCellsHit_survivingFraction);

            fractionOfTotalCellsHitSurviving += fractionOfTotalCellsHit_survivingFraction;
        }

        double varianceFractionOfTotalCellsSurviving = 0.;
        for(int i=0; i<deltaFractionOfTotalCellsSurviving_Vec.size(); i++)
        {
            varianceFractionOfTotalCellsSurviving += std::pow(deltaFractionOfTotalCellsSurviving_Vec[i],2.);
        }

        double deltaFractionOfTotalCellsSurviving = std::sqrt(varianceFractionOfTotalCellsSurviving);
        std::tuple<double,double> fractionSurvived_deltaFractionSurived = std::make_tuple(fractionOfTotalCellsHitSurviving, deltaFractionOfTotalCellsSurviving);
        //------------------------------------------------
        return fractionSurvived_deltaFractionSurived;
    };

    //------------------------------------------------------------------------------------------------
    //      This generates a new graph corresponding to the predicted survivability.
    //      This graph is calculated everytime the parameters (par) are changed.
    auto GenerateGraph_CellSurvivalFraction = [&](std::vector<std::tuple<double,TH1D*,TH1D*>> vec, double *par)
    {
        TGraphErrors gr;

        int graphPointN = 0;

        for(auto & entry : vec)
        {
            double activity = std::get<0>(entry);
            auto hist_eV = std::get<1>(entry);
            auto hist_keV = std::get<2>(entry);

            double cellSurvival = std::get<0>(CalculateCellSurvivalFraction(hist_eV, hist_keV, par));
            double deltaCellSurvival = std::get<1>(CalculateCellSurvivalFraction(hist_eV, hist_keV, par));

            gr.SetPoint(graphPointN, activity, cellSurvival);
            gr.SetPointError(graphPointN, 0.0, deltaCellSurvival);

            graphPointN++;
        }

        return gr;
    };

     //------------------------------------------------------------------------------------------------
    int nParameters = 2;
    double savedParameters[nParameters];
    TGraphErrors gr_cellSurvivability_vs_activitykBqPerMl;

    //------------------------------------------------------------------------------------------------
    auto f_cellSurvivalVsDose_C4_2 = new TF1("f_cellSurvivalVsDose_C4_2",
        [&](double*x, double *p)
        {
            double activity_Bq = x[0];

            //------------------------------------------------
            //      This block of code is to test whether the parameters have changed.
            //      A new graph of survivability is generated only if the parameters are found to have changed.
            //      This makes the code much faster, which is useful if your calculation for cell survivability becomes more computationally heavy
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
                gr_cellSurvivability_vs_activitykBqPerMl = GenerateGraph_CellSurvivalFraction(vec_activities_histograms, p);
            }

            //------------------------------------------------
            double cellSurvival = gr_cellSurvivability_vs_activitykBqPerMl.Eval(x[0], 0, "");

            return cellSurvival;

        }, 0.0, 150.0, nParameters);

    f_cellSurvivalVsDose_C4_2->SetParLimits(0, 0.0, 1.0e+05);
    f_cellSurvivalVsDose_C4_2->SetParLimits(1, 0.0, 1.0e+05);

    f_cellSurvivalVsDose_C4_2->SetNpx(10000);
    f_cellSurvivalVsDose_C4_2->SetParameter(0, 1.0e+00);
    // f_cellSurvivalVsDose_C4_2->SetParameter(1, 1.0e+00);

    if(modelName=="LM")
    {
        f_cellSurvivalVsDose_C4_2->FixParameter(1, 0.0e+00);
    }
    if(modelName=="LQ")
    {
        f_cellSurvivalVsDose_C4_2->SetParameter(1, 1.0e+00);
    }

    std::cout << " Cell Line " + cellSurvivalInstance.GetCellLine() << " Model : " << modelName << " Region : " << regionName << std::endl;
    gr_clonogenicSurvival->Fit(f_cellSurvivalVsDose_C4_2, "", "", 0.0, 150.0);

    double chi_sq = f_cellSurvivalVsDose_C4_2->GetChisquare();

    int deg_freedom = f_cellSurvivalVsDose_C4_2->GetNDF();

    std::cout << "Chi squared: " << chi_sq << std::endl;
    std::cout << "Deg Freedom: " << deg_freedom << std::endl;
    std::cout << "Reduced Chi squared: " << chi_sq/deg_freedom << std::endl;



    //------------------------------------------------------------------------------------------------
    std::string outputName = "Output_AnalyseCellSurvival_" + cellSurvivalInstance.GetCellLine() + "_" + regionName + "_" + modelName + ".root";
    TFile *outputFile = new TFile(outputName.c_str(), "RECREATE");

    std::string titleGraph = "Cell Survival Using Energy Deposition in " + regionName;
    gr_clonogenicSurvival->SetTitle(titleGraph.c_str());
    gr_clonogenicSurvival->GetXaxis()->SetTitle("Activity [kBq/mL]");
    gr_clonogenicSurvival->GetYaxis()->SetTitle("Survival Fraction");
    gr_clonogenicSurvival->Write();
    f_cellSurvivalVsDose_C4_2->Write();
    gr_cellSurvivability_vs_activitykBqPerMl.Write();

    auto legend = new TLegend(0.1,0.7,0.48,0.9);
    legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
    legend->AddEntry(gr_clonogenicSurvival,"Experimental Data","f");
    legend->AddEntry(f_cellSurvivalVsDose_C4_2,"Simulated Data","f");

    legend->Write();


    for(auto & entry : vec_activities_histograms)
    {
        auto hist = std::get<2>(entry);

        hist->Write();
    }

    double alpha = savedParameters[0];


    //------------------------
    // Hit Analysis

    std::vector<std::tuple<double,double>> percentKilledNumberHits[11];


    // Tuple< nHits , Vector< activity, mean dose >>
    // std::tuple<int,std::vector<std::tuple<double,double>>> meandDose_Activity_PerNHits;

    std::vector< std::tuple<double,double>> meanDose_PerHits[11];

    for(int i=0; i<dose_hits_histograms_vec.size(); i++)
    {
        //------------------------
        // Defining activity
        double activity = std::get<0>(dose_hits_histograms_vec[i]);

        // Extracting dose per hits histograms
        TH2D* dose_hits_histogram = std::get<1>(dose_hits_histograms_vec[i]);

        // Define number of bins and max number of hits
        int nBins = dose_hits_histogram->GetNbinsY();
        double hitsMax = 700.;


        //-------------------------
        // Making histogram for survical fraction for number of hits
        std::string hitMultiplicity_Survival_Name = "hHitMultiplicity_" + regionName + "_Survival_" + std::to_string((int)activity) + "kBq";
        TH1D* hitMultiplicity_Survival_Histogram = new TH1D(hitMultiplicity_Survival_Name.c_str(), "Cell Survival per N Number of Hits by Alpha Particle", nBins, 0,hitsMax);

        //--------------------------
        // Making histogram for percent of cells killed for number of hist
        std::string hitMultiplicity_Survival_Percentage_Name = "hHitMultiplicity_" + regionName + "_Death_Percentage_" + std::to_string((int)activity) + "kBq";
        TH1D* hitMultiplicity_Survival_Percentage_Histogram = new TH1D(hitMultiplicity_Survival_Percentage_Name.c_str(), "Cell Death per N Number of Hits by Alpha Particle", nBins, 0,hitsMax);


        //---------------------------
        // Vector to store fraction of cells survived at each number of hits
        std::vector<std::tuple<int,double>> hits_SurvivalFraction_Vec;

        // Loop over y axis (number hits)
        for(int i=0; i<dose_hits_histogram->GetNbinsY(); i++)
        {
            double survivalFractionThisHitNumber = 0.;

            // Loop over x axis (dose deposited)
            for(int j=0; j<dose_hits_histogram->GetNbinsX(); j++)
            {
                double doseDep = (dose_hits_histogram->GetXaxis())->GetBinCenter(j+1);
                double survivalThisEnergy = dose_hits_histogram->GetBinContent(j+1,i+1)*TMath::Exp(-alpha*doseDep);

                survivalFractionThisHitNumber += survivalThisEnergy;
            }

            std::string nameProjectedHist = "doseDelivered_NHits_" + std::to_string(i) + "_" + std::to_string((int)activity) + "kBq";
            TH1D* projectedDose_ForNHits = dose_hits_histogram->ProjectionX(nameProjectedHist.c_str(), i,i+1);
            double meanDose = projectedDose_ForNHits->GetMean(1);
            projectedDose_ForNHits->SetDirectory(0);

            if(i<11)
            {
                meanDose_PerHits[i].push_back(std::make_tuple(activity,meanDose));
            }

            hits_SurvivalFraction_Vec.push_back(std::make_tuple(i,survivalFractionThisHitNumber));
        }


        //---------------------------
        // Add survival fraction per number of hits to histogram
        for(int i=0; i<hits_SurvivalFraction_Vec.size(); i++)
        {
            int hits = std::get<0>(hits_SurvivalFraction_Vec[i]);
            double fraction = std::get<1>(hits_SurvivalFraction_Vec[i]);


            if(fraction>=0.)
            {
                hitMultiplicity_Survival_Histogram->SetBinContent(hits+1, fraction);
            }
        }



        //---------------------------
        // Calculating percent of cells dead at each number of hits

        // Getting histogram with fraction of cells hit per cell hit
        auto hitMultiplicityHist  = std::get<1>(hitMultiplicity_histograms_vec[i]);

        // Bool for finding 99% dead
        // bool found90PercentPoint = false;

        // Looping over number of hits
        for(int i=0; i<hitMultiplicityHist->GetNbinsX(); i++)
        {
            // Fraction of cells hit a number of times
            double fractionHit = hitMultiplicityHist->GetBinContent(i+1);

            // Fraction of cells hit a number of times surviving
            double fractionSurvived = hitMultiplicity_Survival_Histogram->GetBinContent(i+1);

            if(fractionHit>0.)
            {
                double percentKilled = 100. - 100.*fractionSurvived/fractionHit;

                // Adding percent killed to histogram
                hitMultiplicity_Survival_Percentage_Histogram->SetBinContent(i+1, percentKilled);

                // Storing percentage killed, separated by number of hits
                if(i<11)
                {
                    percentKilledNumberHits[i].push_back(std::make_tuple(activity,percentKilled));
                }
            }
        }

        hitMultiplicity_Survival_Histogram->GetXaxis()->SetTitle("N Number of Hits to Cell Component by Alpha Particle");
        hitMultiplicity_Survival_Histogram->GetYaxis()->SetTitle("Fraction of Cells Hit");

        hitMultiplicity_Survival_Percentage_Histogram->GetXaxis()->SetTitle("N Number of Hits to Cell Component by Alpha Particle");
        hitMultiplicity_Survival_Percentage_Histogram->GetYaxis()->SetTitle("Percentage of Cells Killed");


        // hitMultiplicityHist->Write();
    }


    auto GraphDataPerNumberHits = [&](TGraph* gr_DataPerNumberHits, std::vector<std::tuple<double, double>> dataPerNumberHits_vec)
    {
        int graphPointN_data = 0;
        for(auto & entry : dataPerNumberHits_vec)
        {
            double activity = std::get<0>(entry);
            double percentKilled = std::get<1>(entry);

            gr_DataPerNumberHits->SetPoint(graphPointN_data, activity, percentKilled);
            graphPointN_data++;
        }

    };

    for(int i=0; i<11; i++)
    {
        std::string graphName_percentKilled = "gr_percentKilled_" + std::to_string(i) + "_NumberHits";
        TGraph * gr_percentKilled_i = new TGraph();
        gr_percentKilled_i->SetName(graphName_percentKilled.c_str());
        GraphDataPerNumberHits(gr_percentKilled_i, percentKilledNumberHits[i]);
        gr_percentKilled_i->GetXaxis()->SetTitle("Activity [kBq/mL]");
        gr_percentKilled_i->GetYaxis()->SetTitle("Percent of Cells Killed");
        gr_percentKilled_i->Write();

        std::string graphName_meanDose = "gr_meanDose_" + std::to_string(i) + "_NumberHits";
        TGraph* gr_meanDose_i = new TGraph();
        gr_meanDose_i->SetName(graphName_meanDose.c_str());
        GraphDataPerNumberHits(gr_meanDose_i, meanDose_PerHits[i]);
        gr_meanDose_i->GetXaxis()->SetTitle("Activity [kBq/mL]");
        gr_meanDose_i->GetYaxis()->SetTitle("Mean Dose Delivered [Gy]");
        gr_meanDose_i->Write();

    }


    // gr_PercentKilledOneHits->Write();
    // gr_PercentKilledTwoHits->Write();
    // gr_PercentKilledThreeHits->Write();

    outputFile->Write();
    outputFile->Close();


}



//----------------------------
void mainFittingCode()
{
    // std::vector<std::tuple<double,double,double>> data_cellSurvival_PC3_PIP;

    // // data_cellSurvival_PC3_PIP.push_back(make_tuple(0.0,    1.0,    0.05));
    // data_cellSurvival_PC3_PIP.push_back(make_tuple(10.0,    0.630,    std::sqrt(std::pow(0.063,2.) + std::pow(0.1*0.630,2.))));
    // data_cellSurvival_PC3_PIP.push_back(make_tuple(25.0,    0.317,    std::sqrt(std::pow(0.05,2.) + std::pow(0.1*0.317,2.))));
    // data_cellSurvival_PC3_PIP.push_back(make_tuple(50.0,    0.071,    std::sqrt(std::pow(0.05,2.) + std::pow(0.1*0.071,2.))));
    // data_cellSurvival_PC3_PIP.push_back(make_tuple(75.0,    0.032,    std::sqrt(std::pow(0.05,2.) + std::pow(0.1*0.032,2.))));
    // data_cellSurvival_PC3_PIP.push_back(make_tuple(100.0,    0.014,    std::sqrt(std::pow(0.05,2.) + std::pow(0.1*0.014,2.))));
    // data_cellSurvival_PC3_PIP.push_back(make_tuple(150.0,    0.006,    std::sqrt(std::pow(0.05,2.) + std::pow(0.1*0.006,2.))));


    // CellSurvival cellSurvival_PC3_PIP = CellSurvival("PC3_PIP");
    // cellSurvival_PC3_PIP.AddCellSurvivalData(data_cellSurvival_PC3_PIP);
    // FitCellSurvival(cellSurvival_PC3_PIP, "LM", 1);
    // FitCellSurvival(cellSurvival_PC3_PIP, "LM", 2);
    // FitCellSurvival(cellSurvival_PC3_PIP, "LM", 3);
    // FitCellSurvival(cellSurvival_PC3_PIP, "LM", 4);


    // std::vector<std::tuple<double,double,double>> data_cellSurvival_PC3_Flu;
    // // data_cellSurvival_PC3_Flu.push_back(make_tuple(0.0,    1.0,    0.05));
    // data_cellSurvival_PC3_Flu.push_back(make_tuple(10.0,    0.955,    std::sqrt(std::pow(0.0955,2.) + std::pow(0.1*0.955,2.))));
    // data_cellSurvival_PC3_Flu.push_back(make_tuple(25.0,    0.724,    std::sqrt(std::pow(0.0724,2.) + std::pow(0.1*0.724,2.))));
    // data_cellSurvival_PC3_Flu.push_back(make_tuple(50.0,    0.733,    std::sqrt(std::pow(0.0733,2.) + std::pow(0.1*0.733,2.))));
    // data_cellSurvival_PC3_Flu.push_back(make_tuple(75.0,    0.798,    std::sqrt(std::pow(0.0798,2.) + std::pow(0.1*0.798,2.))));
    // data_cellSurvival_PC3_Flu.push_back(make_tuple(100.0,    0.729,    std::sqrt(std::pow(0.0720,2.) + std::pow(0.1*0.729,2.))));
    // data_cellSurvival_PC3_Flu.push_back(make_tuple(150.0,    0.690,    std::sqrt(std::pow(0.0698,2.) + std::pow(0.1*0.690,2.))));


    // CellSurvival cellSurvival_PC3_Flu = CellSurvival("PC3_Flu");
    // cellSurvival_PC3_Flu.AddCellSurvivalData(data_cellSurvival_PC3_Flu);
    // FitCellSurvival(cellSurvival_PC3_Flu, "LM", 1);
    // FitCellSurvival(cellSurvival_PC3_Flu, "LM", 2);
    // FitCellSurvival(cellSurvival_PC3_Flu, "LM", 3);
    //  FitCellSurvival(cellSurvival_PC3_Flu, "LM", 4);

    std::vector<std::tuple<double,double,double>> data_cellSurvival_C4_2;
    // data_cellSurvival_C4_2.push_back(make_tuple(0.0,    1.0,    0.0));
    data_cellSurvival_C4_2.push_back(make_tuple(5.0,    0.86,   std::sqrt(std::pow(0.117,2.) + std::pow(0.1*0.86,2.))));
    data_cellSurvival_C4_2.push_back(make_tuple(10.0,    0.69,    std::sqrt(std::pow(0.014,2.0) + std::pow(0.1*0.69,2.))));
    data_cellSurvival_C4_2.push_back(make_tuple(25.0,    0.51,    std::sqrt(std::pow(0.096,2.) + std::pow(0.1*0.51,2.))));
    data_cellSurvival_C4_2.push_back(make_tuple(50.0,    0.26,   std::sqrt(std::pow(0.035,2.) + std::pow(0.1*0.26,2.))));
    data_cellSurvival_C4_2.push_back(make_tuple(75.0,    0.13,    std::sqrt(std::pow(0.031,2.) + std::pow(0.1*0.13,2.))));
    data_cellSurvival_C4_2.push_back(make_tuple(100.0,    0.08,    std::sqrt(std::pow(0.035,2.) + std::pow(0.1*0.08,2.))));
    data_cellSurvival_C4_2.push_back(make_tuple(150.0,    0.04,    std::sqrt(std::pow(0.024,2.) + std::pow(0.1*0.04,2.))));

    CellSurvival cellSurvival_C4_2 = CellSurvival("C4_2");
    cellSurvival_C4_2.AddCellSurvivalData(data_cellSurvival_C4_2);
    // cellSurvival_C4_2.GraphCellSurvivalData();

    // // std::vector<double,double> C4_2_Parameters;
    // FitCellSurvival(cellSurvival_C4_2, "LM", 1);
    // FitCellSurvival(cellSurvival_C4_2, "LM", 2);
    FitCellSurvival(cellSurvival_C4_2, "LM", 3);
    FitCellSurvival(cellSurvival_C4_2, "LM", 4);

};