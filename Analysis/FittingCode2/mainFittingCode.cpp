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
        std::vector<std::tuple<double, TH2D*>> Get_energyDeps_hits_histograms_vec(){return energyDeps_hits_histograms_vec;};

    private:
        double n;
        std::vector<std::tuple<double, TH1D*, TH1D*>> activity_energyDeps_histograms_vec;

        std::vector<std::tuple<double, TH1D*>> hitMultiplicity_histograms_vec;

        std::vector<std::tuple<double, TH2D*>> energyDeps_hits_histograms_vec;
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
            energyDeps_hits_histograms_vec.push_back(std::make_tuple(0.,hEnergyDeps_hitsAlpha_Nucleus));
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
            energyDeps_hits_histograms_vec.push_back(std::make_tuple(activity,hEnergyDeps_hitsAlpha_Nucleus));
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

    std::vector<std::tuple<double, TH2D*>> energyDeps_hits_histograms_vec = energyDepHistograms.Get_energyDeps_hits_histograms_vec();


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
        double fractionOfTotalCellsSurviving = fractionOfComponentsMissed;

        for(int i=0; i<h_energyDeposition_eV->GetNbinsX();i++)
        {
            double energyDeposition = h_energyDeposition_eV->GetBinCenter(i+1);
            double fractionOfTotalCells = h_energyDeposition_eV->GetBinContent(i+1);

            //----------------------------
            double cellSurvivalfraction = TMath::Exp(-(alpha*energyDeposition + beta*TMath::Power(energyDeposition, 2.0)));

            double fractionOfTotalCells_survivingFraction = cellSurvivalfraction*fractionOfTotalCells;

            fractionOfTotalCellsSurviving += fractionOfTotalCells_survivingFraction;
        }

        for(int i=1000; i<h_energyDeposition_keV->GetNbinsX(); i++)
        {
            double energyDeposition = h_energyDeposition_keV->GetBinCenter(i+1);
            double fractionOfTotalCells = h_energyDeposition_keV->GetBinContent(i+1);

            //----------------------------
            double cellSurvivalfraction = TMath::Exp(-(alpha*energyDeposition + beta*TMath::Power(energyDeposition, 2.0)));

            double fractionOfTotalCells_survivingFraction = cellSurvivalfraction*fractionOfTotalCells;

            fractionOfTotalCellsSurviving += fractionOfTotalCells_survivingFraction;
        }

        //------------------------------------------------
        return fractionOfTotalCellsSurviving;
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

            double cellSurvival = CalculateCellSurvivalFraction(hist_eV, hist_keV, par);

            gr.SetPoint(graphPointN, activity, cellSurvival);

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

    std::vector<std::tuple<double,double>> percentKilledOneHit;
    std::vector<std::tuple<double,double>> percentKilledTwoHit;
    std::vector<std::tuple<double,double>> percentKilledThreeHit;




    for(int i=0; i<energyDeps_hits_histograms_vec.size(); i++)
    {
        double activity = std::get<0>(energyDeps_hits_histograms_vec[i]);
        TH2D* energyDeps_hits_histogram = std::get<1>(energyDeps_hits_histograms_vec[i]);

        int nBins = energyDeps_hits_histogram->GetNbinsY();
        double hitsMax = 700.;

        std::string hitMultiplicity_Survival_Name = "hHitMultiplicity_" + regionName + "_Survival_" + std::to_string((int)activity) + "kBq";
        TH1D* hitMultiplicity_Survival_Histogram = new TH1D(hitMultiplicity_Survival_Name.c_str(), "Cell Survival per N Number of Hits by Alpha Particle", nBins, 0,hitsMax);

        std::string hitMultiplicity_Survival_Percentage_Name = "hHitMultiplicity_" + regionName + "_Death_Percentage_" + std::to_string((int)activity) + "kBq";
        TH1D* hitMultiplicity_Survival_Percentage_Histogram = new TH1D(hitMultiplicity_Survival_Percentage_Name.c_str(), "Cell Death per N Number of Hits by Alpha Particle", nBins, 0,hitsMax);

        std::vector<std::tuple<int,double>> hitSurvivalVec;

        // Loop over y axis (number hits)
        for(int i=0; i<energyDeps_hits_histogram->GetNbinsY(); i++)
        {
            double survivalFractionThisHitNumber = 0.;

            for(int j=0; j<energyDeps_hits_histogram->GetNbinsX(); j++)
            {
                double energyDep = (energyDeps_hits_histogram->GetXaxis())->GetBinCenter(j+1);
                double survivalThisEnergy = energyDeps_hits_histogram->GetBinContent(j+1,i+1)*TMath::Exp(-alpha*energyDep);

                survivalFractionThisHitNumber += survivalThisEnergy;
            }

            hitSurvivalVec.push_back(std::make_tuple(i,survivalFractionThisHitNumber));
        }


        for(int i=0; i<hitSurvivalVec.size(); i++)
        {
            int hits = std::get<0>(hitSurvivalVec[i]);
            double fraction = std::get<1>(hitSurvivalVec[i]);


            if(fraction>=0.)
            {
                hitMultiplicity_Survival_Histogram->SetBinContent(hits+1, fraction);
            }
        }

        // activity = std::get<0>(hitMultiplicity_histograms_vec[i]);
        auto hitMultiplicityHist  = std::get<1>(hitMultiplicity_histograms_vec[i]);

        bool found90PercentPoint = false;

        for(int i=0; i<hitMultiplicityHist->GetNbinsX(); i++)
        {
            double fractionHit = hitMultiplicityHist->GetBinContent(i+1);
            double fractionSurvived = hitMultiplicity_Survival_Histogram->GetBinContent(i+1);

            if(fractionHit>0.)
            {
                double percentKilled = 100. - 100.*fractionSurvived/fractionHit;
                hitMultiplicity_Survival_Percentage_Histogram->SetBinContent(i+1, percentKilled);


                if(percentKilled>99.)
                {
                    if(!found90PercentPoint){
                        std::cout << "For " << activity << "kBq 99p cell death at hit number " << i+1 << std::endl;
                    }
                    found90PercentPoint = true;
                }

                hitMultiplicity_Survival_Percentage_Histogram->SetBinContent(i+1, percentKilled);
                if(i==1)
                {
                    percentKilledOneHit.push_back(std::make_tuple(activity,percentKilled));
                    std::cout << "Activity : " << activity << " Num : " << i << " Percent Killed : " << percentKilled << std::endl;
                }
                else if(i==2)
                {
                    percentKilledTwoHit.push_back(std::make_tuple(activity,percentKilled));
                    std::cout << "Activity : " << activity << " Num : " << i << " Percent Killed : " << percentKilled << std::endl;
                }
                else if(i==3)
                {
                    percentKilledThreeHit.push_back(std::make_tuple(activity,percentKilled));
                    std::cout << "Activity : " << activity << " Num : " << i << " Percent Killed : " << percentKilled << std::endl;
                }
            }
        }

        hitMultiplicity_Survival_Histogram->GetXaxis()->SetTitle("N Number of Hits to Cell Component by Alpha Particle");
        hitMultiplicity_Survival_Histogram->GetYaxis()->SetTitle("Fraction of Cells Hit");

        hitMultiplicity_Survival_Percentage_Histogram->GetXaxis()->SetTitle("N Number of Hits to Cell Component by Alpha Particle");
        hitMultiplicity_Survival_Percentage_Histogram->GetYaxis()->SetTitle("Percentage of Cells Killed");


        // hitMultiplicityHist->Write();
    }

    TGraph *gr_PercentKilledOneHits = new TGraph();
    gr_PercentKilledOneHits->SetName("gr_PercentKilledOneHit");

    TGraph *gr_PercentKilledTwoHits = new TGraph();
    gr_PercentKilledTwoHits->SetName("gr_PercentKilledTwoHits");

    TGraph *gr_PercentKilledThreeHits = new TGraph();
    gr_PercentKilledThreeHits->SetName("gr_PercentKilledThreeHits");

    auto Graph_percentKilledNumberHits = [&](TGraph* gr_percentKilledNumberHits, std::vector<std::tuple<double, double>> percentKilledNumberHits_vec)
    {
        int graphPointN_data = 0;
        for(auto & entry : percentKilledNumberHits_vec)
        {
            double activity = std::get<0>(entry);
            double percentKilled = std::get<1>(entry);

            gr_percentKilledNumberHits->SetPoint(graphPointN_data, activity, percentKilled);
            graphPointN_data++;
        }
    };

    Graph_percentKilledNumberHits(gr_PercentKilledOneHits, percentKilledOneHit);
    Graph_percentKilledNumberHits(gr_PercentKilledTwoHits, percentKilledTwoHit);
    Graph_percentKilledNumberHits(gr_PercentKilledThreeHits, percentKilledThreeHit);

    // TGraph* gr_PercentKilledOneHit = GenerateGraph_PercentKilledForNumberHit(percentKilledOneHit);
    // TGraph* gr_PercentKilledTwoHit = GenerateGraph_PercentKilledForNumberHit(percentKilledTwoHit);
    // TGraph* gr_PercentKilledThreeHit = GenerateGraph_PercentKilledForNumberHit(percentKilledThreeHit);


    gr_PercentKilledOneHits->Write();
    gr_PercentKilledTwoHits->Write();
    gr_PercentKilledThreeHits->Write();

    outputFile->Write();
    outputFile->Close();

}



//----------------------------
void mainFittingCode()
{
    // std::vector<std::tuple<double,double,double>> data_cellSurvival_PC3_PIP;

    // data_cellSurvival_PC3_PIP.push_back(make_tuple(0.0,    1.0,    0.05));
    // data_cellSurvival_PC3_PIP.push_back(make_tuple(10.0,    0.630,    0.063));
    // data_cellSurvival_PC3_PIP.push_back(make_tuple(25.0,    0.317,    0.05));
    // data_cellSurvival_PC3_PIP.push_back(make_tuple(50.0,    0.071,    0.05));
    // data_cellSurvival_PC3_PIP.push_back(make_tuple(75.0,    0.032,    0.05));
    // data_cellSurvival_PC3_PIP.push_back(make_tuple(100.0,    0.014,    0.05));
    // data_cellSurvival_PC3_PIP.push_back(make_tuple(150.0,    0.006,    0.05));


    // CellSurvival cellSurvival_PC3_PIP = CellSurvival("PC3_PIP");
    // cellSurvival_PC3_PIP.AddCellSurvivalData(data_cellSurvival_PC3_PIP);
    // FitCellSurvival(cellSurvival_PC3_PIP, "LM", 1);
    // FitCellSurvival(cellSurvival_PC3_PIP, "LM", 2);
    // FitCellSurvival(cellSurvival_PC3_PIP, "LM", 3);


    // std::vector<std::tuple<double,double,double>> data_cellSurvival_PC3_Flu;
    // data_cellSurvival_PC3_Flu.push_back(make_tuple(0.0,    1.0,    0.05));
    // data_cellSurvival_PC3_Flu.push_back(make_tuple(10.0,    0.955,    0.0955));
    // data_cellSurvival_PC3_Flu.push_back(make_tuple(25.0,    0.724,    0.0724));
    // data_cellSurvival_PC3_Flu.push_back(make_tuple(50.0,    0.733,    0.0733));
    // data_cellSurvival_PC3_Flu.push_back(make_tuple(75.0,    0.798,    0.0798));
    // data_cellSurvival_PC3_Flu.push_back(make_tuple(100.0,    0.729,    0.0720));
    // data_cellSurvival_PC3_Flu.push_back(make_tuple(150.0,    0.690,    0.0698));


    // CellSurvival cellSurvival_PC3_Flu = CellSurvival("PC3_Flu");
    // cellSurvival_PC3_Flu.AddCellSurvivalData(data_cellSurvival_PC3_Flu);
    // FitCellSurvival(cellSurvival_PC3_Flu, "LM", 1);
    // FitCellSurvival(cellSurvival_PC3_Flu, "LM", 2);
    // FitCellSurvival(cellSurvival_PC3_Flu, "LM", 3);

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

    // // std::vector<double,double> C4_2_Parameters;
    FitCellSurvival(cellSurvival_C4_2, "LM", 1);
    // FitCellSurvival(cellSurvival_C4_2, "LM", 2);
    // FitCellSurvival(cellSurvival_C4_2, "LM", 3);

};