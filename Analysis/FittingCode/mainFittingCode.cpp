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
        std::vector<std::tuple<double, TH1D*, TH1D*>> Get_hDose_Activity_Vec(){return hDose_Activity_Vec;};
        std::vector<std::tuple<double, TH1D*>> Get_hHitMultiplicity_vec(){return hHitMultiplicity_vec;};
        std::vector<std::tuple<double, TH2D*>> Get_hDose_HitsCellComponent_Vec(){return hDose_HitsCellComponent_Vec;};

    private:
        double n;
        std::vector<std::tuple<double, TH1D*, TH1D*>> hDose_Activity_Vec;

        std::vector<std::tuple<double, TH1D*>> hHitMultiplicity_vec;

        std::vector<std::tuple<double, TH2D*>> hDose_HitsCellComponent_Vec;
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
    std::string generalFilePath = "../OutputAnalysisCode/Output_" + cellSurvivalInstance.GetCellLine() + "_";

    // Specific filepath
    std::string filePath;

    // Make one for 0kBq activity here
    // Also 1kBq and 5kBq for the ones that are missing uptake at 5kBq

    //---------------------------
    // Making zero activity case by cloning the first histogram and resetting it

    // Making filename
    double firstActivity = std::get<0>(cellSurvivalData[0]);
    filePath = generalFilePath + std::to_string(((int)firstActivity)) + "kBq.root";

    // Read file containing histogram
    TFile* inputFile = new TFile(filePath.c_str(), "READ");

    //-------------------------------------
    // Extracting dose deposition histogram for cell component with mGY binning
    TH1D* hDose_CellComponent_mGy = nullptr;
    std::string histogramName_mGy = "i0_hDose_212Pb_" + cellSurvivalInstance.GetCellLine() + "_" + std::to_string(((int)firstActivity)) +"kBq_" + regionName + "_mGyBinning";
    std::string histogramNewName = "i0_hDose_212Pb_" + cellSurvivalInstance.GetCellLine() + "_" + std::to_string(0) +"kBq_" + regionName + "_mGyBinning";


    inputFile->GetObject(histogramName_mGy.c_str(), hDose_CellComponent_mGy);
    hDose_CellComponent_mGy->SetDirectory(0);
    hDose_CellComponent_mGy->SetName(histogramNewName.c_str());
    hDose_CellComponent_mGy->Reset();


    //-------------------------------------
    // Extracting dose deposition histogram for nucleus with uGyBinning
    TH1D* hDose_CellComponent_uGy = nullptr;
    std::string histogramName_uGy = "i0_hDose_212Pb_" + cellSurvivalInstance.GetCellLine() + "_" + std::to_string(((int)firstActivity)) +"kBq_" + regionName + "_uGyBinning";
    histogramNewName = "i0_hDose_212Pb_" + cellSurvivalInstance.GetCellLine() + "_" + std::to_string(0) +"kBq_" + regionName + "_uGyBinning";

    inputFile->GetObject(histogramName_uGy.c_str(), hDose_CellComponent_uGy);
    hDose_CellComponent_uGy->SetDirectory(0);
    hDose_CellComponent_uGy->SetName(histogramNewName.c_str());
    hDose_CellComponent_uGy->Reset();


    //-------------------------------------
    // Extracting 2D histogram for dose in cell component per number of alpha-particle hits
    TH2D* hDose_hitsAlpha_CellComponent = nullptr;
    std::string histogramName_Dose_HitsAlpha_CellComponent = "i0_hDose_212Pb_" + cellSurvivalInstance.GetCellLine() + "_" + std::to_string(((int)firstActivity)) +"kBq_HitsAlpha_" + regionName;
    histogramNewName = "i0_hDose_212Pb_" + cellSurvivalInstance.GetCellLine() + "_" + std::to_string(0) +"kBq_HitsAlpha_" + regionName;


    inputFile->GetObject(histogramName_Dose_HitsAlpha_CellComponent.c_str(), hDose_hitsAlpha_CellComponent);
    hDose_hitsAlpha_CellComponent->SetDirectory(0);
    hDose_hitsAlpha_CellComponent->SetName(histogramNewName.c_str());
    hDose_hitsAlpha_CellComponent->Reset();


    //-------------------------------------
    // Extracting hit multiplicity histogram for cell component
    TH1D* hHitMultiplicity_CellComponent = nullptr;
    std::string histogramName_hitMultiplicity = "i0_hFractionHitsAlpha_212Pb_" + cellSurvivalInstance.GetCellLine() + "_" + std::to_string(((int)firstActivity)) +"kBq_" + regionName;

    inputFile->GetObject(histogramName_hitMultiplicity.c_str(), hHitMultiplicity_CellComponent);
    hHitMultiplicity_CellComponent->SetDirectory(0);
    hHitMultiplicity_CellComponent->SetName(histogramNewName.c_str());
    hHitMultiplicity_CellComponent->Reset();


    //-------------------------------------
    inputFile->Close();


    hDose_Activity_Vec.push_back(std::make_tuple(0.,hDose_CellComponent_uGy,hDose_CellComponent_mGy));
    hHitMultiplicity_vec.push_back(std::make_tuple(0,hHitMultiplicity_CellComponent));
    hDose_HitsCellComponent_Vec.push_back(std::make_tuple(0.,hDose_hitsAlpha_CellComponent));

    //--------------------------------
    // Looping through all activities other activities
    for(int i=0; i<cellSurvivalData.size(); i++)
    {
        // Extracting activity
        double activity = std::get<0>(cellSurvivalData[i]);

        // Making filename
        filePath = generalFilePath + std::to_string(((int)activity)) + "kBq.root";

        // Read file containing histogram
        TFile* inputFile = new TFile(filePath.c_str(), "READ");


        //-------------------------------------
        // Extracting dose deposition histogram for nucleus with keV binning
        TH1D* hDose_CellComponent_mGy = nullptr;
        std::string histogramName_mGy = "i0_hDose_212Pb_" + cellSurvivalInstance.GetCellLine() + "_" + std::to_string(((int)activity)) +"kBq_" + regionName + "_mGyBinning";

        inputFile->GetObject(histogramName_mGy.c_str(), hDose_CellComponent_mGy);
        hDose_CellComponent_mGy->SetDirectory(0);


        //-------------------------------------
        // Extracting dose deposition histogram for nucleus with eV binning
        TH1D* hDose_CellComponent_uGy = nullptr;
        std::string histogramName_uGy = "i0_hDose_212Pb_" + cellSurvivalInstance.GetCellLine() + "_" + std::to_string(((int)activity)) +"kBq_" + regionName + "_uGyBinning";

        inputFile->GetObject(histogramName_uGy.c_str(), hDose_CellComponent_uGy);
        hDose_CellComponent_uGy->SetDirectory(0);

        //-------------------------------------
        // Extracting 2D histogram for dose in cell component per number of alpha-particle hits
        TH2D* hDose_hitsAlpha_CellComponent = nullptr;
        std::string histogramName_Dose_HitsAlpha_CellComponent = "i0_hDose_212Pb_" + cellSurvivalInstance.GetCellLine() + "_" + std::to_string(((int)activity)) +"kBq_HitsAlpha_" + regionName;

        inputFile->GetObject(histogramName_Dose_HitsAlpha_CellComponent.c_str(), hDose_hitsAlpha_CellComponent);
        hDose_hitsAlpha_CellComponent->SetDirectory(0);

        //-------------------------------------
        // Extracting hit multiplicity histogram for cell component
        TH1D* hHitMultiplicity_CellComponent = nullptr;
        std::string histogramName_hitMultiplicity = "i0_hFractionHitsAlpha_212Pb_" + cellSurvivalInstance.GetCellLine() + "_" + std::to_string(((int)activity)) +"kBq_" + regionName;

        inputFile->GetObject(histogramName_hitMultiplicity.c_str(), hHitMultiplicity_CellComponent);
        hHitMultiplicity_CellComponent->SetDirectory(0);


        //-------------------------------------
        inputFile->Close();


        hDose_Activity_Vec.push_back(std::make_tuple(activity,hDose_CellComponent_uGy,hDose_CellComponent_mGy));
        hHitMultiplicity_vec.push_back(std::make_tuple(activity,hHitMultiplicity_CellComponent));
        hDose_HitsCellComponent_Vec.push_back(std::make_tuple(activity,hDose_hitsAlpha_CellComponent));
    }
}


//----------------------------------
void FitCellSurvival(CellSurvival cellSurvivalInstance, std::string modelName, int region)
{
    //--------------------------------
    // Function fitting the cell survival data provided with a model, either:
    // LM : Linear model, beta set to zero automatically
    // LQ : Linear-quadratic model
    // The region is
    // Returns a vector of tuples containing the fit parameters and their uncertainties

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


    std::vector<std::tuple<double,TH1D*,TH1D*>> vec_activities_histograms = energyDepHistograms.Get_hDose_Activity_Vec();

    // std::vector<std::tuple<double, TH1D*>> hHitMultiplicity_vec = energyDepHistograms.Get_hHitMultiplicity_vec();

    // std::vector<std::tuple<double, TH2D*>> hDose_HitsCellComponent_Vec = energyDepHistograms.Get_hDose_HitsCellComponent_Vec();

    // Factor histograms have been scaled by
    double scalingFactorHistogram = ((double)get<1>(vec_activities_histograms[1])->GetEntries())/((double) get<1>(vec_activities_histograms[1])->Integral());

    //------------------------------
    auto CalculateCellSurvivalFraction = [&](TH1D *hDose_uGy, TH1D *hDose_mGy, double *par)
    {
        double alpha = par[0];
        double beta = par[1];


        //------------------------------------------------
        double fractionOfComponentsHit = hDose_mGy->Integral();
        double fractionOfComponentsMissed = (1.0 - fractionOfComponentsHit);

        //------------------------------------------------
        //      The fraction of missed cells obviously all survive so that is immediately added.
        double fractionOfTotalCellsHitSurviving = fractionOfComponentsMissed;

        // Looping over uGy binned histogram
        for(int i=0; i<hDose_uGy->GetNbinsX();i++)
        {
            double dose = hDose_uGy->GetBinCenter(i+1);
            double fractionOfTotalCellsHit = hDose_uGy->GetBinContent(i+1);

            //----------------------------
            double cellSurvivalfraction = TMath::Exp(-(alpha*dose + beta*TMath::Power(dose, 2.0)));
            double fractionOfTotalCellsHit_survivingFraction = cellSurvivalfraction*fractionOfTotalCellsHit;


            fractionOfTotalCellsHitSurviving += fractionOfTotalCellsHit_survivingFraction;
        }

        // Looping over mGy binned histogram
        for(int i=1000; i<hDose_mGy->GetNbinsX(); i++)
        {
            double dose = hDose_mGy->GetBinCenter(i+1);
            double fractionOfTotalCellsHit = hDose_mGy->GetBinContent(i+1);

            //----------------------------
            double cellSurvivalfraction = TMath::Exp(-(alpha*dose + beta*TMath::Power(dose, 2.0)));
            double fractionOfTotalCellsHit_survivingFraction = cellSurvivalfraction*fractionOfTotalCellsHit;

            fractionOfTotalCellsHitSurviving += fractionOfTotalCellsHit_survivingFraction;
        }

        //------------------------------------------------
        return fractionOfTotalCellsHitSurviving;
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
            auto histDose_uGy = std::get<1>(entry);
            auto histDose_mGy = std::get<2>(entry);

            double cellSurvival = CalculateCellSurvivalFraction(histDose_uGy, histDose_mGy, par);

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
    auto f_cellSurvivalVsDose = new TF1("f_cellSurvivalVsDose",
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

    f_cellSurvivalVsDose->SetParLimits(0, 0.0, 1.0e+05);
    f_cellSurvivalVsDose->SetParLimits(1, 0.0, 1.0e+05);

    f_cellSurvivalVsDose->SetNpx(10000);
    f_cellSurvivalVsDose->SetParameter(0, 1.0e+00);

    if(modelName=="LM")
    {
        f_cellSurvivalVsDose->FixParameter(1, 0.0e+00);
    }
    if(modelName=="LQ")
    {
        f_cellSurvivalVsDose->SetParameter(1, 1.0e+00);
    }

    std::cout << " Cell Line " + cellSurvivalInstance.GetCellLine() << " Model : " << modelName << " Region : " << regionName << std::endl;
    gr_clonogenicSurvival->Fit(f_cellSurvivalVsDose, "", "", 0.0, 150.0);

    double chi_sq = f_cellSurvivalVsDose->GetChisquare();

    int deg_freedom = f_cellSurvivalVsDose->GetNDF();

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
    f_cellSurvivalVsDose->Write();
    gr_cellSurvivability_vs_activitykBqPerMl.Write();

    auto legend = new TLegend(0.1,0.7,0.48,0.9);
    legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
    legend->AddEntry(gr_clonogenicSurvival,"Experimental Data","f");
    legend->AddEntry(f_cellSurvivalVsDose,"Simulated Data","f");

    legend->Write();


    for(auto & entry : vec_activities_histograms)
    {
        auto hist = std::get<2>(entry);

        hist->Write();
    }

    //-----------------------
    // Extracting fit parameters

    std::vector<std::tuple<double,double>> parametersAndUncertainties_Vec;

    double alpha = f_cellSurvivalVsDose->GetParameter(0);
    double dAlpha = f_cellSurvivalVsDose->GetParError(0);

    parametersAndUncertainties_Vec.push_back(std::make_tuple(alpha,dAlpha));

    if(modelName=="LQ")
    {
        double beta = f_cellSurvivalVsDose->GetParameter(1);
        double dBeta = f_cellSurvivalVsDose->GetParError(1);

        parametersAndUncertainties_Vec.push_back(std::make_tuple(beta,dBeta));
    }



    //------------------------
    // Hit Analysis



    //--------------------------------------------
    auto Fill_hHitMultiplicity_SurvivalFraction_OneActivity = [&](double activity, TH2D* hDose_hitsAlpha_CellComponent_OneActivity, TH1D* hHitMultiplicity_SurvivalFraction_OneActivity, TH1D* hDose_mGy)
    {

        //---------------------------------
        // Function to calculate the error in the survival fraction of one bin at a certain dose
        auto CalculateUncertainty_HitMultiplicity_SurvivalFraction_OneBin = [&](double doseInBin, double fractionSurvivedDoseBin)
        {
            double dFractionSurvivedDoseBin = doseInBin*fractionSurvivedDoseBin*dAlpha;
            return dFractionSurvivedDoseBin;
        };

        //---------------------------
        // Vector to store fraction of cells survived at each number of hits
        // <nHits, survivalFractionNHits, uncertainty survivalFractionNHits>
        std::vector<std::tuple<int,double,double>> hitMultiplicity_SurvivalFraction_Vec;

        double fractionCellsNotHit;

        if(activity<=0.)
        {
            fractionCellsNotHit = 1. - hDose_mGy->Integral();
        }
        else
        {
            fractionCellsNotHit = 0.;
        }

        //----------------------------
        // Loop over number of hits (y-axis)
        for(int i=0; i<hDose_hitsAlpha_CellComponent_OneActivity->GetNbinsY(); i++)
        {
            // Those not hit always survive
            double survivalFractionThisHitNumber = fractionCellsNotHit;

            // Vector to store the uncertainty in the survival fraction at every bin
            std::vector<double> dSurvivaFractionThisDose_Vec;

            //--------------------------
            // Loop over x axis (dose deposited)
            for(int j=0; j<hDose_hitsAlpha_CellComponent_OneActivity->GetNbinsX(); j++)
            {
                double doseDep = (hDose_hitsAlpha_CellComponent_OneActivity->GetXaxis())->GetBinCenter(j+1);
                double fractionHitThisDose = hDose_hitsAlpha_CellComponent_OneActivity->GetBinContent(j+1,i+1);

                // Only calculate survival if some fraction has dose deposited
                if(fractionHitThisDose>0.)
                {
                    double survivaFractionThisDose = fractionHitThisDose*TMath::Exp(-alpha*doseDep);
                    double dSurvivaFractionThisDose = CalculateUncertainty_HitMultiplicity_SurvivalFraction_OneBin(doseDep,survivaFractionThisDose);
                    dSurvivaFractionThisDose_Vec.push_back(dSurvivaFractionThisDose);

                    // Adding survival fraction to total survival
                    survivalFractionThisHitNumber += survivaFractionThisDose;
                }
            }

            //------------------------------
            // Calculating uncertainty in survival fraction for hit number
            double dSurvivalFractionThisHitNumber_squared = 0.;
            for(auto& entry : dSurvivaFractionThisDose_Vec)
            {
                dSurvivalFractionThisHitNumber_squared += std::pow(entry,2.);
            }
            double dSurvivalFractionThisHitNumber = std::sqrt(dSurvivalFractionThisHitNumber_squared);
            hitMultiplicity_SurvivalFraction_Vec.push_back(std::make_tuple(i,survivalFractionThisHitNumber, dSurvivalFractionThisHitNumber));
        }

        //-------------------------
        // Adding survival fraction to histogram
        for(auto& entry : hitMultiplicity_SurvivalFraction_Vec)
        {
            int nHits = get<0>(entry);
            double survivalFraction = get<1>(entry);
            double dSurvivalFraction = get<2>(entry);

            hHitMultiplicity_SurvivalFraction_OneActivity->SetBinContent(nHits+1,survivalFraction);
            hHitMultiplicity_SurvivalFraction_OneActivity->SetBinError(nHits+1, dSurvivalFraction);
        }

        return hHitMultiplicity_SurvivalFraction_OneActivity;
    };



    //-----------------------------
    auto Fill_hHitMultiplicity_PercentKilled_OneActivity = [&](TH1D* hHitMultiplicity_SurvivalFraction_OneActivity, TH1D* hHitMultiplicity_OneActivity, TH1D* hHitMultiplicity_PercentKilled_OneActivity)
    {
        //----------------------------
        auto CalculateUncertainty_PercentKilled_NHits = [&](double fractionHitNHits, double dFractionSurvivedNHits)
        {
            double dPercentDeathNHits = (100./fractionHitNHits)*dFractionSurvivedNHits;

            return dPercentDeathNHits;
        };

        //-------------------------
        // Looping over number of hits to cell component
        for(int i=0; i<hHitMultiplicity_SurvivalFraction_OneActivity->GetNbinsX(); i++)
        {
            double fractionHitThisNHits = hHitMultiplicity_OneActivity->GetBinContent(i+1);

            if(fractionHitThisNHits>0)
            {
                double fractionSurvivedThisNHits = hHitMultiplicity_SurvivalFraction_OneActivity->GetBinContent(i+1);
                double dFractionSurvivedThisNHits = hHitMultiplicity_SurvivalFraction_OneActivity->GetBinError(i+1);

                double percentKilled = 100.*(fractionHitThisNHits-fractionSurvivedThisNHits)/fractionHitThisNHits;
                double dPercentKilled = CalculateUncertainty_PercentKilled_NHits(fractionHitThisNHits,dFractionSurvivedThisNHits);

                hHitMultiplicity_PercentKilled_OneActivity->SetBinContent(i+1, percentKilled);
                hHitMultiplicity_PercentKilled_OneActivity->SetBinError(i+1, dPercentKilled);
            }
        }

        return hHitMultiplicity_PercentKilled_OneActivity;
    };


    //----------------------------
    // Function returning average dose per NHits, for activites
    auto Fill_grMeanDose_PerHits_OneActivity = [&](double activity, TH2D* hDose_hitsAlpha_CellComponent_OneActivity, TGraphErrors* grMeanDose_PerOneHit)
    {
        double intNHits = 0.;
        int graphPointN = 0;

        //---------------------------
        // Looping over number of hits (y-axis)
        for(int i=0; i<hDose_hitsAlpha_CellComponent_OneActivity->GetNbinsY();i++)
        {
            //---------------------------
            // Extract dose histogram for i number of hits
            std::string nameProjectedHist = "doseDelivered_NHits_" + std::to_string(i) + "_" + std::to_string((int)activity) + "kBq";
            TH1D* projectedDose_ForNHits = hDose_hitsAlpha_CellComponent_OneActivity->ProjectionX(nameProjectedHist.c_str(),i+1,i+1);

            double intDose_ThisN = projectedDose_ForNHits->Integral();

            if(intDose_ThisN>0)
            {
                // Normalizing dose distribution
                projectedDose_ForNHits->Scale(1./intDose_ThisN);

                // expectation value of dose, and dose squared
                double ei_di = 0.;
                double ei_di_squared = 0.;

                // Looping over all doses
                for(int j=0; j<projectedDose_ForNHits->GetNbinsX();j++)
                {
                    double di = projectedDose_ForNHits->GetBinCenter(j+1);
                    double pi = projectedDose_ForNHits->GetBinContent(j+1);
                    ei_di += di*pi;
                    ei_di_squared += std::pow(di,2.)*pi;
                }

                double ei_di_variance = ei_di_squared - std::pow(ei_di,2.);
                double ei_di_std;
                if(ei_di_variance<0.){ei_di_std = std::sqrt(ei_di_variance);}
                else{ei_di_std = 0.05;}


                grMeanDose_PerOneHit->SetPoint(graphPointN,((double)i),ei_di);
                grMeanDose_PerOneHit->SetPointError(graphPointN,0.0,ei_di_std);
                graphPointN++;

            }
            projectedDose_ForNHits->SetDirectory(0);
        }
        return grMeanDose_PerOneHit;
    };


    //----------------------------
    // Function returning average dose per NHits, for activites
    auto Fill_grMeanDose_PerHits_OneActivity_ASym = [&](double activity, TH2D* hDose_hitsAlpha_CellComponent_OneActivity, TGraphAsymmErrors* grMeanDose_PerOneHit_ASym)
    {
        int graphPointN = 0;

        //---------------------------
        // Looping over number of hits (y-axis)
        for(int i=0; i<hDose_hitsAlpha_CellComponent_OneActivity->GetNbinsY();i++)
        {
            //---------------------------
            // Extract dose histogram for i number of hits
            std::string nameProjectedHist = "doseDelivered_NHits_" + std::to_string(i) + "_" + std::to_string((int)activity) + "kBq";
            TH1D* projectedDose_ForNHits = hDose_hitsAlpha_CellComponent_OneActivity->ProjectionX(nameProjectedHist.c_str(),i+1,i+1);


            double intDose_ThisN = projectedDose_ForNHits->Integral();

            if(intDose_ThisN>0.)
            {
                // double meanDose = projectedDose_ForNHits->GetMean(1);
                // int binAtMean = projectedDose_ForNHits->FindBin(meanDose);
                // double confInterval = 0.5*intDose_ThisN;

                // //-----------------------
                // // Finding quartiles of dose distribution
                // double doseAtIntervalPlus;
                // double doseAtIntervalMinus;
                // double binW = projectedDose_ForNHits->GetBinWidth(1);

                // double integralPlusSide = 0.;
                // int k = 0;
                // while(integralPlusSide<(confInterval/2.))
                // {
                //     if(k==0){integralPlusSide = 0.5*projectedDose_ForNHits->Integral(binAtMean,binAtMean+k);}
                //     else{integralPlusSide = projectedDose_ForNHits->Integral(binAtMean,binAtMean+k);}
                //     doseAtIntervalPlus = projectedDose_ForNHits->GetBinCenter(binAtMean+k);
                //     k++;
                // }

                // double integralMinusSide = 0.;
                // k = 0;
                // while(integralMinusSide<(confInterval/2.))
                // {
                //     if(k==0){integralMinusSide = 0.5*projectedDose_ForNHits->Integral(binAtMean,binAtMean-k);}
                //     else{integralMinusSide = projectedDose_ForNHits->Integral(binAtMean,binAtMean-k);}
                //     doseAtIntervalMinus = meanDose - k*binW;
                //     k++;
                // }

                // double uncertaintyHigh = doseAtIntervalPlus - meanDose;
                // double uncertaintyLow = meanDose - doseAtIntervalMinus;

                // if(uncertaintyHigh<=0.)
                // {
                //     uncertaintyHigh = 0.02*meanDose;
                //     uncertaintyLow = 0.01*meanDose;
                // }

                // // std::cout << "N : " << i << " M : " << meanDose << "+" << uncertaintyHigh << "-" << uncertaintyLow << std::endl;

                int n = 3;
                double x_q[3];
                double y_q[3] = {0.25, 0.50, 0.75}; // quartile positions


                projectedDose_ForNHits->GetQuantiles(n, x_q, y_q);
                double doseAtQuartileMinus = x_q[0];
                double doseAtMean = x_q[1];
                double doseAtQuartilePlus = x_q[2];
                //-------------------
                // Setting points
                grMeanDose_PerOneHit_ASym->SetPoint(graphPointN, ((double)i), doseAtMean);
                grMeanDose_PerOneHit_ASym->SetPointEYhigh(graphPointN, doseAtQuartilePlus - doseAtMean);
                grMeanDose_PerOneHit_ASym->SetPointEYlow(graphPointN, doseAtMean - doseAtQuartileMinus);
                graphPointN++;
            }
            projectedDose_ForNHits->SetDirectory(0);
        }
        return grMeanDose_PerOneHit_ASym;
    };


    //---------------------------
    // Filling histograms


    std::vector<std::tuple<double, TH1D*>> hHitMultiplicity_vec = energyDepHistograms.Get_hHitMultiplicity_vec();

    std::vector<std::tuple<double, TH2D*>> hDose_HitsCellComponent_Vec = energyDepHistograms.Get_hDose_HitsCellComponent_Vec();

    std::vector<std::tuple<double,TH1D*>> hHitMultiplicity_PercentKilled_Vec;

    std::vector<std::tuple<double,TGraphErrors*>> grMeanDose_PerHits_Vec;

    //-----------------
    // Looping over every activity
    for(int i = 0; i<hHitMultiplicity_vec.size(); i++)
    {
        double activity = std::get<0>(hHitMultiplicity_vec[i]);
        TH1D* hHitMultiplicity = std::get<1>(hHitMultiplicity_vec[i]);
        TH2D* hDose_Hits = std::get<1>(hDose_HitsCellComponent_Vec[i]);
        TH1D* hDose_mGy = std::get<2>(vec_activities_histograms[i]);


        //-----------------------------
        // Define number of bins and max number of hits
        int nBins = hDose_Hits->GetNbinsY();
        // double hitsMax = 700.;

        double hitsMax = hHitMultiplicity->GetXaxis()->GetXmax();

        //--------------------------
        // Making histogram for cell survival fraction for hit multiplicity
        std::string hHitMultiplicity_SurvivalFraction_ThisActivity_Name = "hHitMultiplicity_" + regionName + "_SurvivalFraction_" + std::to_string((int)activity) + "kBq";
        TH1D* hHitMultiplicity_SurvivalFraction_ThisActivity = new TH1D(hHitMultiplicity_SurvivalFraction_ThisActivity_Name.c_str(), "Cell Survival Fraction per N Number of Hits by Alpha Particle", nBins, 0.,hitsMax);
        hHitMultiplicity_SurvivalFraction_ThisActivity->GetXaxis()->SetTitle("N Number of Hits to Cell Component by Alpha Particle");
        hHitMultiplicity_SurvivalFraction_ThisActivity->GetYaxis()->SetTitle("Fraction of Total Cells Survived");

        //--------------------------
        // Making histogram for percent of cells killed for hit multiplicity
        std::string hHitMultiplicity_PercentKilled_Name = "hHitMultiplicity_" + regionName + "_PercentKilled_" + std::to_string((int)activity) + "kBq";
        TH1D* hHitMultiplicity_PercentKilled = new TH1D(hHitMultiplicity_PercentKilled_Name.c_str(), "Cell Death per N Number of Hits by Alpha Particle", nBins, 0.,hitsMax);
        hHitMultiplicity_PercentKilled->GetXaxis()->SetTitle("N Number of Hits to Cell Component by Alpha Particle");
        hHitMultiplicity_PercentKilled->GetYaxis()->SetTitle("Percentage of Cells Killed");

        //----------------------------
        // Making graph for mean dose per number of hits for one activity
        std::string grMeanDose_PerHits_OneActivity_Name = "grMeanDose_PerNHits_" + std::to_string((int)activity) +"kBq";
        TGraphErrors* grMeanDose_PerHits_OneActivity = new TGraphErrors();
        grMeanDose_PerHits_OneActivity->SetTitle("Mean Dose Delivered to Cell Component Per Number of Alpha-Particle Hits");
        grMeanDose_PerHits_OneActivity->SetName(grMeanDose_PerHits_OneActivity_Name.c_str());
        grMeanDose_PerHits_OneActivity->GetXaxis()->SetTitle("N Number of Hits to Cell Component by Alpha Particle");
        grMeanDose_PerHits_OneActivity->GetYaxis()->SetTitle("Mean Dose Delivered to Cell Component");

        //----------------------------
        // Making graph for mean dose per number of hits for one activity, asymetrical errors
        std::string grMeanDose_PerHits_OneActivity_ASym_Name = "grMeanDose_PerNHits_ASym_" + std::to_string((int)activity) +"kBq";
        TGraphAsymmErrors* grMeanDose_PerHits_OneActivity_ASym = new TGraphAsymmErrors();
        grMeanDose_PerHits_OneActivity_ASym->SetTitle("Mean Dose Delivered to Cell Component Per Number of Alpha-Particle Hits");
        grMeanDose_PerHits_OneActivity_ASym->SetName(grMeanDose_PerHits_OneActivity_ASym_Name.c_str());
        grMeanDose_PerHits_OneActivity_ASym->GetXaxis()->SetTitle("N Number of Hits to Cell Component by Alpha Particle");
        grMeanDose_PerHits_OneActivity_ASym->GetYaxis()->SetTitle("Mean Dose Delivered to Cell Component");


        //------------------------------------
        hHitMultiplicity_SurvivalFraction_ThisActivity = Fill_hHitMultiplicity_SurvivalFraction_OneActivity(activity, hDose_Hits, hHitMultiplicity_SurvivalFraction_ThisActivity, hDose_mGy);

        //-----------------------------------
        hHitMultiplicity_PercentKilled = Fill_hHitMultiplicity_PercentKilled_OneActivity(hHitMultiplicity_SurvivalFraction_ThisActivity, hHitMultiplicity, hHitMultiplicity_PercentKilled);
        hHitMultiplicity_PercentKilled_Vec.push_back(std::make_tuple(activity,hHitMultiplicity_PercentKilled));

        //-------------------------------
        grMeanDose_PerHits_OneActivity = Fill_grMeanDose_PerHits_OneActivity(activity, hDose_Hits, grMeanDose_PerHits_OneActivity);
        grMeanDose_PerHits_Vec.push_back(std::make_tuple(activity,grMeanDose_PerHits_OneActivity));
        grMeanDose_PerHits_OneActivity->Write();


        //----------------------------
        grMeanDose_PerHits_OneActivity_ASym = Fill_grMeanDose_PerHits_OneActivity_ASym(activity, hDose_Hits, grMeanDose_PerHits_OneActivity_ASym);
        grMeanDose_PerHits_OneActivity_ASym->Write();



        //-----------------------------
        // 10 kBq case, saving the dose deposition histograms for different NHit
        if(activity==25.)
        {
            for(int j=0; j<21; j++)
            {
                std::string nameProjectedHist = "hDoseDelivered_NHits_" + std::to_string(j) + "_" + std::to_string((int)activity) + "kBq";
                TH1D* projectedDose_ForNHits = hDose_Hits->ProjectionX(nameProjectedHist.c_str(),j+1,j+1);
                projectedDose_ForNHits->SetFillColorAlpha(kGreen+2, 0.3);
                if(projectedDose_ForNHits->GetNbinsX()>0)
                {
                    projectedDose_ForNHits->Write();
                }
                projectedDose_ForNHits->SetDirectory(0);
            }
        }
    }


    //---------------------------
    auto Fill_grPercentKilled_PerOneNHit = [&](int NHits, TGraphErrors* grPercentKilled_PerOneNHit)
    {
        std::vector<std::tuple<double,double,double>> percentKilled_Activity_ThisNHit_Vec;
        for(int i=0; i<hHitMultiplicity_PercentKilled_Vec.size(); i++)
        {
            double activity = std::get<0>(hHitMultiplicity_PercentKilled_Vec[i]);
            TH1D* hHitMultiplicity_PercentKilled_OneActivity = std::get<1>(hHitMultiplicity_PercentKilled_Vec[i]);

            if(hHitMultiplicity_PercentKilled_OneActivity->GetNbinsX()>0)
            {
                double percentKilled_NHits = hHitMultiplicity_PercentKilled_OneActivity->GetBinContent(NHits+1);
                double dPercentKilled_NHits = hHitMultiplicity_PercentKilled_OneActivity->GetBinError(NHits+1);
                // std::cout << percentKilled_NHits << std::endl;

                if(percentKilled_NHits>0.){percentKilled_Activity_ThisNHit_Vec.push_back(std::make_tuple(activity,percentKilled_NHits,dPercentKilled_NHits));}
            }
        }

        for(int i=0; i<percentKilled_Activity_ThisNHit_Vec.size(); i++)
        {
            double activity = std::get<0>(percentKilled_Activity_ThisNHit_Vec[i]);
            double percentKilled_NHits = std::get<1>(percentKilled_Activity_ThisNHit_Vec[i]);
            double dPercentKilled_NHits = std::get<2>(percentKilled_Activity_ThisNHit_Vec[i]);

            grPercentKilled_PerOneNHit->SetPoint(i, activity, percentKilled_NHits);
            grPercentKilled_PerOneNHit->SetPointError(i, 0., dPercentKilled_NHits);
        }
        return percentKilled_Activity_ThisNHit_Vec;
    };


    //---------------------------------
    auto Fill_grMeanDose_PerOneHit = [&](int NHits, TGraphErrors* grMeanDose_PerOneHit)
    {
        std::vector<std::tuple<double,double,double>> meanDose_Activity_ThisHit_Vec;
        for(int i=0; i<grMeanDose_PerHits_Vec.size();i++)
        {
            double activity = std::get<0>(grMeanDose_PerHits_Vec[i]);
            TGraphErrors* grMeanDose_OneActivity = std::get<1>(grMeanDose_PerHits_Vec[i]);

            if(grMeanDose_OneActivity->GetN()>0)
            {
                double *hits = grMeanDose_OneActivity->GetX();
                double *meanDose = grMeanDose_OneActivity->GetY();
                double meanDoseAtNHits;
                double dMeanDoseAtNHits;
                for(int j=0; j<grMeanDose_OneActivity->GetN();j++)
                {
                    if(hits[j]==((double)NHits))
                    {
                        meanDoseAtNHits = meanDose[j];
                        dMeanDoseAtNHits = grMeanDose_OneActivity->GetErrorY(j);
                        break;
                    }
                }
                double tolerance = std::pow(10.,-14.);
                if(meanDoseAtNHits>tolerance){meanDose_Activity_ThisHit_Vec.push_back(std::make_tuple(activity,meanDoseAtNHits,dMeanDoseAtNHits));}
            }

        }

        for(int i=0; i<meanDose_Activity_ThisHit_Vec.size();i++)
        {
            double activity = std::get<0>(meanDose_Activity_ThisHit_Vec[i]);
            double meanDose = std::get<1>(meanDose_Activity_ThisHit_Vec[i]);
            double dMeanDose = std::get<2>(meanDose_Activity_ThisHit_Vec[i]);

            grMeanDose_PerOneHit->SetPoint(i, activity, meanDose);
            grMeanDose_PerOneHit->SetPointError(i, 0.0, dMeanDose);
        }

        return meanDose_Activity_ThisHit_Vec;
    };


    //--------------------------------
    // Making graph for average probability of death at N number of hits
    TGraphErrors* grUWA_percentDeath_perN = new TGraphErrors();
    grUWA_percentDeath_perN->SetName("grUWA_percentDeath");
    grUWA_percentDeath_perN->GetXaxis()->SetTitle("Number of Hits by Alpha Particle");
    grUWA_percentDeath_perN->GetYaxis()->SetTitle("Average Probability of Death at N Number of Hits");

    //--------------------------------
    // Making graph for average mean dose at N Number of hits
    TGraphErrors* grUWA_meanDose_perN = new TGraphErrors();
    grUWA_meanDose_perN->SetName("grUWA_meanDose_perN");
    grUWA_meanDose_perN->GetXaxis()->SetTitle("Number of Hits by Alpha Particle");
    grUWA_meanDose_perN->GetYaxis()->SetTitle("Average Mean Dose at N Number of Hits");


    //--------------------------

    double doseMax_histograms = (std::get<1>(hDose_HitsCellComponent_Vec[1]))->GetXaxis()->GetXmax();
    double nBins_dose = (std::get<1>(hDose_HitsCellComponent_Vec[1]))->GetNbinsX();


    std::vector<std::tuple<int,TH1D*>> overlappedDoses_PerNHit_Vec;

    for(int i=0; i<41; i++)
    {

        //--------------------------------
        // Making graph for mean dose for activity, for NHits
        std::string grMeanDose_i_Name = "grMeanDose_" + std::to_string(i) + "_NumberHitsToCellComponent";
        TGraphErrors* grMeanDose_i = new TGraphErrors();
        grMeanDose_i->SetName(grMeanDose_i_Name.c_str());
        grMeanDose_i->GetXaxis()->SetTitle("Activity [kBq/mL]");
        grMeanDose_i->GetYaxis()->SetTitle("Mean Dose Delivered to Cell Component");

        std::vector<std::tuple<double,double,double>> meanDose_Activity_Vec = Fill_grMeanDose_PerOneHit(i,grMeanDose_i);

        if(meanDose_Activity_Vec.size()>0)
        {
            grMeanDose_i->Write();

            //----------------------------------
            // Calculating mean dose averaged over all activities per number of hits
            double sum_w_i = 0.;
            double sum_mu_i_w = 0.;

            for(auto& entry : meanDose_Activity_Vec)
            {
                double meanDoseOneActivity = std::get<1>(entry);
                double dMeanDoseOneActivity = std::get<2>(entry);

                double w_i = 1./std::pow(dMeanDoseOneActivity,2.);
                double mu_i_w = w_i*meanDoseOneActivity;

                sum_w_i += w_i;
                sum_mu_i_w += mu_i_w;
            }

            double UWA_meanDose_perN = sum_mu_i_w/sum_w_i;
            double dUWA_meanDose_perN = 1./std::sqrt(sum_w_i);

            grUWA_meanDose_perN->SetPoint(i, ((double)i), UWA_meanDose_perN);
            grUWA_meanDose_perN->SetPointError(i, 0.0, dUWA_meanDose_perN);
        }


        //----------------------------
        std::string gr_percentDeath_i_Name = "grProbabilityDeath_" + std::to_string(i) + "_NumberHitsToCellComponent";
        TGraphErrors* gr_percentDeath_i = new TGraphErrors();
        gr_percentDeath_i->SetName(gr_percentDeath_i_Name.c_str());
        gr_percentDeath_i->GetXaxis()->SetTitle("Activity [kBq/mL]");
        gr_percentDeath_i->GetYaxis()->SetTitle("Probability of Cell Death");

        // Filling percent death for activity per number of hits
        std::vector<std::tuple<double,double,double>> percentKilled_Activity_Vec = Fill_grPercentKilled_PerOneNHit(i,gr_percentDeath_i);
        if(percentKilled_Activity_Vec.size()>0)
        {
            gr_percentDeath_i->Write();

            //----------------------------------
            // Calculating percent death averaged over all activities per number of hits
            double sum_w_i = 0.;
            double sum_mu_i_w = 0.;

            for(auto & entry : percentKilled_Activity_Vec)
            {
                double percentDeathOneActivity = std::get<1>(entry);
                double dPercentDeathOneActivity = std::get<2>(entry);
                // std::cout << percentDeathOneActivity << " " << dPercentDeathOneActivity << std::endl;

                double w_i = 1./std::pow(dPercentDeathOneActivity,2.);
                double mu_i_w = w_i*percentDeathOneActivity;

                sum_w_i += w_i;
                sum_mu_i_w += mu_i_w;
            }

            double UWA_percentDeath_perN = sum_mu_i_w/sum_w_i;
            // std::cout << sum_w_i << " " << sum_mu_i_w <<  std::endl;
            double dUWA_percentDeath_perN = 1./std::sqrt(sum_w_i);

            grUWA_percentDeath_perN->SetPoint(i, ((double)i), UWA_percentDeath_perN);
            grUWA_percentDeath_perN->SetPointError(i, 0., dUWA_percentDeath_perN);
        }

        std::string hDose_ThisNHits_OverlappedAllActivities_Name = "hDose_NHits_" + std::to_string(i) + "_Overlapped_AllActivities";
        TH1D* hDose_ThisNHits_OverlappedAllActivities = new TH1D(hDose_ThisNHits_OverlappedAllActivities_Name.c_str(),"Dose Delivered for one NHit, All Activities Overlapped", nBins_dose, 0., doseMax_histograms);
        hDose_ThisNHits_OverlappedAllActivities->GetXaxis()->SetTitle("Dose Delivered [Gy]");
        hDose_ThisNHits_OverlappedAllActivities->GetXaxis()->SetTitle("Fraction of Cells Hit with Dose");


        int nHistogramsAdded = 0;
        //------------------------------------
        for(int j=0; j<hDose_HitsCellComponent_Vec.size(); j++)
        {
            double activity = std::get<0>(hDose_HitsCellComponent_Vec[j]);
            TH2D* hDose_Hits = std::get<1>(hDose_HitsCellComponent_Vec[j]);

            std::string hDose_ThisNHits_Name = "hDose_ThisNHits";
            TH1D* hDose_ThisNHits = hDose_Hits->ProjectionX(hDose_ThisNHits_Name.c_str(),i+1,i+1);
            if(hDose_ThisNHits->Integral()>0.)
            {
                hDose_ThisNHits_OverlappedAllActivities->Add(hDose_ThisNHits);
                nHistogramsAdded++;
            }
            hDose_ThisNHits->SetDirectory(0);
        }

        hDose_ThisNHits_OverlappedAllActivities->Scale(1./((double)nHistogramsAdded));

        overlappedDoses_PerNHit_Vec.push_back(std::make_tuple(i,hDose_ThisNHits_OverlappedAllActivities));
        hDose_ThisNHits_OverlappedAllActivities->Write();
        hDose_ThisNHits_OverlappedAllActivities->SetDirectory(0);

    }

    grUWA_percentDeath_perN->SetTitle("Uncertainty Weighted Average of Probability of Cell Death For N Number of Alpha-Particle Hits");
    grUWA_percentDeath_perN->Write();
    grUWA_meanDose_perN->SetTitle("Uncertainty Weighted Average of Mean Dose For N Number of Alpha-Particle Hits");
    grUWA_meanDose_perN->Write();


    //---------------------------------
    std::string gr_MeanDose_PerNHit_OverlappedActivities_Name = "grUWA_meanDose_OverlappedActivities";
    TGraphAsymmErrors* gr_MeanDose_PerNHit_OverlappedActivities = new TGraphAsymmErrors();
    gr_MeanDose_PerNHit_OverlappedActivities->SetName(gr_MeanDose_PerNHit_OverlappedActivities_Name.c_str());
    gr_MeanDose_PerNHit_OverlappedActivities->SetTitle("Mean Dose Per Alpha-Particle Hit to Cell Component, Averaged Over All Activities");
    gr_MeanDose_PerNHit_OverlappedActivities->GetXaxis()->SetTitle("N Number of Alpha-Particle Hits");
    gr_MeanDose_PerNHit_OverlappedActivities->GetYaxis()->SetTitle("Mean Dose Delivered [Gy]");

    int NGraphPoint = 0;
    for(int i=0; i<overlappedDoses_PerNHit_Vec.size(); i++)
    {
        int NHits = std::get<0>(overlappedDoses_PerNHit_Vec[i]);
        TH1D* hOverlappedDose = std::get<1>(overlappedDoses_PerNHit_Vec[i]);

        int n = 3;
        double x_q[3];
        double y_q[3] = {0.25, 0.50, 0.75}; // quartile positions

        hOverlappedDose->GetQuantiles(n, x_q, y_q);
        double doseAtQuartileMinus = x_q[0];
        double doseAtMean = x_q[1];
        double doseAtQuartilePlus = x_q[2];

        if(doseAtMean>0)
        {
            gr_MeanDose_PerNHit_OverlappedActivities->SetPoint(NGraphPoint, ((double)NHits), doseAtMean);
            gr_MeanDose_PerNHit_OverlappedActivities->SetPointEYhigh(NGraphPoint, doseAtQuartilePlus - doseAtMean);
            gr_MeanDose_PerNHit_OverlappedActivities->SetPointEYlow(NGraphPoint, doseAtMean - doseAtQuartileMinus);
            NGraphPoint++;
        }

    }

    gr_MeanDose_PerNHit_OverlappedActivities->Write();


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
    // // FitCellSurvival(cellSurvival_PC3_PIP, "LM", 1);
    // // FitCellSurvival(cellSurvival_PC3_PIP, "LM", 2);
    // FitCellSurvival(cellSurvival_PC3_PIP, "LM", 3);
    // // FitCellSurvival(cellSurvival_PC3_PIP, "LM", 4);


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
    // // FitCellSurvival(cellSurvival_PC3_Flu, "LM", 1);
    // // FitCellSurvival(cellSurvival_PC3_Flu, "LM", 2);
    // FitCellSurvival(cellSurvival_PC3_Flu, "LM", 3);
    // // FitCellSurvival(cellSurvival_PC3_Flu, "LM", 4);

    std::vector<std::tuple<double,double,double>> data_cellSurvival_C4_2;
    // data_cellSurvival_C4_2.push_back(make_tuple(0.0,    1.0,    0.05));
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

    // std::vector<double,double> C4_2_Parameters;
    // FitCellSurvival(cellSurvival_C4_2, "LM", 1);
    // FitCellSurvival(cellSurvival_C4_2, "LM", 2);
    FitCellSurvival(cellSurvival_C4_2, "LM", 3);
    // FitCellSurvival(cellSurvival_C4_2, "LM", 4);

};