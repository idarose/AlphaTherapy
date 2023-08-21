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


void AnalyseCellSurvival()
{
    //------------------------------------------------
    //      Define the biological data <activity_kBqPerMl, cellSurvival, cellSurvivalUncertainty_stdDev>
    std::vector<std::tuple<double, double, double>> data_cellSurvival_C4_2;
    data_cellSurvival_C4_2.push_back(make_tuple(0.0,    1.0,    0.0));
    data_cellSurvival_C4_2.push_back(make_tuple(5.0,    0.86,   0.117));
    data_cellSurvival_C4_2.push_back(make_tuple(10.0,   0.69,   0.014));
    data_cellSurvival_C4_2.push_back(make_tuple(25.0,   0.51,   0.096));
    data_cellSurvival_C4_2.push_back(make_tuple(50.0,   0.26,   0.035));
    data_cellSurvival_C4_2.push_back(make_tuple(75.0,   0.13,   0.031));
    data_cellSurvival_C4_2.push_back(make_tuple(100.0,  0.08,   0.035));
    data_cellSurvival_C4_2.push_back(make_tuple(150.0,  0.04,   0.024));

    //------------------------------------------------------------------------------------------------
    //      Construct the TGraphErrors object corresponding to the biological data
    TGraphErrors *gr_clonogenicSurvival_C4_2 = new TGraphErrors();
    gr_clonogenicSurvival_C4_2->SetName("gr_clonogenicSurvival_C4_2");

    int graphPointN_data = 0;
    for(auto & entry : data_cellSurvival_C4_2)
    {
        double activity = std::get<0>(entry);
        double cellSurvival = std::get<1>(entry);
        double cellSurvivalUncertainty = std::get<2>(entry);

        gr_clonogenicSurvival_C4_2->SetPoint(graphPointN_data, activity, cellSurvival);
        gr_clonogenicSurvival_C4_2->SetPointError(graphPointN_data, 0.0, cellSurvivalUncertainty);

        graphPointN_data++;
    }

    //------------------------------------------------------------------------------------------------
    //      Here, you shall read in the different histograms of energy depositions.
    //      These histograms are produced from your analysis of your GEANT4 simulation
    TFile *inputFile = new TFile("outputMainAnalysisCode.root", "READ");

    TH1D* hEnergyDeps = nullptr;
    inputFile->GetObject("hEnergyDepsNucleus", hEnergyDeps);
    // inputFile->GetObject("hEnergyDepsCytoplasm", hEnergyDeps);

    hEnergyDeps->SetDirectory(0);

    inputFile->Close();

    //------------------------------------------------------------------------------------------------
    //      This is just to make some dummy histograms for different activities
    //      You will not need this block of code when you have the histograms for each activity.

    hEnergyDeps->SetBinContent(1, 0.0);
    // hEnergyDeps->Scale(1.0/15000.0);

    double integral = hEnergyDeps->Integral();
    std::cout << "integral: " << integral << std::endl;

    double globalScaling_histograms = 9.5e+01/2.3/100.0; // This is an arbitrary scaling factor so that for 150 kBq, almost all cells are hit

    std::vector<std::tuple<double, TH1D*>> vec_activities_histograms;

    auto hEnergyDeps_0kBq = (TH1D*)hEnergyDeps->Clone();
    hEnergyDeps_0kBq->Reset();
    hEnergyDeps_0kBq->SetName("hEnergyDeps_0kBq");

    auto hEnergyDeps_5kBq = (TH1D*)hEnergyDeps->Clone();
    hEnergyDeps_5kBq->Scale(globalScaling_histograms*5.0/10.0);
    hEnergyDeps_5kBq->SetName("hEnergyDeps_5kBq");

    auto hEnergyDeps_10kBq = (TH1D*)hEnergyDeps->Clone();
    hEnergyDeps_10kBq->Scale(globalScaling_histograms*10.0/10.0);
    hEnergyDeps_10kBq->SetName("hEnergyDeps_10kBq");

    auto hEnergyDeps_25kBq = (TH1D*)hEnergyDeps->Clone();
    hEnergyDeps_25kBq->Scale(globalScaling_histograms*25.0/10.0);
    hEnergyDeps_25kBq->SetName("hEnergyDeps_25kBq");

    auto hEnergyDeps_50kBq = (TH1D*)hEnergyDeps->Clone();
    hEnergyDeps_50kBq->Scale(globalScaling_histograms*50.0/10.0);
    hEnergyDeps_50kBq->SetName("hEnergyDeps_50kBq");

    auto hEnergyDeps_75kBq = (TH1D*)hEnergyDeps->Clone();
    hEnergyDeps_75kBq->Scale(globalScaling_histograms*75.0/10.0);
    hEnergyDeps_75kBq->SetName("hEnergyDeps_75kBq");

    auto hEnergyDeps_100kBq = (TH1D*)hEnergyDeps->Clone();
    hEnergyDeps_100kBq->Scale(globalScaling_histograms*100.0/10.0);
    hEnergyDeps_100kBq->SetName("hEnergyDeps_100kBq");

    auto hEnergyDeps_150kBq = (TH1D*)hEnergyDeps->Clone();
    hEnergyDeps_150kBq->Scale(globalScaling_histograms*150.0/10.0);
    hEnergyDeps_150kBq->SetName("hEnergyDeps_150kBq");

    vec_activities_histograms.push_back(std::make_tuple(0.0, hEnergyDeps_0kBq));
    vec_activities_histograms.push_back(std::make_tuple(5.0, hEnergyDeps_5kBq));
    vec_activities_histograms.push_back(std::make_tuple(10.0, hEnergyDeps_10kBq));
    vec_activities_histograms.push_back(std::make_tuple(25.0, hEnergyDeps_25kBq));
    vec_activities_histograms.push_back(std::make_tuple(50.0, hEnergyDeps_50kBq));
    vec_activities_histograms.push_back(std::make_tuple(75.0, hEnergyDeps_75kBq));
    vec_activities_histograms.push_back(std::make_tuple(100.0, hEnergyDeps_100kBq));
    vec_activities_histograms.push_back(std::make_tuple(150.0, hEnergyDeps_150kBq));

    for(auto & entry : vec_activities_histograms)
    {
        double activity = std::get<0>(entry);
        auto hist = std::get<1>(entry);

        double integral = hist->Integral();

        std::cout << "activity: " << activity << ",\t integral: " << integral << std::endl;
    }

    //------------------------------------------------------------------------------------------------
    //      This is the function which calculates the cell survival (fraction) for a particular
    //      In terms of programming, this is a lambda function/expression which captures by reference "[&]" (I shall explain).
    //      You might have some experience with this already—it's a very useful feature of >=C++11
    auto CalculateCellSurvivalFraction = [&](TH1D *h_energyDeposition, double *par)
    {
        double alpha = par[0];
        double beta = par[1];

        // std::cout << "alpha: " << alpha << std::endl;
        // std::cout << "beta: " << beta << std::endl;

        //------------------------------------------------
        double fractionOfComponentsHit = h_energyDeposition->Integral();
        double fractionOfComponentsMissed = (1.0 - fractionOfComponentsHit);

        //------------------------------------------------
        //      The fraction of missed cells obviously all survive so that is immediately added.
        double fractionOfTotalCellsSurviving = fractionOfComponentsMissed;

        for(int i=0; i<h_energyDeposition->GetNbinsX(); i++)
        {
            double energyDeposition = h_energyDeposition->GetBinCenter(i+1);
            double fractionOfTotalCells = h_energyDeposition->GetBinContent(i+1);

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
    auto GenerateGraph_CellSurvivalFraction = [&](std::vector<std::tuple<double, TH1D*>> vec, double *par)
    {
        TGraphErrors gr;

        int graphPointN = 0;

        for(auto & entry : vec)
        {
            double activity = std::get<0>(entry);
            auto hist = std::get<1>(entry);

            double cellSurvival = CalculateCellSurvivalFraction(hist, par);

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
            double cellSurvival = gr_cellSurvivability_vs_activitykBqPerMl.Eval(x[0], 0, "S");

            return cellSurvival;

        }, 0.0, 150.0, nParameters);

    f_cellSurvivalVsDose_C4_2->SetParLimits(0, 0.0, 1.0e+05);
    f_cellSurvivalVsDose_C4_2->SetParLimits(1, 0.0, 1.0e+05);

    f_cellSurvivalVsDose_C4_2->SetNpx(10000);
    f_cellSurvivalVsDose_C4_2->SetParameter(0, 1.0e+00);
    f_cellSurvivalVsDose_C4_2->SetParameter(1, 1.0e+00);

    gr_clonogenicSurvival_C4_2->Fit(f_cellSurvivalVsDose_C4_2, "", "", 0.0, 150.0);


    //------------------------------------------------------------------------------------------------
    TFile *outputFile = new TFile("Output_AnalyseCellSurvival.root", "RECREATE");

    gr_clonogenicSurvival_C4_2->Write();
    f_cellSurvivalVsDose_C4_2->Write();

    for(auto & entry : vec_activities_histograms)
    {
        auto hist = std::get<1>(entry);

        hist->Write();
    }

    outputFile->Write();
}






















