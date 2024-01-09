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
    std::vector<std::tuple<double, double, double>> data_cellSurvival_PC3_PIP;
    std::vector<std::tuple<double, double, double>> data_cellSurvival_PC3_Flu;


    C4-2
    data_cellSurvival_C4_2.push_back(make_tuple(0.0,    1.0,    0.0));
    data_cellSurvival_C4_2.push_back(make_tuple(5.0,    0.86,   0.117));
    data_cellSurvival_C4_2.push_back(make_tuple(10.0,    0.69,    0.014));
    data_cellSurvival_C4_2.push_back(make_tuple(25.0,    0.51,    0.096));
    data_cellSurvival_C4_2.push_back(make_tuple(50.0,    0.26,    0.035));
    data_cellSurvival_C4_2.push_back(make_tuple(75.0,    0.13,    0.031));
    data_cellSurvival_C4_2.push_back(make_tuple(100.0,    0.08,    0.035));
    data_cellSurvival_C4_2.push_back(make_tuple(150.0,    0.04,    0.024));


    // //PC3 Flu
    // // data_cellSurvival_C4_2.push_back(make_tuple(0.0,    1.0,    0.0));
    // // data_cellSurvival_C4_2.push_back(make_tuple(5.0,    0.86,   0.117));
    // data_cellSurvival_C4_2.push_back(make_tuple(10.0,    0.955,    0.0955));
    // data_cellSurvival_C4_2.push_back(make_tuple(25.0,    0.724,    0.0724));
    // data_cellSurvival_C4_2.push_back(make_tuple(50.0,    0.733,    0.0733));
    // data_cellSurvival_C4_2.push_back(make_tuple(75.0,    0.798,    0.0798));
    // data_cellSurvival_C4_2.push_back(make_tuple(100.0,    0.720,    0.0720));
    // data_cellSurvival_C4_2.push_back(make_tuple(150.0,    0.680,    0.0698));

    // //PC3 PIP
    // // data_cellSurvival_C4_2.push_back(make_tuple(0.0,    1.0,    0.0));
    // // data_cellSurvival_C4_2.push_back(make_tuple(5.0,    0.86,   0.117));
    // data_cellSurvival_C4_2.push_back(make_tuple(10.0,    0.630,    0.063));
    // data_cellSurvival_C4_2.push_back(make_tuple(25.0,    0.317,    0.05));
    // data_cellSurvival_C4_2.push_back(make_tuple(50.0,    0.071,    0.05));
    // data_cellSurvival_C4_2.push_back(make_tuple(75.0,    0.032,    0.05));
    // data_cellSurvival_C4_2.push_back(make_tuple(100.0,    0.014,    0.05));
    // data_cellSurvival_C4_2.push_back(make_tuple(150.0,    0.006,    0.05));




    //------------------------------------------------------------------------------------------------
    //      Construct the TGraphErrors object corresponding to the biological data
    TGraphErrors *gr_clonogenicSurvival_C4_2 = new TGraphErrors();
    gr_clonogenicSurvival_C4_2->SetName("gr_clonogenicSurvival_C4_2");

    // TGraphErrors *gr_clonogenicSurvival_PC3_PIP = new TGraphErrors();
    // gr_clonogenicSurvival_PC3_PIP->SetName("gr_clonogenicSurvival_PC3_PIP");

    // TGraphErrors *gr_clonogenicSurvival_PC3_Flu = new TGraphErrors();
    // gr_clonogenicSurvival_PC3_Flu->SetName("gr_clonogenicSurvival_PC3_Flu");

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

    GraphData(gr_clonogenicSurvival_C4_2, data_cellSurvival_C4_2);
    // GraphData(gr_clonogenicSurvival_PC3_PIP, data_cellSurvival_PC3_PIP);
    // GraphData(gr_clonogenicSurvival_PC3_Flu, data_cellSurvival_PC3_Flu);

    // TGraphErrors *gr_clonogenicSurvival_PC3_PIP = new TGraphErrors();
    // gr_clonogenicSurvival_PC3_PIP->SetName("gr_clonogenicSurvival_PC3_PIP");

    // int graphPointN_data = 0;
    // for(auto & entry : data_cellSurvival_PC3_PIP)
    // {
    //     double activity = std::get<0>(entry);
    //     double cellSurvival = std::get<1>(entry);
    //     double cellSurvivalUncertainty = std::get<2>(entry);

    //     gr_clonogenicSurvival_PC3_PIP->SetPoint(graphPointN_data, activity, cellSurvival);
    //     gr_clonogenicSurvival_PC3_PIP->SetPointError(graphPointN_data, 0.0, cellSurvivalUncertainty);

    //     graphPointN_data++;
    // }

    // TGraphErrors *gr_clonogenicSurvival_PC3_Flu = new TGraphErrors();
    // gr_clonogenicSurvival_PC3_Flu->SetName("gr_clonogenicSurvival_PC3_Flu");

    // int graphPointN_data = 0;
    // for(auto & entry : data_cellSurvival_PC3_Flu)
    // {
    //     double activity = std::get<0>(entry);
    //     double cellSurvival = std::get<1>(entry);
    //     double cellSurvivalUncertainty = std::get<2>(entry);

    //     gr_clonogenicSurvival_PC3_Flu->SetPoint(graphPointN_data, activity, cellSurvival);
    //     gr_clonogenicSurvival_PC3_Flu->SetPointError(graphPointN_data, 0.0, cellSurvivalUncertainty);

    //     graphPointN_data++;
    // }

    //------------------------------------------------------------------------------------------------
    //      Here, you shall read in the different histograms of energy depositions.
    //      These histograms are produced from your analysis of your GEANT4 simulation
    // TFile *inputFile_5kBq = new TFile("../../OutputAnalysisCode/Output_C4_2_5kBq.root", "READ");
    // TH1D* hEnergyDeps_212Pb_C4_2_5kBq_Nucleus = nullptr;
    // inputFile_5kBq->GetObject("hEnergyDeps_212Pb_C4_2_5kBq_Nucleus", hEnergyDeps_212Pb_C4_2_5kBq_Nucleus);
    // hEnergyDeps_212Pb_C4_2_5kBq_Nucleus->SetDirectory(0);
    // inputFile_5kBq->Close();

    TFile *inputFile_10kBq = new TFile("../../OutputAnalysisCode/Output_PC3_PIP_10kBq.root", "READ");
    TH1D* hEnergyDeps_212Pb_C4_2_10kBq_Nucleus = nullptr;
    inputFile_10kBq->GetObject("hEnergyDeps_212Pb_PC3_PIP_10kBq_Nucleus", hEnergyDeps_212Pb_C4_2_10kBq_Nucleus);
    hEnergyDeps_212Pb_C4_2_10kBq_Nucleus->SetDirectory(0);
    inputFile_10kBq->Close();

    TFile *inputFile_25kBq = new TFile("../../OutputAnalysisCode/Output_PC3_PIP_25kBq.root", "READ");
    TH1D* hEnergyDeps_212Pb_C4_2_25kBq_Nucleus = nullptr;
    inputFile_25kBq->GetObject("hEnergyDeps_212Pb_PC3_PIP_25kBq_Nucleus", hEnergyDeps_212Pb_C4_2_25kBq_Nucleus);
    hEnergyDeps_212Pb_C4_2_25kBq_Nucleus->SetDirectory(0);
    inputFile_25kBq->Close();

    TFile *inputFile_50kBq = new TFile("../../OutputAnalysisCode/Output_PC3_PIP_50kBq.root", "READ");
    TH1D* hEnergyDeps_212Pb_C4_2_50kBq_Nucleus = nullptr;
    inputFile_50kBq->GetObject("hEnergyDeps_212Pb_PC3_PIP_50kBq_Nucleus", hEnergyDeps_212Pb_C4_2_50kBq_Nucleus);
    hEnergyDeps_212Pb_C4_2_50kBq_Nucleus->SetDirectory(0);
    inputFile_50kBq->Close();

    TFile *inputFile_75kBq = new TFile("../../OutputAnalysisCode/Output_PC3_PIP_75kBq.root", "READ");
    TH1D* hEnergyDeps_212Pb_C4_2_75kBq_Nucleus = nullptr;
    inputFile_75kBq->GetObject("hEnergyDeps_212Pb_PC3_PIP_75kBq_Nucleus", hEnergyDeps_212Pb_C4_2_75kBq_Nucleus);
    hEnergyDeps_212Pb_C4_2_75kBq_Nucleus->SetDirectory(0);
    inputFile_75kBq->Close();

    TFile *inputFile_100kBq = new TFile("../../OutputAnalysisCode/Output_PC3_PIP_100kBq.root", "READ");
    TH1D* hEnergyDeps_212Pb_C4_2_100kBq_Nucleus = nullptr;
    inputFile_100kBq->GetObject("hEnergyDeps_212Pb_PC3_PIP_100kBq_Nucleus", hEnergyDeps_212Pb_C4_2_100kBq_Nucleus);
    hEnergyDeps_212Pb_C4_2_100kBq_Nucleus->SetDirectory(0);
    inputFile_100kBq->Close();

    TFile *inputFile_150kBq = new TFile("../../OutputAnalysisCode/Output_PC3_PIP_150kBq.root", "READ");
    TH1D* hEnergyDeps_212Pb_C4_2_150kBq_Nucleus = nullptr;
    inputFile_150kBq->GetObject("hEnergyDeps_212Pb_PC3_PIP_150kBq_Nucleus", hEnergyDeps_212Pb_C4_2_150kBq_Nucleus);
    hEnergyDeps_212Pb_C4_2_150kBq_Nucleus->SetDirectory(0);
    inputFile_150kBq->Close();




    // auto hEnergyDeps_212Pb_C4_2_0kBq_Nucleus = (TH1D*)hEnergyDeps_212Pb_C4_2_10kBq_Nucleus->Clone();
    // hEnergyDeps_212Pb_C4_2_0kBq_Nucleus->Reset();
    // hEnergyDeps_212Pb_C4_2_0kBq_Nucleus->SetName("hEnergyDeps_212Pb_PC3_PIP_0kBq_Nucleus");


    //------------------------------------------------------------------------------------------------
    //      This is just to make some dummy histograms for different activities
    //      You will not need this block of code when you have the histograms for each activity.

    std::vector<std::tuple<double, TH1D*>> vec_activities_histograms;

    // vec_activities_histograms.push_back(std::make_tuple(0.0, hEnergyDeps_212Pb_C4_2_0kBq_Nucleus));
    // vec_activities_histograms.push_back(std::make_tuple(5.0, hEnergyDeps_212Pb_C4_2_5kBq_Nucleus));
    vec_activities_histograms.push_back(std::make_tuple(10.0, hEnergyDeps_212Pb_C4_2_10kBq_Nucleus));
    vec_activities_histograms.push_back(std::make_tuple(25.0, hEnergyDeps_212Pb_C4_2_25kBq_Nucleus));
    vec_activities_histograms.push_back(std::make_tuple(50.0, hEnergyDeps_212Pb_C4_2_50kBq_Nucleus));
    vec_activities_histograms.push_back(std::make_tuple(75.0, hEnergyDeps_212Pb_C4_2_75kBq_Nucleus));
    vec_activities_histograms.push_back(std::make_tuple(100.0, hEnergyDeps_212Pb_C4_2_100kBq_Nucleus));
    vec_activities_histograms.push_back(std::make_tuple(150.0, hEnergyDeps_212Pb_C4_2_150kBq_Nucleus));


    for(auto & entry : vec_activities_histograms)
    {
        double activity = std::get<0>(entry);
        auto hist = std::get<1>(entry);

        double integral = hist->Integral();

        std::cout << std::setprecision(15) << activity << ",\t integral: " << integral << std::endl;
        // std::cout << "activity: " << activity << ",\t integral: " << integral << " mean: " << hist->GetMean() << std::endl;
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
        double fractionOfComponentsHit = h_energyDeposition->Integral(0.,26.);
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
            double cellSurvival = gr_cellSurvivability_vs_activitykBqPerMl.Eval(x[0], 0, "");

            return cellSurvival;

        }, 0.0, 150.0, nParameters);

    f_cellSurvivalVsDose_C4_2->SetParLimits(0, 0.0, 1.0e+05);
    f_cellSurvivalVsDose_C4_2->SetParLimits(1, 0.0, 1.0e+05);

    f_cellSurvivalVsDose_C4_2->SetNpx(10000);
    f_cellSurvivalVsDose_C4_2->SetParameter(0, 1.0e+00);
    f_cellSurvivalVsDose_C4_2->SetParameter(1, 1.0e+00);

    gr_clonogenicSurvival_C4_2->Fit(f_cellSurvivalVsDose_C4_2, "", "", 0.0, 150.0);


    //------------------------------------------------------------------------------------------------
    TFile *outputFile = new TFile("Output_AnalyseCellSurvival_PC3_PC3.root", "RECREATE");

    gr_clonogenicSurvival_C4_2->Write();
    f_cellSurvivalVsDose_C4_2->Write();
    gr_cellSurvivability_vs_activitykBqPerMl.Write();

    for(auto & entry : vec_activities_histograms)
    {
        auto hist = std::get<1>(entry);

        hist->Write();
    }

    outputFile->Write();
}






















