#include <tuple>
#include <vector>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TMath.h>
#include <TChain.h>
#include <TGraph.h>

#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>
#include <ctime>
#include <tuple>
#include <future>
#include <thread>

void MakeBraggPeakCurves()
{
    TFile *fileRun1 = TFile::Open("BraggPeakSimulation-build/Output_Run1.root");
    auto treeRun1= fileRun1->Get<TTree>("B4");
    TTreeReader myReaderRun1(treeRun1);


    // TFile *fileRun2 = TFile::Open("BraggPeakSimulation-build/Output_Run2.root");
    // auto treeRun2= fileRun2->Get<TTree>("B4");
    // TTreeReader myReaderRun2(treeRun2);


    TTreeReaderArray<double> energyDepRun1(myReaderRun1, "EnergyDep");
    TTreeReaderArray<double> positionZRun1(myReaderRun1, "PositionZ");
    TTreeReaderArray<double> stepLengthRun1(myReaderRun1, "StepLength");


    // TTreeReaderArray<double> energyDepRun2(myReaderRun2, "EnergyDep");
    // TTreeReaderArray<double> positionZRun2(myReaderRun2, "PositionZ");
    // TTreeReaderArray<double> stepLengthRun2(myReaderRun2, "StepLength");



    int NBins = 1000;
    TH2D* histRun1 = new TH2D("hEnergyLoss_Lengths_Run1", "Energy loss for lengths travelled", NBins, 0., 2., NBins, 0.,1000.);
    histRun1->GetXaxis()->SetTitle("Energy Deposition [MeV");
    histRun1->GetYaxis()->SetTitle("Position z [um]");

    TH1D* histMaxZ = new TH1D("hMaxZPosition", "Maximum z position", NBins, 0., 1000.);
    histMaxZ->GetXaxis()->SetTitle("Position z [um]");

    TH1D* histEnergyDep = new TH1D("hEnergyDep", "Energy deposition for range", NBins, 0., 1000.);
    histMaxZ->GetXaxis()->SetTitle("Position z [um]");

    std::cout << 100./((double)NBins) << std::endl;

    int scale = 0;
    // std::vector<std::tuple<double,double>> vec;



    while(myReaderRun1.Next())
    {
        double maxZ = 0.;
        // if(scale>1)
        // {
        //     break;
        // }

        double z = 0.;
        double dz = 1.;

        for(int i=0; i<101; i++)
        {
            z = i;

            double energyAtZ = 0;

            for(int i=0; i<energyDepRun1.GetSize(); i++)
            {
                double step = stepLengthRun1[i];
                double positionZ = positionZRun1[i];

                if(positionZ>(z-dz) && positionZ < z)
                {
                    double energyDep = energyDepRun1[i];

                    double S = energyDep*1000.;

                    energyAtZ+=S;
                }

                // histRun1->Fill(S,positionZ);
                // scale++;

            }
            histEnergyDep->Fill(z, energyAtZ);
        }

        scale++;

        histMaxZ->Fill(maxZ);
    }
    histEnergyDep->Scale(1./((double)scale));

    // histEnergyDep->Rebin(10);

    TH1D* projectionY = histRun1->ProjectionY("projectionY", 1, histRun1->GetNbinsX());
    projectionY->SetDirectory(0);
    projectionY->SetName("hProj");

    projectionY->Scale(1./((double)scale));
    projectionY->Rebin(10);

    TGraph* grRun1 = new TGraph();
    grRun1->SetName("grRun1");
    grRun1->SetTitle("Energy Loss of 6.09 MeV Alpha");
    grRun1->GetYaxis()->SetTitle("Energy Loss [MeV]");
    grRun1->GetXaxis()->SetTitle("Length Travelled z [um]");

    int NPoint = 0;
    for(int i=0; i<histEnergyDep->GetNbinsX(); i++)
    {
        double len = (histEnergyDep->GetXaxis())->GetBinCenter(i+1);
        // TH1D* projected = histRun1->ProjectionX("projRun1",i+1,i+1);
        // double meanS = projected->GetMean();

        // projected->SetDirectory(0);
        // // double len = std::get<0>(vec[i]);
        // // double meanS = std::get<1>(vec[i]);

        // grRun1->SetPoint(NPoint, len, meanS);

        double energyDep = histEnergyDep->GetBinContent(i+1);


        // std::cout << len << " . " << energyDep << std::endl;
        grRun1->SetPoint(NPoint,len,energyDep);

        NPoint++;
    }

    // TH2D* histRun2 = new TH2D("hEnergyLoss_Lengths_Run2", "Energy loss for lengths travelled", NBins, 0., 20., NBins, 0.,20.);
    // histRun1->GetXaxis()->SetTitle("Energy loss [MeV");
    // histRun1->GetYaxis()->SetTitle("Lengths [cm]");

    // scale = 0;
    // while(myReaderRun2.Next())
    // {
    //     double totalStep = 0.;

    //     for(int i=0; i<energyLossRun2.GetSize(); i++)
    //     {
    //         totalStep += stepLengthRun2[i];
    //         double energyloss = energyLossRun2[i];

    //         histRun2->Fill(energyloss,totalStep);

    //     }
    //     scale++;
    // }
    // histRun2->Scale(((double)scale));

    std::string outputName = "Output_Bragg_1_electrons.root";
    auto output = new TFile(outputName.c_str(), "RECREATE");
    histRun1->Write();
    projectionY->Write();
    // histRun2->Write();
    histMaxZ->Write();
    histEnergyDep->Write();
    grRun1->Write();
    output->Write();
    output->Close();

}