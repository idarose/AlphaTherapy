#include "TFile.h"
#include "TTree.h"
#include <iostream>


void Analysis()
{
  // gROOT->Reset();
  // gROOT->SetStyle("Plain");

    std::unique_ptr<TFile> myFile(TFile::Open("../B4a-build/B4.root", "READ"));
    auto tree = myFile->Get<TTree>("B4");

    //std::vector<double> en = tree->GetEntries("EnergyDeps");


    TH1D *hEnergyDepositions = new TH1D("hEnergyDepositions", "Energy Deposition", 10000, 0.0, 20.0);
    TH1D *hEnergyDepositionsNucleus = new TH1D("hEnergyDepositionsNucleus", "Energy Deposition", 10000, 0.0, 20.0);

    tree->Draw("EnergyDeps>>hEnergyDepositions", "EnergyDeps>0.", "");

    tree->Draw("EnergyDeps>>hEnergyDepositionsNucleus", "VolumeTypes == 1", "");
    //h1->Draw();

    auto outputFile = new TFile("Output.root", "RECREATE");
    // auto outputFile2 = new TFile("Output2.root", "RECREATE");

    //----
    outputFile->cd();
    hEnergyDepositions->Write();
    hEnergyDepositionsNucleus->Write();

    //----
    // outputFile2->cd();
    // h2->Write();

    //----
    outputFile->Write();
    outputFile->Close();
    // h1->Fill(tree->GetEntries("EnergyDeps"));

  // int variable;
  // tree->SetBranchAddress("Energydeps", &variable);

  // for(int iEntry = 0; tree->LoadTree(iEntry) >=0; iEntry++)
  // {
  //   LoadTree->GetEntry(iEntry);
  //   printf("%d\n", variable);
  // }


}

