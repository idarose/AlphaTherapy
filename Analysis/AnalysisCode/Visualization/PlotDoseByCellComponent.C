

void PlotHistogramsDoseDifferentComponentsCell(std::string cellLine, int activity)
{

    // i0_hDose_212Pb_C4_2_5kBq_Cytoplasm_uGyBinning


    //------------------------------------------------
    std::string fileName = "../../OutputAnalysisCode/Output_" + cellLine + "_" + std::to_string(activity) + "kBq.root";
    auto inputFile = std::unique_ptr<TFile>(TFile::Open(fileName.c_str()));

    std::string histogramName_Dose_TotalCell = "i0_hEnergyDeps_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_TotalCell_keVBinning";
    std::string histogramName_Dose_Nucleus = "i0_hEnergyDeps_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_Nucleus_keVBinning";
    std::string histogramName_Dose_Membrane = "i0_hEnergyDeps_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_Membrane_keVBinning";
    std::string histogramName_Dose_Cytoplasm = "i0_hEnergyDeps_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_Cytoplasm_keVBinning";

    //--------------------------------
    std::string OutputFile = "Output_PlotHistograms_DoseDifferentCellComponent_" + cellLine + "_"  + std::to_string(activity) + "kBq.root";
    auto output = new TFile(OutputFile.c_str(), "RECREATE");

    //-----------------------------
    TH1D* histogram_Dose_TotalCell = 0;
    inputFile->GetObject(histogramName_Dose_TotalCell.c_str(), histogram_Dose_TotalCell);
    histogram_Dose_TotalCell->SetDirectory(0);

    TH1D* histogram_Dose_Membrane = 0;
    inputFile->GetObject(histogramName_Dose_Membrane.c_str(), histogram_Dose_Membrane);
    histogram_Dose_Membrane->SetDirectory(0);

    TH1D* histogram_Dose_Cytoplasm = 0;
    inputFile->GetObject(histogramName_Dose_Cytoplasm.c_str(), histogram_Dose_Cytoplasm);
    histogram_Dose_Cytoplasm->SetDirectory(0);

    TH1D* histogram_Dose_Nucleus = 0;
    inputFile->GetObject(histogramName_Dose_Nucleus.c_str(), histogram_Dose_Nucleus);
    histogram_Dose_Nucleus->SetDirectory(0);

    //------------------------------
    // Create the canvas
    auto c1 = new TCanvas("c1", "Energy Depositions in Different Components of the Cell", 600, 400);


    histogram_Dose_TotalCell->SetLineColor(kMagenta);
    histogram_Dose_TotalCell->SetFillColor(kMagenta);

    histogram_Dose_Membrane->SetLineColor(kRed);
    histogram_Dose_Membrane->SetFillColor(kRed);

    histogram_Dose_Cytoplasm->SetLineColor(kBlue);
    histogram_Dose_Cytoplasm->SetFillColor(kBlue);

    histogram_Dose_Nucleus->SetLineColor(kGreen);
    histogram_Dose_Nucleus->SetFillColor(kGreen);


    histogram_Dose_TotalCell->Draw("HIST");
    histogram_Dose_Cytoplasm->Draw("HIST SAME");
    histogram_Dose_Nucleus->Draw("HIST SAME");
    histogram_Dose_Membrane->Draw("HIST SAME");

    auto legend = new TLegend(0.1,0.7,0.48,0.9);
    legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
    legend->AddEntry(histogram_Dose_TotalCell,"Total Cell");
    legend->AddEntry(histogram_Dose_Membrane,"Membrane");
    legend->AddEntry(histogram_Dose_Cytoplasm,"Cytoplasm");
    legend->AddEntry(histogram_Dose_Nucleus,"Nucleus");
    legend->Draw();


    std::string figureName = "DoseInDifferentComponents_" + std::to_string(activity) + "kBq_" + "ByOriginDecay.pdf";

    c1->Update();
    c1->SaveAs(figureName.c_str());

    histogram_Dose_Cytoplasm->Write();
    histogram_Dose_TotalCell->Write();
    histogram_Dose_Membrane->Write();
    histogram_Dose_Nucleus->Write();

    output->Write();
    output->Close();
}