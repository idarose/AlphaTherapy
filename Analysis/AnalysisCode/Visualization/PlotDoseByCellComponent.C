#include <TH1.h>
#include <TCanvas.h>

// Function to find the first and last non-empty bins
std::pair<int, int> FindNonEmptyBinRange(TH1D* hist)
{
    const int nBins = hist->GetNbinsX();
    int firstNonEmptyBin = 1;
    int lastNonEmptyBin = nBins;

    // Find the first non-empty bin
    while(firstNonEmptyBin <= nBins && hist->GetBinContent(firstNonEmptyBin) == 0) {
        firstNonEmptyBin++;
    }

    // Find the last non-empty bin
    while(lastNonEmptyBin >= 1 && hist->GetBinContent(lastNonEmptyBin) == 0) {
        lastNonEmptyBin--;
    }

    return {firstNonEmptyBin, lastNonEmptyBin};
}

void MakePlots(std::string cellLine, std::string cellGeometry, int activity, int doseOrEnergy)
{
    std::vector<int> colours;
    colours.push_back(kCyan+2);
    colours.push_back(kGreen-2);
    colours.push_back(kOrange-3);
    colours.push_back(kViolet+1);
    colours.push_back(kRed);
    colours.push_back(kBlue);
    colours.push_back(kOrange+2);
    colours.push_back(kGreen+2);
    colours.push_back(kViolet);


    std::string histogramName_TotalCell;
    std::string histogramName_Nucleus;
    std::string histogramName_Membrane;
    std::string histogramName_Cytoplasm;

    std::string type;

    if(doseOrEnergy==0)
    {
        histogramName_TotalCell = "i0_hEnergyDeps_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_TotalCell_keVBinning";
        histogramName_Nucleus = "i0_hEnergyDeps_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_Nucleus_keVBinning";
        histogramName_Membrane = "i0_hEnergyDeps_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_Membrane_keVBinning";
        histogramName_Cytoplasm = "i0_hEnergyDeps_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_Cytoplasm_keVBinning";

    }
    if(doseOrEnergy==1)
    {
        histogramName_TotalCell = "i0_hDose_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_TotalCell_mGyBinning";
        histogramName_Nucleus = "i0_hDose_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_Nucleus_mGyBinning";
        histogramName_Membrane = "i0_hDose_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_Membrane_mGyBinning";
        histogramName_Cytoplasm = "i0_hDose_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_Cytoplasm_mGyBinning";
    }

    //------------------------------------------------
    std::string fileName = "../Output_" + cellGeometry + "/Output_"+ cellLine + "_" + std::to_string(activity) + "kBq.root";
    auto inputFile = std::unique_ptr<TFile>(TFile::Open(fileName.c_str()));


    //-----------------------------
    TH1D* histogram_TotalCell = 0;
    inputFile->GetObject(histogramName_TotalCell.c_str(), histogram_TotalCell);
    histogram_TotalCell->SetDirectory(0);

    TH1D* histogram_Membrane = 0;
    inputFile->GetObject(histogramName_Membrane.c_str(), histogram_Membrane);
    histogram_Membrane->SetDirectory(0);

    TH1D* histogram_Cytoplasm = 0;
    inputFile->GetObject(histogramName_Cytoplasm.c_str(), histogram_Cytoplasm);
    histogram_Cytoplasm->SetDirectory(0);

    TH1D* histogram_Nucleus = 0;
    inputFile->GetObject(histogramName_Nucleus.c_str(), histogram_Nucleus);
    histogram_Nucleus->SetDirectory(0);

    std::vector<TH1D*> hists = {histogram_TotalCell,histogram_Membrane,histogram_Cytoplasm,histogram_Nucleus};


    //--------------------------------
    std::string cellLine_Name;
    if(cellLine=="C4_2"){cellLine_Name = "C4-2";}
    if(cellLine=="PC3_PIP"){cellLine_Name = "PC3-PIP";}
    if(cellLine=="PC3_Flu"){cellLine_Name = "PC3-Flu";}

    int diameter;
    if(cellGeometry=="D12RP"||cellGeometry=="D12CP")
    {
        diameter = 12;
    }
    if(cellGeometry=="D5RP"||cellGeometry=="D5CP")
    {
        diameter = 5;
    }

    std::string nucleiDist;
    if(cellGeometry=="D12RP"||cellGeometry=="D5RP")
    {
        nucleiDist ="Uniformly";
    }
    if(cellGeometry=="D12CP"||cellGeometry=="D5CP")
    {
        nucleiDist ="Centrally";
    }

    std::string xAxisName;
    std::string yAxisName;
    if(doseOrEnergy==0)
    {
        xAxisName = "Energy Deposited [MeV]";
        if(cellLine=="PC3_Flu")
        {
            yAxisName = "Fraction of cells in sample / 10 keV bin";
        }
        else
        {
            yAxisName = "Fraction of cells in sample / 20 keV bin";
        }
    }
    if(doseOrEnergy==1)
    {
        xAxisName = "Dose Delivered [Gy]";

        if(cellLine=="PC3_Flu")
        {
            yAxisName = "Fraction of cells in sample / 1 mGy bin";
        }
        else
        {
            yAxisName = "Fraction of cells in sample / 2 mGy bin";
        }
    }

    std::string title = cellLine_Name + ", " + nucleiDist + " Dist. Nuclei, d_{nuc} = " + Form("%d", diameter) + " #mum, Activity = " + Form("%d", activity)   + "kBq / 1mL";

    //------------------------
    histogram_TotalCell->GetXaxis()->CenterTitle(true);
    histogram_TotalCell->GetYaxis()->CenterTitle(true);
    histogram_TotalCell->GetXaxis()->SetTitle(xAxisName.c_str());
    histogram_TotalCell->GetYaxis()->SetTitle(yAxisName.c_str());
    histogram_TotalCell->SetTitle(title.c_str());

    histogram_Nucleus->GetXaxis()->CenterTitle(true);
    histogram_Nucleus->GetYaxis()->CenterTitle(true);
    histogram_Nucleus->GetXaxis()->SetTitle(xAxisName.c_str());
    histogram_Nucleus->GetYaxis()->SetTitle(yAxisName.c_str());
    histogram_Nucleus->SetTitle(title.c_str());

    histogram_Membrane->GetXaxis()->CenterTitle(true);
    histogram_Membrane->GetYaxis()->CenterTitle(true);
    histogram_Membrane->GetXaxis()->SetTitle(xAxisName.c_str());
    histogram_Membrane->GetYaxis()->SetTitle(yAxisName.c_str());
    histogram_Membrane->SetTitle(title.c_str());

    histogram_Cytoplasm->GetXaxis()->CenterTitle(true);
    histogram_Cytoplasm->GetYaxis()->CenterTitle(true);
    histogram_Cytoplasm->GetXaxis()->SetTitle(xAxisName.c_str());
    histogram_Cytoplasm->GetYaxis()->SetTitle(yAxisName.c_str());
    histogram_Cytoplasm->SetTitle(title.c_str());

    //------------------------------
    // Create the canvas
    auto c1 = new TCanvas("c1", title.c_str(), 600,400);

    gStyle->SetOptStat(0);

    //------------------------------
    int reBin;
    if(cellLine=="PC3_Flu")
    {
        reBin = 1;
    }
    else
    {
        reBin = 2;
    }

    histogram_TotalCell->Rebin(reBin);
    histogram_Cytoplasm->Rebin(reBin);
    histogram_Nucleus->Rebin(reBin);
    histogram_Membrane->Rebin(reBin);


    //-------------------------
    double xMax = 0.;
    double xMin = 0.;
    double maxY;

    for(auto & entry : hists)
    {
        double Min;
        double Max;


        double maxF = entry->GetMaximum() + 0.05*entry->GetMaximum();
        if(maxF > maxY)
        {
            maxY = maxF;
        }

        auto [firstBin, lastBin] = FindNonEmptyBinRange(entry);
        if (firstBin <= lastBin) {
            Min = entry->GetBinLowEdge(firstBin);
            Max = entry->GetBinLowEdge(lastBin) + entry->GetBinWidth(lastBin);

        } else {
            std::cerr << "The histogram is empty!" << std::endl;
        }

        if(Max>xMax)
        {
            xMax = Max;
            xMin = Min;
        }
    }

    if(cellLine=="C4_2")
    {
        // maxY = 1./2.*maxY;
        if(diameter==12)
        {
            maxY = 1./2.*maxY;
        }
        else
        {
            maxY = 1./4.*maxY;
        }
    }

    if(cellLine=="PC3_Flu")
    {
        xMax = 1./4.*xMax;
    }



    histogram_TotalCell->GetXaxis()->SetRangeUser(0.,xMax);
    histogram_TotalCell->GetYaxis()->SetRangeUser(0.,maxY);

    histogram_Cytoplasm->GetXaxis()->SetRangeUser(0.,xMax);
    histogram_Cytoplasm->GetYaxis()->SetRangeUser(0.,maxY);

    histogram_Nucleus->GetXaxis()->SetRangeUser(0.,xMax);
    histogram_Nucleus->GetYaxis()->SetRangeUser(0.,maxY);

    histogram_Membrane->GetXaxis()->SetRangeUser(0.,xMax);
    histogram_Membrane->GetYaxis()->SetRangeUser(0.,maxY);


    histogram_TotalCell->SetLineColor(colours[0]);
    histogram_TotalCell->SetFillColorAlpha(colours[0], 0.5);

    histogram_Membrane->SetLineColor(colours[1]);
    histogram_Membrane->SetFillColorAlpha(colours[1], 0.5);

    histogram_Cytoplasm->SetLineColor(colours[2]);
    histogram_Cytoplasm->SetFillColorAlpha(colours[2], 0.5);

    histogram_Nucleus->SetLineColor(colours[3]);
    histogram_Nucleus->SetFillColorAlpha(colours[3], 0.5);


    histogram_TotalCell->Draw("HIST");
    histogram_Membrane->Draw("HIST SAME");
    histogram_Cytoplasm->Draw("HIST SAME");
    histogram_Nucleus->Draw("HIST SAME");

    auto legend = new TLegend(0.67,0.60,0.9,0.9);
    legend->SetHeader("Cell Component",title.c_str());
    TLegendEntry *header = (TLegendEntry*)legend->GetListOfPrimitives()->First();
    header->SetTextAlign(22);
    header->SetTextSize(.05);

    legend->AddEntry(histogram_TotalCell,"Total Cell")->SetTextSize(0.05);
    legend->AddEntry(histogram_Membrane,"Membrane")->SetTextSize(0.05);
    legend->AddEntry(histogram_Cytoplasm,"Cytoplasm")->SetTextSize(0.05);
    legend->AddEntry(histogram_Nucleus,"Nucleus")->SetTextSize(0.05);
    legend->Draw();

    std::string figureName;
    if(doseOrEnergy==0)
    {
        figureName =  "/Users/idarosenqvist/Desktop/Academics/MasterThesis/Thesis/figures/Results/" + cellGeometry + "/" + cellLine + "/EnergyDeposition_" + cellLine + "_" + std::to_string(activity) + "kBq.pdf";
    }
    if(doseOrEnergy==1)
    {
        figureName = "/Users/idarosenqvist/Desktop/Academics/MasterThesis/Thesis/figures/Results/" + cellGeometry + "/" + cellLine + "/DoseDelivered_" + cellLine + "_" + std::to_string(activity) + "kBq.pdf";
    }

    c1->Update();
    c1->SaveAs(figureName.c_str());


}

void MakePlotsLog(std::string cellLine, std::string cellGeometry, int activity, int doseOrEnergy)
{
    std::vector<int> colours;
    colours.push_back(kCyan+1);
    colours.push_back(kGreen-2);
    colours.push_back(kOrange-3);
    colours.push_back(kViolet+1);
    colours.push_back(kRed);
    colours.push_back(kBlue);
    colours.push_back(kOrange+2);
    colours.push_back(kGreen+2);
    colours.push_back(kViolet);


    std::string histogramName_TotalCell;
    std::string histogramName_Nucleus;
    std::string histogramName_Membrane;
    std::string histogramName_Cytoplasm;

    std::string type;

    if(doseOrEnergy==0)
    {
        histogramName_TotalCell = "i0_hEnergyDeps_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_TotalCell_keVBinning";
        histogramName_Nucleus = "i0_hEnergyDeps_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_Nucleus_keVBinning";
        histogramName_Membrane = "i0_hEnergyDeps_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_Membrane_keVBinning";
        histogramName_Cytoplasm = "i0_hEnergyDeps_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_Cytoplasm_keVBinning";

    }
    if(doseOrEnergy==1)
    {
        histogramName_TotalCell = "i0_hDose_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_TotalCell_mGyBinning";
        histogramName_Nucleus = "i0_hDose_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_Nucleus_mGyBinning";
        histogramName_Membrane = "i0_hDose_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_Membrane_mGyBinning";
        histogramName_Cytoplasm = "i0_hDose_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_Cytoplasm_mGyBinning";
    }

    //------------------------------------------------
    std::string fileName = "../Output_" + cellGeometry + "/Output_"+ cellLine + "_" + std::to_string(activity) + "kBq.root";
    auto inputFile = std::unique_ptr<TFile>(TFile::Open(fileName.c_str()));


    //-----------------------------
    TH1D* histogram_TotalCell = 0;
    inputFile->GetObject(histogramName_TotalCell.c_str(), histogram_TotalCell);
    histogram_TotalCell->SetDirectory(0);

    TH1D* histogram_Membrane = 0;
    inputFile->GetObject(histogramName_Membrane.c_str(), histogram_Membrane);
    histogram_Membrane->SetDirectory(0);

    TH1D* histogram_Cytoplasm = 0;
    inputFile->GetObject(histogramName_Cytoplasm.c_str(), histogram_Cytoplasm);
    histogram_Cytoplasm->SetDirectory(0);

    TH1D* histogram_Nucleus = 0;
    inputFile->GetObject(histogramName_Nucleus.c_str(), histogram_Nucleus);
    histogram_Nucleus->SetDirectory(0);

    std::vector<TH1D*> hists = {histogram_TotalCell,histogram_Membrane,histogram_Cytoplasm,histogram_Nucleus};


    //--------------------------------
    std::string cellLine_Name;
    if(cellLine=="C4_2"){cellLine_Name = "C4-2";}
    if(cellLine=="PC3_PIP"){cellLine_Name = "PC3-PIP";}
    if(cellLine=="PC3_Flu"){cellLine_Name = "PC3-Flu";}

    int diameter;
    if(cellGeometry=="D12RP"||cellGeometry=="D12CP")
    {
        diameter = 12;
    }
    if(cellGeometry=="D5RP"||cellGeometry=="D5CP")
    {
        diameter = 5;
    }

    std::string nucleiDist;
    if(cellGeometry=="D12RP"||cellGeometry=="D5RP")
    {
        nucleiDist ="Uniformly";
    }
    if(cellGeometry=="D12CP"||cellGeometry=="D5CP")
    {
        nucleiDist ="Centrally";
    }

    std::string xAxisName;
    std::string yAxisName;
    if(doseOrEnergy==0)
    {
        xAxisName = "Energy Deposited [MeV]";
        yAxisName = "Fraction of cells in sample / 20 keV bin";
    }
    if(doseOrEnergy==1)
    {
        xAxisName = "Dose Delivered [Gy]";
        yAxisName = "Fraction of cells in sample / 2 mGy bin";
    }

    std::string title = cellLine_Name + ", " + nucleiDist + " Dist. Nuclei, d_{nuc} = " + Form("%d", diameter) + " #mum, Activity = " + Form("%d", activity)   + "kBq / 1mL";

    // TGaxis::SetMaxDigits(6);

    //------------------------
    histogram_TotalCell->GetXaxis()->CenterTitle(true);
    histogram_TotalCell->GetYaxis()->CenterTitle(true);
    histogram_TotalCell->GetXaxis()->SetTitle(xAxisName.c_str());
    histogram_TotalCell->GetYaxis()->SetTitle(yAxisName.c_str());
    histogram_TotalCell->SetTitle(title.c_str());

    histogram_Nucleus->GetXaxis()->CenterTitle(true);
    histogram_Nucleus->GetYaxis()->CenterTitle(true);
    histogram_Nucleus->GetXaxis()->SetTitle(xAxisName.c_str());
    histogram_Nucleus->GetYaxis()->SetTitle(yAxisName.c_str());
    histogram_Nucleus->SetTitle(title.c_str());

    histogram_Membrane->GetXaxis()->CenterTitle(true);
    histogram_Membrane->GetYaxis()->CenterTitle(true);
    histogram_Membrane->GetXaxis()->SetTitle(xAxisName.c_str());
    histogram_Membrane->GetYaxis()->SetTitle(yAxisName.c_str());
    histogram_Membrane->SetTitle(title.c_str());

    histogram_Cytoplasm->GetXaxis()->CenterTitle(true);
    histogram_Cytoplasm->GetYaxis()->CenterTitle(true);
    histogram_Cytoplasm->GetXaxis()->SetTitle(xAxisName.c_str());
    histogram_Cytoplasm->GetYaxis()->SetTitle(yAxisName.c_str());
    histogram_Cytoplasm->SetTitle(title.c_str());


    //------------------------------
    // Create the canvas
    auto c1 = new TCanvas("c1", "Energy Depositions in Different Components of the Cell", 600,400);

    gStyle->SetOptStat(0);

    //------------------------------
    int reBin = 2;

    histogram_TotalCell->Rebin(reBin);
    histogram_Cytoplasm->Rebin(reBin);
    histogram_Nucleus->Rebin(reBin);
    histogram_Membrane->Rebin(reBin);


    //-------------------------
    double xMax = 0.;
    double xMin = 0.;
    double maxY;
    // double minY;

    for(auto & entry : hists)
    {
        double Min;
        double Max;


        double maxF = entry->GetMaximum() + 0.05*entry->GetMaximum();
        if(maxF > maxY)
        {
            maxY = maxF;
        }

        auto [firstBin, lastBin] = FindNonEmptyBinRange(entry);
        if (firstBin <= lastBin) {
            Min = entry->GetBinLowEdge(firstBin);
            Max = entry->GetBinLowEdge(lastBin) + entry->GetBinWidth(lastBin);

        } else {
            std::cerr << "The histogram is empty!" << std::endl;
        }

        if(Max>xMax)
        {
            xMax = Max;
            xMin = Min;
        }
    }


    histogram_TotalCell->GetXaxis()->SetRangeUser(0.,xMax);
    histogram_TotalCell->GetYaxis()->SetRangeUser(1.e-6,1.);

    histogram_Cytoplasm->GetXaxis()->SetRangeUser(0.,xMax);
    histogram_Cytoplasm->GetYaxis()->SetRangeUser(1.e-6,1.);

    histogram_Nucleus->GetXaxis()->SetRangeUser(0.,xMax);
    histogram_Nucleus->GetYaxis()->SetRangeUser(1.e-6,1.);

    histogram_Membrane->GetXaxis()->SetRangeUser(0.,xMax);
    histogram_Membrane->GetYaxis()->SetRangeUser(1.e-6,1.);


    histogram_TotalCell->SetLineColor(colours[0]);
    histogram_TotalCell->SetFillColorAlpha(colours[0], 0.3);

    histogram_Membrane->SetLineColor(colours[1]);
    histogram_Membrane->SetFillColorAlpha(colours[1], 0.3);

    histogram_Cytoplasm->SetLineColor(colours[2]);
    histogram_Cytoplasm->SetFillColorAlpha(colours[2], 0.3);

    histogram_Nucleus->SetLineColor(colours[3]);
    histogram_Nucleus->SetFillColorAlpha(colours[3], 0.3);


    histogram_TotalCell->Draw("HIST");
    histogram_Cytoplasm->Draw("HIST SAME");
    histogram_Nucleus->Draw("HIST SAME");
    histogram_Membrane->Draw("HIST SAME");

    auto legend = new TLegend(0.67,0.60,0.9,0.9);
    legend->SetHeader("Cell Component","C");
    TLegendEntry *header = (TLegendEntry*)legend->GetListOfPrimitives()->First();
    header->SetTextAlign(22);
    header->SetTextSize(.05);

    legend->AddEntry(histogram_TotalCell,"Total Cell")->SetTextSize(0.05);
    legend->AddEntry(histogram_Membrane,"Membrane")->SetTextSize(0.05);
    legend->AddEntry(histogram_Cytoplasm,"Cytoplasm")->SetTextSize(0.05);
    legend->AddEntry(histogram_Nucleus,"Nucleus")->SetTextSize(0.05);
    legend->Draw();


    std::string figureName;
    if(doseOrEnergy==0)
    {
        figureName =  "/Users/idarosenqvist/Desktop/Academics/MasterThesis/Thesis/figures/Results/" + cellGeometry + "/" + cellLine + "/EnergyDeposition_" + cellLine + "_" + std::to_string(activity) + "kBq_logy.pdf";
    }
    if(doseOrEnergy==1)
    {
        figureName = "/Users/idarosenqvist/Desktop/Academics/MasterThesis/Thesis/figures/Results/" + cellGeometry + "/" + cellLine + "/DoseDelivered_" + cellLine + "_" + std::to_string(activity) + "kBq_logy.pdf";
    }


    c1->SetLogy();

    c1->SetLogy();
    c1->SaveAs(figureName.c_str());


}

void PlotDoseByCellComponent()
{
    // std::string cellLine = "C4_2";
    std::string cellLine = "PC3_PIP";
    // std::string cellLine = "PC3_Flu";

    // MakePlots(cellLine, "D12RP", 25, 1);
    // MakePlotsLog(cellLine, "D12RP", 25, 1);

    // MakePlots(cellLine, "D12RP", 25, 0);
    // MakePlotsLog(cellLine, "D12RP", 25, 0);


    // MakePlots(cellLine, "D12CP", 25, 1);
    // MakePlotsLog(cellLine, "D12CP", 25, 1);

    MakePlots(cellLine, "D12CP", 25, 0);
    // MakePlotsLog(cellLine, "D12CP", 25, 0);

    // MakePlots(cellLine, "D5RP", 25, 1);
    // MakePlotsLog(cellLine, "D5RP", 25, 1);

    // MakePlots(cellLine, "D5RP", 25, 0);
    // MakePlotsLog(cellLine, "D5RP", 25, 0);

    // MakePlots(cellLine, "D5CP", 25, 1);
    // MakePlotsLog(cellLine, "D5CP", 25, 1);

    // MakePlots(cellLine, "D5CP", 25, 0);
    // MakePlotsLog(cellLine, "D5CP", 25, 0);

}