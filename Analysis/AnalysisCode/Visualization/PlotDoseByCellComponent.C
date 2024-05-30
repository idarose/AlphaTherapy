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
        if(cellLine=="PC3_Flu")
        {
            histogramName_TotalCell = "i0_hEnergyDeps_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_TotalCell_eVBinning";
            histogramName_Nucleus = "i0_hEnergyDeps_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_Nucleus_eVBinning";
            histogramName_Membrane = "i0_hEnergyDeps_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_Membrane_eVBinning";
            histogramName_Cytoplasm = "i0_hEnergyDeps_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_Cytoplasm_eVBinning";
        }
        else
        {
            histogramName_TotalCell = "i0_hEnergyDeps_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_TotalCell_keVBinning";
            histogramName_Nucleus = "i0_hEnergyDeps_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_Nucleus_keVBinning";
            histogramName_Membrane = "i0_hEnergyDeps_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_Membrane_keVBinning";
            histogramName_Cytoplasm = "i0_hEnergyDeps_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_Cytoplasm_keVBinning";
        }

    }
    if(doseOrEnergy==1)
    {
        if(cellLine=="PC3_Flu")
        {
            histogramName_TotalCell = "i0_hDose_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_TotalCell_uGyBinning";
            histogramName_Nucleus = "i0_hDose_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_Nucleus_uGyBinning";
            histogramName_Membrane = "i0_hDose_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_Membrane_uGyBinning";
            histogramName_Cytoplasm = "i0_hDose_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_Cytoplasm_uGyBinning";
        }
        else
        {
            histogramName_TotalCell = "i0_hDose_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_TotalCell_mGyBinning";
            histogramName_Nucleus = "i0_hDose_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_Nucleus_mGyBinning";
            histogramName_Membrane = "i0_hDose_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_Membrane_mGyBinning";
            histogramName_Cytoplasm = "i0_hDose_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_Cytoplasm_mGyBinning";
        }
    }

    //------------------------------------------------
    std::string fileName = "../Output_" + cellGeometry + "/Output_"+ cellLine + "_" + std::to_string(activity) + "kBq.root";
    auto inputFile = std::unique_ptr<TFile>(TFile::Open(fileName.c_str()));

    double titleSize = 0.045;

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
    std::string generalTitle;
    if(doseOrEnergy==0)
    {
        xAxisName = "Energy Deposited [MeV]";
        generalTitle = "Energy Deposition";
        if(cellLine=="PC3_Flu")
        {
            yAxisName = "Fraction of cells in sample / 2 keV bin";
        }
        else
        {
            yAxisName = "Fraction of cells in sample / 40 keV bin";
        }
    }
    if(doseOrEnergy==1)
    {
        xAxisName = "Dose Delivered [Gy]";
        generalTitle = "Absorbed Dose";
        if(cellLine=="PC3_Flu")
        {
            yAxisName = "Fraction of cells in sample / 1 mGy bin";
        }
        else
        {
            yAxisName = "Fraction of cells in sample / 4 mGy bin";
        }
    }

    std::string title = cellLine_Name + ", " + nucleiDist + " Dist. Nuclei, d_{nuc} = " + Form("%d", diameter) + " #mum, " + Form("%d", activity)   + " kBq / mL";

    //------------------------
    histogram_TotalCell->GetXaxis()->CenterTitle(true);
    histogram_TotalCell->GetYaxis()->CenterTitle(true);
    histogram_TotalCell->GetXaxis()->SetTitle(xAxisName.c_str());
    histogram_TotalCell->GetYaxis()->SetTitle(yAxisName.c_str());
    histogram_TotalCell->SetTitle("");

    histogram_Nucleus->GetXaxis()->CenterTitle(true);
    histogram_Nucleus->GetYaxis()->CenterTitle(true);
    histogram_Nucleus->GetXaxis()->SetTitle(xAxisName.c_str());
    histogram_Nucleus->GetYaxis()->SetTitle(yAxisName.c_str());
    histogram_Nucleus->SetTitle("");

    histogram_Membrane->GetXaxis()->CenterTitle(true);
    histogram_Membrane->GetYaxis()->CenterTitle(true);
    histogram_Membrane->GetXaxis()->SetTitle(xAxisName.c_str());
    histogram_Membrane->GetYaxis()->SetTitle(yAxisName.c_str());
    histogram_Membrane->SetTitle("");

    histogram_Cytoplasm->GetXaxis()->CenterTitle(true);
    histogram_Cytoplasm->GetYaxis()->CenterTitle(true);
    histogram_Cytoplasm->GetXaxis()->SetTitle(xAxisName.c_str());
    histogram_Cytoplasm->GetYaxis()->SetTitle(yAxisName.c_str());
    histogram_Cytoplasm->SetTitle("");


    histogram_TotalCell->GetXaxis()->SetTitleSize(titleSize);
    histogram_Nucleus->GetXaxis()->SetTitleSize(titleSize);
    histogram_Membrane->GetXaxis()->SetTitleSize(titleSize);
    histogram_Cytoplasm->GetXaxis()->SetTitleSize(titleSize);

    histogram_TotalCell->GetYaxis()->SetTitleSize(titleSize);
    histogram_Nucleus->GetYaxis()->SetTitleSize(titleSize);
    histogram_Membrane->GetYaxis()->SetTitleSize(titleSize);
    histogram_Cytoplasm->GetYaxis()->SetTitleSize(titleSize);

    //------------------------------
    // Create the canvas
    auto c1 = new TCanvas("c1", title.c_str(), 600,500);

    gStyle->SetOptStat(0);

    //------------------------------
    int reBin;
    if(cellLine=="PC3_Flu")
    {
        reBin = 2000;
    }
    else
    {
        reBin = 10;
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

    auto legend = new TLegend(0.67,0.55,0.95,0.85);
    legend->SetHeader("Cell Component",title.c_str());
    TLegendEntry *header = (TLegendEntry*)legend->GetListOfPrimitives()->First();
    header->SetTextAlign(22);
    header->SetTextSize(.05);

    legend->AddEntry(histogram_TotalCell,"Total Cell")->SetTextSize(0.05);
    legend->AddEntry(histogram_Membrane,"Membrane")->SetTextSize(0.05);
    legend->AddEntry(histogram_Cytoplasm,"Cytoplasm")->SetTextSize(0.05);
    legend->AddEntry(histogram_Nucleus,"Nucleus")->SetTextSize(0.05);
    legend->Draw();

    c1->SetTopMargin(0.15);    // Set the top margin (5% of the canvas height)
    c1->SetBottomMargin(0.15); // Set the bottom margin (15% of the canvas height)
    c1->SetLeftMargin(0.15);   // Set the left margin (15% of the canvas width)
    c1->SetRightMargin(0.05);  // Set the right margin (5% of the canvas width)

    TLatex *latex1 = new TLatex();
    latex1->SetNDC();
    latex1->SetTextSize(0.04); // Adjust to appropriate subtitle size
    latex1->DrawLatex(0.15, 0.89,title.c_str()); // Position the subtitle. Adjust x, y to fit well.


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
        if(cellLine=="PC3_Flu")
        {
            histogramName_TotalCell = "i0_hEnergyDeps_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_TotalCell_eVBinning";
            histogramName_Nucleus = "i0_hEnergyDeps_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_Nucleus_eVBinning";
            histogramName_Membrane = "i0_hEnergyDeps_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_Membrane_eVBinning";
            histogramName_Cytoplasm = "i0_hEnergyDeps_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_Cytoplasm_eVBinning";
        }
        else
        {
            histogramName_TotalCell = "i0_hEnergyDeps_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_TotalCell_keVBinning";
            histogramName_Nucleus = "i0_hEnergyDeps_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_Nucleus_keVBinning";
            histogramName_Membrane = "i0_hEnergyDeps_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_Membrane_keVBinning";
            histogramName_Cytoplasm = "i0_hEnergyDeps_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_Cytoplasm_keVBinning";
        }

    }
    if(doseOrEnergy==1)
    {
        if(cellLine=="PC3_Flu")
        {
            histogramName_TotalCell = "i0_hDose_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_TotalCell_uGyBinning";
            histogramName_Nucleus = "i0_hDose_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_Nucleus_uGyBinning";
            histogramName_Membrane = "i0_hDose_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_Membrane_uGyBinning";
            histogramName_Cytoplasm = "i0_hDose_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_Cytoplasm_uGyBinning";
        }
        else
        {
            histogramName_TotalCell = "i0_hDose_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_TotalCell_mGyBinning";
            histogramName_Nucleus = "i0_hDose_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_Nucleus_mGyBinning";
            histogramName_Membrane = "i0_hDose_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_Membrane_mGyBinning";
            histogramName_Cytoplasm = "i0_hDose_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_Cytoplasm_mGyBinning";
        }
    }

    //------------------------------------------------
    std::string fileName = "../Output_" + cellGeometry + "/Output_"+ cellLine + "_" + std::to_string(activity) + "kBq.root";
    auto inputFile = std::unique_ptr<TFile>(TFile::Open(fileName.c_str()));

    double titleSize = 0.045;

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
    std::string generalTitle;


    //------------------------------
    int reBin;
    if(cellLine=="PC3_Flu")
    {
        if(diameter==12)
        {
            reBin = 1000;
        }
        if(diameter==5)
        {
            reBin = 1000;
        }
    }
    if(cellLine=="C4_2")
    {
        if(diameter==12)
        {
            reBin = 4;
        }
        if(diameter==5)
        {
            reBin = 8;
        }
    }
    if(cellLine=="PC3_PIP")
    {
        if(diameter==12)
        {
            reBin = 20;
        }
        if(diameter==5)
        {
            reBin = 20;
        }
    }

    if(doseOrEnergy==0)
    {
        xAxisName = "Energy Deposited [MeV]";
        generalTitle = "Energy Deposition";
        if(cellLine=="PC3_Flu")
        {
            yAxisName = "Fraction of cells in sample / 1 keV bin";
        }
        if(cellLine=="C4_2")
        {
            if(diameter==12)
            {
                yAxisName = "Fraction of cells in sample / 40 keV bin";
            }
            if(diameter==5)
            {
                yAxisName = "Fraction of cells in sample / 80 keV bin";
            }
        }
        if(cellLine=="PC3_PIP")
        {
            if(diameter==12)
            {
                yAxisName = "Fraction of cells in sample / 200 keV bin";
            }
            if(diameter==5)
            {
                yAxisName = "Fraction of cells in sample / 200 keV bin";
            }
        }
    }
    if(doseOrEnergy==1)
    {
        xAxisName = "Dose Delivered [Gy]";
        generalTitle = "Absorbed Dose";
        if(cellLine=="PC3_Flu")
        {
            yAxisName = "Fraction of cells in sample / 1 mGy bin";
        }
        else
        {
            yAxisName = "Fraction of cells in sample / 4 mGy bin";
        }
    }

    std::string title = cellLine_Name + ", " + nucleiDist + " Dist. Nuclei, d_{nuc} = " + Form("%d", diameter) + " #mum, " + Form("%d", activity)   + " kBq / mL";

    //------------------------
    histogram_TotalCell->GetXaxis()->CenterTitle(true);
    histogram_TotalCell->GetYaxis()->CenterTitle(true);
    histogram_TotalCell->GetXaxis()->SetTitle(xAxisName.c_str());
    histogram_TotalCell->GetYaxis()->SetTitle(yAxisName.c_str());
    histogram_TotalCell->SetTitle("");

    histogram_Nucleus->GetXaxis()->CenterTitle(true);
    histogram_Nucleus->GetYaxis()->CenterTitle(true);
    histogram_Nucleus->GetXaxis()->SetTitle(xAxisName.c_str());
    histogram_Nucleus->GetYaxis()->SetTitle(yAxisName.c_str());
    histogram_Nucleus->SetTitle("");

    histogram_Membrane->GetXaxis()->CenterTitle(true);
    histogram_Membrane->GetYaxis()->CenterTitle(true);
    histogram_Membrane->GetXaxis()->SetTitle(xAxisName.c_str());
    histogram_Membrane->GetYaxis()->SetTitle(yAxisName.c_str());
    histogram_Membrane->SetTitle("");

    histogram_Cytoplasm->GetXaxis()->CenterTitle(true);
    histogram_Cytoplasm->GetYaxis()->CenterTitle(true);
    histogram_Cytoplasm->GetXaxis()->SetTitle(xAxisName.c_str());
    histogram_Cytoplasm->GetYaxis()->SetTitle(yAxisName.c_str());
    histogram_Cytoplasm->SetTitle("");


    histogram_TotalCell->GetXaxis()->SetTitleSize(titleSize);
    histogram_Nucleus->GetXaxis()->SetTitleSize(titleSize);
    histogram_Membrane->GetXaxis()->SetTitleSize(titleSize);
    histogram_Cytoplasm->GetXaxis()->SetTitleSize(titleSize);

    histogram_TotalCell->GetYaxis()->SetTitleSize(titleSize);
    histogram_Nucleus->GetYaxis()->SetTitleSize(titleSize);
    histogram_Membrane->GetYaxis()->SetTitleSize(titleSize);
    histogram_Cytoplasm->GetYaxis()->SetTitleSize(titleSize);

    //------------------------------
    // Create the canvas
    auto c1 = new TCanvas("c1", title.c_str(), 600,500);

    gStyle->SetOptStat(0);

    histogram_TotalCell->Rebin(reBin);
    histogram_Cytoplasm->Rebin(reBin);
    histogram_Nucleus->Rebin(reBin);
    histogram_Membrane->Rebin(reBin);

    std::cout << histogram_Nucleus->GetBinWidth(0) << std::endl;


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
        xMax = 1./6.*xMax;
    }



    histogram_TotalCell->GetXaxis()->SetRangeUser(0.,xMax);
    histogram_TotalCell->GetYaxis()->SetRangeUser(1.e-6,10.);

    histogram_Cytoplasm->GetXaxis()->SetRangeUser(0.,xMax);
    histogram_Cytoplasm->GetYaxis()->SetRangeUser(1.e-6,10.);

    histogram_Nucleus->GetXaxis()->SetRangeUser(0.,xMax);
    histogram_Nucleus->GetYaxis()->SetRangeUser(1.e-6,10.);

    histogram_Membrane->GetXaxis()->SetRangeUser(0.,xMax);
    histogram_Membrane->GetYaxis()->SetRangeUser(1.e-6,10.);


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

    auto legend = new TLegend(0.67,0.55,0.95,0.85);
    legend->SetHeader("Cell Component",title.c_str());
    TLegendEntry *header = (TLegendEntry*)legend->GetListOfPrimitives()->First();
    header->SetTextAlign(22);
    header->SetTextSize(.04);

    legend->AddEntry(histogram_TotalCell,"Total Cell")->SetTextSize(0.04);
    legend->AddEntry(histogram_Membrane,"Membrane")->SetTextSize(0.04);
    legend->AddEntry(histogram_Cytoplasm,"Cytoplasm")->SetTextSize(0.04);
    legend->AddEntry(histogram_Nucleus,"Nucleus")->SetTextSize(0.04);
    legend->Draw();

    c1->SetTopMargin(0.15);    // Set the top margin (5% of the canvas height)
    c1->SetBottomMargin(0.15); // Set the bottom margin (15% of the canvas height)
    c1->SetLeftMargin(0.15);   // Set the left margin (15% of the canvas width)
    c1->SetRightMargin(0.05);  // Set the right margin (5% of the canvas width)

    TLatex *latex1 = new TLatex();
    latex1->SetNDC();
    latex1->SetTextSize(0.04); // Adjust to appropriate subtitle size
    latex1->DrawLatex(0.15, 0.89,title.c_str()); // Position the subtitle. Adjust x, y to fit well.


    std::string figureName;
    if(doseOrEnergy==0)
    {
        figureName =  "/Users/idarosenqvist/Desktop/Academics/MasterThesis/Thesis/figures/Results/" + cellGeometry + "/" + cellLine + "/EnergyDeposition_" + cellLine + "_" + std::to_string(activity) + "kBq_logy.pdf";
    }
    if(doseOrEnergy==1)
    {
        figureName = "/Users/idarosenqvist/Desktop/Academics/MasterThesis/Thesis/figures/Results/" + cellGeometry + "/" + cellLine + "/DoseDelivered_" + cellLine + "_" + std::to_string(activity) + "kBq_logy.pdf";
    }

    c1->Update();
    c1->SetLogy();
    c1->Update();
    c1->SaveAs(figureName.c_str());

}

void PlotDoseByCellComponent()
{
    std::string cellLine = "C4_2";
    // std::string cellLine = "PC3_PIP";
    // std::string cellLine = "PC3_Flu";

    MakePlots(cellLine, "D12RP", 25, 1);
    MakePlotsLog(cellLine, "D12RP", 25, 1);

    MakePlots(cellLine, "D12RP", 25, 0);
    MakePlotsLog(cellLine, "D12RP", 25, 0);


    MakePlots(cellLine, "D12CP", 25, 1);
    MakePlotsLog(cellLine, "D12CP", 25, 1);

    MakePlots(cellLine, "D12CP", 25, 0);
    MakePlotsLog(cellLine, "D12CP", 25, 0);

    MakePlots(cellLine, "D5RP", 25, 1);
    MakePlotsLog(cellLine, "D5RP", 25, 1);

    MakePlots(cellLine, "D5RP", 25, 0);
    MakePlotsLog(cellLine, "D5RP", 25, 0);

    MakePlots(cellLine, "D5CP", 25, 1);
    MakePlotsLog(cellLine, "D5CP", 25, 1);

    MakePlots(cellLine, "D5CP", 25, 0);
    MakePlotsLog(cellLine, "D5CP", 25, 0);

}