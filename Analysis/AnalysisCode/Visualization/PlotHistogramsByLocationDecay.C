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

void PlotHistogramsLocationDecay(std::string cellLine, std::string cellGeometry, std::string cellComponent, int activity)
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

    //------------------------------------------------
    std::string fileName = "../Output_" + cellGeometry + "/Output_"+ cellLine + "_" + std::to_string(activity) + "kBq.root";
    auto inputFile = std::unique_ptr<TFile>(TFile::Open(fileName.c_str()));



    if(cellLine=="PC3_Flu")
    {
        //------------------------------
        std::string histogramName_Dose = "i0_hDose_212Pb_"+ cellLine + "_" + std::to_string(activity) + "kBq_" + cellComponent + "_uGyBinning";
        std::string histogramName_Dose_FractionSolution = "i0_hDose_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_" + cellComponent + "_FractionFromSolution";

        //-----------------------------
        TH1D* histogramDose = 0;
        inputFile->GetObject(histogramName_Dose.c_str(), histogramDose);
        histogramDose->SetDirectory(0);

        //-----------------------------
        TH1D* histogramDose_FractionFromSolution = 0;
        inputFile->GetObject(histogramName_Dose_FractionSolution.c_str(), histogramDose_FractionFromSolution);
        histogramDose_FractionFromSolution->SetDirectory(0);

        gStyle->SetOptStat(0);

        //------------------------------
        std::string region;
        if(cellComponent=="Nucleus"){region = "nucleus";}
        if(cellComponent=="Cytoplasm"){region = "cytoplasm";}
        if(cellComponent=="Membrane"){region = "membrane";}
        if(cellComponent=="TotalCell"){region = "cell";}

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


        //------------------------------
        // Create the canvas
        std::string canvasName = "c1_" + cellLine + "_" + std::to_string(activity) + "kBq";
        auto c1 = new TCanvas(canvasName.c_str(), "Dose Distribution Split by Location of 212-Pb Decay", 600, 400);


        histogramDose_FractionFromSolution->SetLineColor(colours[0]);
        histogramDose_FractionFromSolution->SetFillColor(colours[0]);

        histogramDose->SetLineColor(kBlack);
        histogramDose->SetLineWidth(2);

        //--------------------------------

        double doseMin;
        double doseMax;

        auto [firstBin, lastBin] = FindNonEmptyBinRange(histogramDose);
        if (firstBin <= lastBin) {
            doseMin = histogramDose->GetBinLowEdge(firstBin);
            doseMax = histogramDose->GetBinLowEdge(lastBin) + histogramDose->GetBinWidth(lastBin);

        } else {
            std::cerr << "The histogram is empty!" << std::endl;
        }


        int reBin;
        if(cellComponent=="Nucleus")
        {
            if(diameter==5)
            {
                doseMax = doseMax*1./3.;
                reBin = 1000;
            }
            if(diameter==12)
            {
                doseMax = doseMax*1./3.;
                reBin = 100;
            }
        }
        if(cellComponent=="TotalCell")
        {
            if(diameter==5)
            {
                doseMax = doseMax*1./3.;
                reBin = 1000;
            }
            if(diameter==12)
            {
                doseMax = doseMax*1./3.;
                reBin = 100;
            }
        }

        // //------------------------------

        histogramDose_FractionFromSolution->Rebin(reBin);
        histogramDose->Rebin(reBin);

        double maxY = histogramDose->GetMaximum() + 0.05*histogramDose->GetMaximum();
        double minY = 1.e-7;

        histogramDose_FractionFromSolution->GetXaxis()->SetRangeUser(doseMin,doseMax);
        histogramDose->GetXaxis()->SetRangeUser(doseMin,doseMax);

        histogramDose_FractionFromSolution->GetYaxis()->SetRangeUser(minY,maxY);
        histogramDose->GetYaxis()->SetRangeUser(minY,maxY);


        std::string xAxisName = "Dose delivered in " + region + " [Gy]";
        std::string yAxisName = "Fraction of cells in sample / " + std::to_string(reBin) + " mGy bin";
        std::string title = cellLine_Name + ", " + nucleiDist + " Dist. Nuclei, d_{nuc} = " + Form("%d", diameter) + " #mum, Activity = " + Form("%d", activity)   + "kBq / 1mL";

        histogramDose_FractionFromSolution->GetXaxis()->CenterTitle(true);
        histogramDose_FractionFromSolution->GetYaxis()->CenterTitle(true);
        histogramDose_FractionFromSolution->GetXaxis()->SetTitle(xAxisName.c_str());
        histogramDose_FractionFromSolution->GetYaxis()->SetTitle(yAxisName.c_str());
        histogramDose_FractionFromSolution->SetTitle(title.c_str());

        histogramDose->GetXaxis()->CenterTitle(true);
        histogramDose->GetYaxis()->CenterTitle(true);
        histogramDose->GetXaxis()->SetTitle(xAxisName.c_str());
        histogramDose->GetYaxis()->SetTitle(yAxisName.c_str());
        histogramDose->SetTitle(title.c_str());

        histogramDose_FractionFromSolution->Draw("HIST");
        histogramDose->Draw("HIST SAME");

        auto legend = new TLegend(0.60,0.75,0.9,0.9);
        legend->SetTextSize(0.035);

        legend->SetHeader("Location of ^{212}Pb Decay","C"); // option "C" allows to center the header
        TLegendEntry *header = (TLegendEntry*)legend->GetListOfPrimitives()->First();
        header->SetTextAlign(22);
        header->SetTextSize(.04);

        legend->AddEntry(histogramDose,"From all locations");
        legend->AddEntry(histogramDose_FractionFromSolution,"Solution");
        legend->Draw();

        c1->Update();

        std::string output = "/Users/idarosenqvist/Desktop/Academics/MasterThesis/Thesis/figures/Results/" + cellGeometry + "/" + cellLine + "/LocationDecays_" + cellLine + "_" + cellComponent + "_" + std::to_string(activity) + "kBq.pdf";
        c1->SaveAs(output.c_str());


        // //--------------------------------------
        // histogramDose_FractionFromSolution->GetYaxis()->SetRangeUser(minY,1.);
        // histogramDose->GetYaxis()->SetRangeUser(minY,1.);

        // c1->SetLogy();
        // c1->Update();

        // output = "/Users/idarosenqvist/Desktop/Academics/MasterThesis/Thesis/figures/Results/" + cellGeometry + "/" + cellLine + "/LocationDecays_" + cellLine + "_" + cellComponent + "_" + std::to_string(activity) + "kBq_logy.pdf";
        // c1->SaveAs(output.c_str());
    }

    if(cellLine=="C4_2"||cellLine=="PC3_PIP")
    {
        //------------------------------
        std::string histogramName_Dose = "i0_hDose_212Pb_"+ cellLine + "_" + std::to_string(activity) + "kBq_" + cellComponent + "_uGyBinning";
        std::string histogramName_Dose_FractionSolution = "i0_hDose_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_" + cellComponent + "_FractionFromSolution";
        std::string histogramName_Dose_FractionMembrane = "i0_hDose_212Pb_"+ cellLine + "_" + std::to_string(activity) + "kBq_" + cellComponent + "_FractionFromMembrane";
        std::string histogramName_Dose_FractionCytoplasm = "i0_hDose_212Pb_"+ cellLine + "_" + std::to_string(activity) + "kBq_" + cellComponent + "_FractionFromCytoplasm";

        std::string histogramName_Dose_SolutionPlusMembrane = "i0_hDose_212Pb_"+ cellLine + "_"+ std::to_string(activity) + "kBq_" + cellComponent + "_FractionFromSolutionPlusMembrane";
        std::string histogramName_Dose_SolutionPlusMembranePlusCytoplasm = "i0_hDose_212Pb_"+ cellLine + "_"+ std::to_string(activity) + "kBq_" + cellComponent + "_FractionFromSolutionPlusMembranePlusCytoplasm";


        std::string region;
        if(cellComponent=="Nucleus"){region = "nucleus";}
        if(cellComponent=="Cytoplasm"){region = "cytoplasm";}
        if(cellComponent=="Membrane"){region = "membrane";}
        if(cellComponent=="TotalCell"){region = "cell";}

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


        //-----------------------------
        TH1D* histogramDose = 0;
        inputFile->GetObject(histogramName_Dose.c_str(), histogramDose);
        histogramDose->SetDirectory(0);

        TH1D* histogramDose_FractionFromSolution = 0;
        inputFile->GetObject(histogramName_Dose_FractionSolution.c_str(), histogramDose_FractionFromSolution);
        histogramDose_FractionFromSolution->SetDirectory(0);

        TH1D* histogramDose_FractionFromMembrane = 0;
        inputFile->GetObject(histogramName_Dose_FractionMembrane.c_str(), histogramDose_FractionFromMembrane);
        histogramDose_FractionFromMembrane->SetDirectory(0);

        TH1D* histogramDose_FractionFromCytoplasm = 0;
        inputFile->GetObject(histogramName_Dose_FractionCytoplasm.c_str(), histogramDose_FractionFromCytoplasm);
        histogramDose_FractionFromCytoplasm->SetDirectory(0);


        //-----------------------------
        TH1D* sumHist_SolutionPlusMembrane = 0;
        inputFile->GetObject(histogramName_Dose_FractionSolution.c_str(), sumHist_SolutionPlusMembrane);
        sumHist_SolutionPlusMembrane->SetDirectory(0);
        sumHist_SolutionPlusMembrane->Reset();

        sumHist_SolutionPlusMembrane->Add(histogramDose_FractionFromMembrane);
        sumHist_SolutionPlusMembrane->Add(histogramDose_FractionFromSolution);

        // //--------------------------------
        TH1D* sumHist_SolutionPlusMembranePlusCytoplasm = 0;
        inputFile->GetObject(histogramName_Dose_FractionSolution.c_str(), sumHist_SolutionPlusMembranePlusCytoplasm);
        sumHist_SolutionPlusMembranePlusCytoplasm->SetDirectory(0);
        sumHist_SolutionPlusMembranePlusCytoplasm->Reset();


        sumHist_SolutionPlusMembranePlusCytoplasm->Add(histogramDose_FractionFromSolution);
        sumHist_SolutionPlusMembranePlusCytoplasm->Add(histogramDose_FractionFromMembrane);
        sumHist_SolutionPlusMembranePlusCytoplasm->Add(histogramDose_FractionFromCytoplasm);


        gStyle->SetOptStat(0);

        //------------------------------
        // Create the canvas
        std::string canvasName = "c1_" + cellLine + "_" + std::to_string(activity) + "kBq";
        auto c1 = new TCanvas(canvasName.c_str(), "Dose Distribution Split by Location of 212-Pb Decay", 600, 400);


        histogramDose_FractionFromSolution->SetLineColor(colours[0]);
        histogramDose_FractionFromSolution->SetFillColor(colours[0]);

        sumHist_SolutionPlusMembrane->SetLineColor(colours[1]);
        sumHist_SolutionPlusMembrane->SetFillColor(colours[1]);

        sumHist_SolutionPlusMembranePlusCytoplasm->SetLineColor(colours[2]);
        sumHist_SolutionPlusMembranePlusCytoplasm->SetFillColor(colours[2]);

        histogramDose->SetLineColor(kBlack);
        histogramDose->SetLineWidth(2);

        double doseMin;
        double doseMax;

        auto [firstBin, lastBin] = FindNonEmptyBinRange(histogramDose);
        if (firstBin <= lastBin) {
            doseMin = histogramDose->GetBinLowEdge(firstBin);
            doseMax = histogramDose->GetBinLowEdge(lastBin) + histogramDose->GetBinWidth(lastBin);

        } else {
            std::cerr << "The histogram is empty!" << std::endl;
        }

        int reBin;
        if(cellLine=="C4_2")
        {
            if(diameter==5 && cellComponent=="Nucleus")
            {
                doseMax = doseMax*1./6.;
                reBin = 500;
            }
            if(diameter==12 && cellComponent=="Nucleus")
            {
                doseMax = doseMax*1./6.;
                reBin = 500;
            }
            if(diameter==12 && cellComponent=="TotalCell")
            {
                reBin = 1000;
            }
            if(diameter==5 && cellComponent=="TotalCell")
            {
                reBin = 1000;
            }
        }
        if(cellLine=="PC3_PIP")
        {
            if(cellComponent=="Nucleus")
            {
                reBin = 4000;
            }
            if(cellComponent=="TotalCell")
            {
                reBin = 4000;
            }
        }

        histogramDose_FractionFromSolution->Rebin(reBin);
        sumHist_SolutionPlusMembrane->Rebin(reBin);
        sumHist_SolutionPlusMembranePlusCytoplasm->Rebin(reBin);
        histogramDose->Rebin(reBin);

        double maxY = histogramDose->GetMaximum() + 0.05*histogramDose->GetMaximum();
        double minY = 1.e-7;


        histogramDose_FractionFromSolution->GetXaxis()->SetRangeUser(doseMin,doseMax);
        sumHist_SolutionPlusMembrane->GetXaxis()->SetRangeUser(doseMin,doseMax);
        sumHist_SolutionPlusMembranePlusCytoplasm->GetXaxis()->SetRangeUser(doseMin,doseMax);
        histogramDose->GetXaxis()->SetRangeUser(doseMin,doseMax);

        histogramDose_FractionFromSolution->GetYaxis()->SetRangeUser(minY,maxY);
        sumHist_SolutionPlusMembrane->GetYaxis()->SetRangeUser(minY,maxY);
        sumHist_SolutionPlusMembranePlusCytoplasm->GetYaxis()->SetRangeUser(minY,maxY);
        histogramDose->GetYaxis()->SetRangeUser(minY,maxY);



        std::string xAxisName = "Dose delivered in " + region + " [Gy]";
        double unitBinGray = ((double)reBin)*1.e-4;
        std::string yAxisName = "Fraction of cells in sample / " + std::string(Form("%.1f", unitBinGray)) + " mGy bin";
        std::string title = cellLine_Name + ", " + nucleiDist + " Dist. Nuclei, d_{nuc} = " + Form("%d", diameter) + " #mum, Activity = " + Form("%d", activity)   + "kBq / 1mL";

        histogramDose_FractionFromSolution->GetXaxis()->CenterTitle(true);
        histogramDose_FractionFromSolution->GetYaxis()->CenterTitle(true);
        histogramDose_FractionFromSolution->GetXaxis()->SetTitle(xAxisName.c_str());
        histogramDose_FractionFromSolution->GetYaxis()->SetTitle(yAxisName.c_str());
        histogramDose_FractionFromSolution->SetTitle(title.c_str());

        sumHist_SolutionPlusMembrane->GetXaxis()->CenterTitle(true);
        sumHist_SolutionPlusMembrane->GetYaxis()->CenterTitle(true);
        sumHist_SolutionPlusMembrane->GetXaxis()->SetTitle(xAxisName.c_str());
        sumHist_SolutionPlusMembrane->GetYaxis()->SetTitle(yAxisName.c_str());
        sumHist_SolutionPlusMembrane->SetTitle(title.c_str());

        sumHist_SolutionPlusMembranePlusCytoplasm->GetXaxis()->CenterTitle(true);
        sumHist_SolutionPlusMembranePlusCytoplasm->GetYaxis()->CenterTitle(true);
        sumHist_SolutionPlusMembranePlusCytoplasm->GetXaxis()->SetTitle(xAxisName.c_str());
        sumHist_SolutionPlusMembranePlusCytoplasm->GetYaxis()->SetTitle(yAxisName.c_str());
        sumHist_SolutionPlusMembranePlusCytoplasm->SetTitle(title.c_str());

        histogramDose->GetXaxis()->CenterTitle(true);
        histogramDose->GetYaxis()->CenterTitle(true);
        histogramDose->GetXaxis()->SetTitle(xAxisName.c_str());
        histogramDose->GetYaxis()->SetTitle(yAxisName.c_str());
        histogramDose->SetTitle(title.c_str());


        sumHist_SolutionPlusMembranePlusCytoplasm->Draw("HIST");
        sumHist_SolutionPlusMembrane->Draw("HIST SAME");
        histogramDose_FractionFromSolution->Draw("HIST SAME");
        histogramDose->Draw("HIST SAME");

        // sumHist_SolutionPlusMembrane->Draw("HIST");

        auto legend = new TLegend(0.67,0.7,0.9,0.9);

        legend->SetHeader("Location of ^{212}Pb Decay","C"); // option "C" allows to center the header
        TLegendEntry *header = (TLegendEntry*)legend->GetListOfPrimitives()->First();
        header->SetTextAlign(22);
        header->SetTextSize(.032);

        legend->AddEntry(histogramDose,"From all locations")->SetTextSize(0.03);
        legend->AddEntry(histogramDose_FractionFromSolution,"Solution")->SetTextSize(0.03);
        legend->AddEntry(sumHist_SolutionPlusMembrane,"Membrane")->SetTextSize(0.03);
        legend->AddEntry(sumHist_SolutionPlusMembranePlusCytoplasm,"Cytoplasm")->SetTextSize(0.03);
        legend->Draw();

        c1->Update();

        std::string output = "/Users/idarosenqvist/Desktop/Academics/MasterThesis/Thesis/figures/Results/" + cellGeometry + "/" + cellLine + "/LocationDecays_" + cellLine + "_" + cellComponent + "_" + std::to_string(activity) + "kBq.pdf";
        c1->SaveAs(output.c_str());



        // histogramDose_FractionFromSolution->GetYaxis()->SetRangeUser(minY,1.);
        // sumHist_SolutionPlusMembrane->GetYaxis()->SetRangeUser(minY,1.);
        // sumHist_SolutionPlusMembranePlusCytoplasm->GetYaxis()->SetRangeUser(minY,1.);
        // histogramDose->GetYaxis()->SetRangeUser(minY,1.);

        // c1->SetLogy();
        // c1->Update();

        // output = "/Users/idarosenqvist/Desktop/Academics/MasterThesis/Thesis/figures/Results/" + cellGeometry + "/" + cellLine + "/LocationDecays_" + cellLine + "_" + cellComponent + "_" + std::to_string(activity) + "kBq_logy.pdf";
        // c1->SaveAs(output.c_str());
    }


}

void PlotHistogramsLocationDecayLog(std::string cellLine, std::string cellGeometry, std::string cellComponent, int activity)
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

    //------------------------------------------------
    std::string fileName = "../Output_" + cellGeometry + "/Output_"+ cellLine + "_" + std::to_string(activity) + "kBq.root";
    auto inputFile = std::unique_ptr<TFile>(TFile::Open(fileName.c_str()));

    if(cellLine=="PC3_Flu")
    {
        //------------------------------
        std::string histogramName_Dose = "i0_hDose_212Pb_"+ cellLine + "_" + std::to_string(activity) + "kBq_" + cellComponent + "_mGyBinning";
        std::string histogramName_Dose_FractionSolution = "i0_hDose_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_" + cellComponent + "_FractionFromSolution";

        //-----------------------------
        TH1D* histogramDose = 0;
        inputFile->GetObject(histogramName_Dose.c_str(), histogramDose);
        histogramDose->SetDirectory(0);

        //-----------------------------
        TH1D* histogramDose_FractionFromSolution = 0;
        inputFile->GetObject(histogramName_Dose_FractionSolution.c_str(), histogramDose_FractionFromSolution);
        histogramDose_FractionFromSolution->SetDirectory(0);

        gStyle->SetOptStat(0);

        //------------------------------
        // Create the canvas
        std::string canvasName = "c1_" + cellLine + "_" + std::to_string(activity) + "kBq";
        auto c1 = new TCanvas(canvasName.c_str(), "Dose Distribution Split by Location of 212-Pb Decay", 600, 400);


        histogramDose_FractionFromSolution->SetLineColor(colours[0]);
        histogramDose_FractionFromSolution->SetFillColor(colours[0]);

        histogramDose->SetLineColor(kBlack);
        histogramDose->SetLineWidth(2);

        //--------------------------------

        double doseMin;
        double doseMax;

        auto [firstBin, lastBin] = FindNonEmptyBinRange(histogramDose);
        if (firstBin <= lastBin) {
            doseMin = histogramDose->GetBinLowEdge(firstBin);
            doseMax = histogramDose->GetBinLowEdge(lastBin) + histogramDose->GetBinWidth(lastBin);

        } else {
            std::cerr << "The histogram is empty!" << std::endl;
        }

        double maxY = histogramDose->GetMaximum() + 0.05*histogramDose->GetMaximum();
        double minY = 1.e-7;

        // //------------------------------
        int reBin = 1;

        // histogramDose_FractionFromSolution->Rebin(reBin);
        // histogramDose->Rebin(reBin);

        histogramDose_FractionFromSolution->GetXaxis()->SetRangeUser(doseMin,doseMax);
        histogramDose->GetXaxis()->SetRangeUser(doseMin,doseMax);

        histogramDose_FractionFromSolution->GetYaxis()->SetRangeUser(minY,maxY);
        histogramDose->GetYaxis()->SetRangeUser(minY,maxY);

        //------------------------------
        std::string region;
        if(cellComponent=="Nucleus"){region = "nucleus";}
        if(cellComponent=="Cytoplasm"){region = "cytoplasm";}
        if(cellComponent=="Membrane"){region = "membrane";}
        if(cellComponent=="TotalCell"){region = "cell";}

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

        std::string xAxisName = "Dose delivered in " + region + " [Gy]";
        std::string yAxisName = "Fraction of cells in sample / " + std::to_string(reBin) + " mGy bin";
        std::string title = cellLine_Name + ", " + nucleiDist + " Dist. Nuclei, d_{nuc} = " + Form("%d", diameter) + " #mum, Activity = " + Form("%d", activity)   + "kBq / 1mL";

        histogramDose_FractionFromSolution->GetXaxis()->CenterTitle(true);
        histogramDose_FractionFromSolution->GetYaxis()->CenterTitle(true);
        histogramDose_FractionFromSolution->GetXaxis()->SetTitle(xAxisName.c_str());
        histogramDose_FractionFromSolution->GetYaxis()->SetTitle(yAxisName.c_str());
        histogramDose_FractionFromSolution->SetTitle(title.c_str());

        histogramDose->GetXaxis()->CenterTitle(true);
        histogramDose->GetYaxis()->CenterTitle(true);
        histogramDose->GetXaxis()->SetTitle(xAxisName.c_str());
        histogramDose->GetYaxis()->SetTitle(yAxisName.c_str());
        histogramDose->SetTitle(title.c_str());

        histogramDose_FractionFromSolution->Draw("HIST");
        histogramDose->Draw("HIST SAME");

        auto legend = new TLegend(0.60,0.75,0.9,0.9);
        legend->SetTextSize(0.035);

        legend->SetHeader("Location of ^{212}Pb Decay","C"); // option "C" allows to center the header
        TLegendEntry *header = (TLegendEntry*)legend->GetListOfPrimitives()->First();
        header->SetTextAlign(22);
        header->SetTextSize(.04);

        legend->AddEntry(histogramDose,"From all locations");
        legend->AddEntry(histogramDose_FractionFromSolution,"Solution");
        legend->Draw();


        //--------------------------------------
        histogramDose_FractionFromSolution->GetYaxis()->SetRangeUser(minY,1.);
        histogramDose->GetYaxis()->SetRangeUser(minY,1.);

        c1->SetLogy();
        c1->Update();

        std::string output = "/Users/idarosenqvist/Desktop/Academics/MasterThesis/Thesis/figures/Results/" + cellGeometry + "/" + cellLine + "/LocationDecays_" + cellLine + "_" + cellComponent + "_" + std::to_string(activity) + "kBq_logy.pdf";
        c1->SaveAs(output.c_str());
    }

    if(cellLine=="C4_2"||cellLine=="PC3_PIP")
    {
        //------------------------------
        std::string histogramName_Dose = "i0_hDose_212Pb_"+ cellLine + "_" + std::to_string(activity) + "kBq_" + cellComponent + "_mGyBinning";
        std::string histogramName_Dose_FractionSolution = "i0_hDose_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_" + cellComponent + "_FractionFromSolution";
        std::string histogramName_Dose_FractionMembrane = "i0_hDose_212Pb_"+ cellLine + "_" + std::to_string(activity) + "kBq_" + cellComponent + "_FractionFromMembrane";
        std::string histogramName_Dose_FractionCytoplasm = "i0_hDose_212Pb_"+ cellLine + "_" + std::to_string(activity) + "kBq_" + cellComponent + "_FractionFromCytoplasm";

        std::string histogramName_Dose_SolutionPlusMembrane = "i0_hDose_212Pb_"+ cellLine + "_"+ std::to_string(activity) + "kBq_" + cellComponent + "_FractionFromSolutionPlusMembrane";
        std::string histogramName_Dose_SolutionPlusMembranePlusCytoplasm = "i0_hDose_212Pb_"+ cellLine + "_"+ std::to_string(activity) + "kBq_" + cellComponent + "_FractionFromSolutionPlusMembranePlusCytoplasm";

        //-----------------------------
        TH1D* histogramDose = 0;
        inputFile->GetObject(histogramName_Dose.c_str(), histogramDose);
        histogramDose->SetDirectory(0);

        TH1D* histogramDose_FractionFromSolution = 0;
        inputFile->GetObject(histogramName_Dose_FractionSolution.c_str(), histogramDose_FractionFromSolution);
        histogramDose_FractionFromSolution->SetDirectory(0);

        TH1D* histogramDose_FractionFromMembrane = 0;
        inputFile->GetObject(histogramName_Dose_FractionMembrane.c_str(), histogramDose_FractionFromMembrane);
        histogramDose_FractionFromMembrane->SetDirectory(0);

        TH1D* histogramDose_FractionFromCytoplasm = 0;
        inputFile->GetObject(histogramName_Dose_FractionCytoplasm.c_str(), histogramDose_FractionFromCytoplasm);
        histogramDose_FractionFromCytoplasm->SetDirectory(0);



        //-----------------------------
        TH1D* sumHist_SolutionPlusMembrane = 0;
        inputFile->GetObject(histogramName_Dose.c_str(), sumHist_SolutionPlusMembrane);
        sumHist_SolutionPlusMembrane->Reset();
        sumHist_SolutionPlusMembrane->SetName(histogramName_Dose_SolutionPlusMembrane.c_str());


        sumHist_SolutionPlusMembrane->Add(histogramDose_FractionFromSolution);
        sumHist_SolutionPlusMembrane->Add(histogramDose_FractionFromMembrane);


        //--------------------------------
        TH1D* sumHist_SolutionPlusMembranePlusCytoplasm = 0;
        inputFile->GetObject(histogramName_Dose.c_str(), sumHist_SolutionPlusMembranePlusCytoplasm);
        sumHist_SolutionPlusMembranePlusCytoplasm->Reset();
        sumHist_SolutionPlusMembranePlusCytoplasm->SetName(histogramName_Dose_SolutionPlusMembranePlusCytoplasm.c_str());

        sumHist_SolutionPlusMembranePlusCytoplasm->Add(histogramDose_FractionFromSolution);
        sumHist_SolutionPlusMembranePlusCytoplasm->Add(histogramDose_FractionFromMembrane);
        sumHist_SolutionPlusMembranePlusCytoplasm->Add(histogramDose_FractionFromCytoplasm);



        sumHist_SolutionPlusMembrane->SetDirectory(0);
        sumHist_SolutionPlusMembranePlusCytoplasm->SetDirectory(0);


        gStyle->SetOptStat(0);

        //------------------------------
        // Create the canvas
        std::string canvasName = "c1_" + cellLine + "_" + std::to_string(activity) + "kBq";
        auto c1 = new TCanvas(canvasName.c_str(), "Dose Distribution Split by Location of 212-Pb Decay", 600, 400);


        histogramDose_FractionFromSolution->SetLineColor(colours[0]);
        histogramDose_FractionFromSolution->SetFillColor(colours[0]);

        sumHist_SolutionPlusMembrane->SetLineColor(colours[1]);
        sumHist_SolutionPlusMembrane->SetFillColor(colours[1]);

        sumHist_SolutionPlusMembranePlusCytoplasm->SetLineColor(colours[2]);
        sumHist_SolutionPlusMembranePlusCytoplasm->SetFillColor(colours[2]);

        histogramDose->SetLineColor(kBlack);
        histogramDose->SetLineWidth(2);

        double doseMin;
        double doseMax;

        auto [firstBin, lastBin] = FindNonEmptyBinRange(histogramDose);
        if (firstBin <= lastBin) {
            doseMin = histogramDose->GetBinLowEdge(firstBin);
            doseMax = histogramDose->GetBinLowEdge(lastBin) + histogramDose->GetBinWidth(lastBin);

        } else {
            std::cerr << "The histogram is empty!" << std::endl;
        }

        double maxY = histogramDose->GetMaximum() + 0.05*histogramDose->GetMaximum();
        double minY = 1.e-7;

        int reBin = 1;


        histogramDose_FractionFromSolution->GetXaxis()->SetRangeUser(doseMin,doseMax);
        sumHist_SolutionPlusMembrane->GetXaxis()->SetRangeUser(doseMin,doseMax);
        sumHist_SolutionPlusMembranePlusCytoplasm->GetXaxis()->SetRangeUser(doseMin,doseMax);
        histogramDose->GetXaxis()->SetRangeUser(doseMin,doseMax);

        histogramDose_FractionFromSolution->GetYaxis()->SetRangeUser(minY,maxY);
        sumHist_SolutionPlusMembrane->GetYaxis()->SetRangeUser(minY,maxY);
        sumHist_SolutionPlusMembranePlusCytoplasm->GetYaxis()->SetRangeUser(minY,maxY);
        histogramDose->GetYaxis()->SetRangeUser(minY,maxY);

        std::string region;
        if(cellComponent=="Nucleus"){region = "nucleus";}
        if(cellComponent=="Cytoplasm"){region = "cytoplasm";}
        if(cellComponent=="Membrane"){region = "membrane";}
        if(cellComponent=="TotalCell"){region = "cell";}

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

        std::string xAxisName = "Dose delivered in " + region + " [Gy]";
        std::string yAxisName = "Fraction of cells in sample / " + std::to_string(reBin) + " mGy bin";
        std::string title = cellLine_Name + ", " + nucleiDist + " Dist. Nuclei, d_{nuc} = " + Form("%d", diameter) + " #mum, Activity = " + Form("%d", activity)   + "kBq / 1mL";

        histogramDose_FractionFromSolution->GetXaxis()->CenterTitle(true);
        histogramDose_FractionFromSolution->GetYaxis()->CenterTitle(true);
        histogramDose_FractionFromSolution->GetXaxis()->SetTitle(xAxisName.c_str());
        histogramDose_FractionFromSolution->GetYaxis()->SetTitle(yAxisName.c_str());
        histogramDose_FractionFromSolution->SetTitle(title.c_str());

        sumHist_SolutionPlusMembrane->GetXaxis()->CenterTitle(true);
        sumHist_SolutionPlusMembrane->GetYaxis()->CenterTitle(true);
        sumHist_SolutionPlusMembrane->GetXaxis()->SetTitle(xAxisName.c_str());
        sumHist_SolutionPlusMembrane->GetYaxis()->SetTitle(yAxisName.c_str());
        sumHist_SolutionPlusMembrane->SetTitle(title.c_str());

        sumHist_SolutionPlusMembranePlusCytoplasm->GetXaxis()->CenterTitle(true);
        sumHist_SolutionPlusMembranePlusCytoplasm->GetYaxis()->CenterTitle(true);
        sumHist_SolutionPlusMembranePlusCytoplasm->GetXaxis()->SetTitle(xAxisName.c_str());
        sumHist_SolutionPlusMembranePlusCytoplasm->GetYaxis()->SetTitle(yAxisName.c_str());
        sumHist_SolutionPlusMembranePlusCytoplasm->SetTitle(title.c_str());

        histogramDose->GetXaxis()->CenterTitle(true);
        histogramDose->GetYaxis()->CenterTitle(true);
        histogramDose->GetXaxis()->SetTitle(xAxisName.c_str());
        histogramDose->GetYaxis()->SetTitle(yAxisName.c_str());
        histogramDose->SetTitle(title.c_str());


        sumHist_SolutionPlusMembranePlusCytoplasm->Draw("HIST");
        sumHist_SolutionPlusMembrane->Draw("HIST SAME");
        histogramDose_FractionFromSolution->Draw("HIST SAME");
        histogramDose->Draw("HIST SAME");

        auto legend = new TLegend(0.67,0.7,0.9,0.9);

        legend->SetHeader("Location of ^{212}Pb Decay","C"); // option "C" allows to center the header
        TLegendEntry *header = (TLegendEntry*)legend->GetListOfPrimitives()->First();
        header->SetTextAlign(22);
        header->SetTextSize(.032);

        legend->AddEntry(histogramDose,"From all locations")->SetTextSize(0.03);
        legend->AddEntry(histogramDose_FractionFromSolution,"Solution")->SetTextSize(0.03);
        legend->AddEntry(sumHist_SolutionPlusMembrane,"Membrane")->SetTextSize(0.03);
        legend->AddEntry(sumHist_SolutionPlusMembranePlusCytoplasm,"Cytoplasm")->SetTextSize(0.03);
        legend->Draw();


        histogramDose_FractionFromSolution->GetYaxis()->SetRangeUser(minY,1.);
        sumHist_SolutionPlusMembrane->GetYaxis()->SetRangeUser(minY,1.);
        sumHist_SolutionPlusMembranePlusCytoplasm->GetYaxis()->SetRangeUser(minY,1.);
        histogramDose->GetYaxis()->SetRangeUser(minY,1.);

        c1->SetLogy();
        c1->Update();

        std::string output = "/Users/idarosenqvist/Desktop/Academics/MasterThesis/Thesis/figures/Results/" + cellGeometry + "/" + cellLine + "/LocationDecays_" + cellLine + "_" + cellComponent + "_" + std::to_string(activity) + "kBq_logy.pdf";
        c1->SaveAs(output.c_str());
    }


}


void PlotHistogramsByLocationDecay()
{
    // std::string cell_component = "Nucleus";
    std::string cell_component = "TotalCell";

    // PlotHistogramsLocationDecay("C4_2", "D12RP", cell_component, 25);
    // PlotHistogramsLocationDecay("C4_2", "D12CP", cell_component, 25);
    // PlotHistogramsLocationDecay("C4_2", "D5RP", cell_component, 25);
    // PlotHistogramsLocationDecay("C4_2", "D5CP", cell_component, 25);

    // PlotHistogramsLocationDecay("PC3_PIP", "D12RP", cell_component, 25);
    // PlotHistogramsLocationDecay("PC3_PIP", "D12CP", cell_component, 25);
    // PlotHistogramsLocationDecay("PC3_PIP", "D5RP", cell_component, 25);
    // PlotHistogramsLocationDecay("PC3_PIP", "D5CP", cell_component, 25);

    PlotHistogramsLocationDecay("PC3_Flu", "D12RP", cell_component, 25);
    // PlotHistogramsLocationDecay("PC3_Flu", "D12CP", cell_component, 25);
    // PlotHistogramsLocationDecay("PC3_Flu", "D5RP", cell_component, 25);
    // PlotHistogramsLocationDecay("PC3_Flu", "D5CP", cell_component, 25);

    // PlotHistogramsLocationDecayLog("C4_2", "D12RP", cell_component, 25);
    // PlotHistogramsLocationDecayLog("C4_2", "D12CP", cell_component, 25);
    // PlotHistogramsLocationDecayLog("C4_2", "D5RP", cell_component, 25);
    // PlotHistogramsLocationDecayLog("C4_2", "D5CP", cell_component, 25);

    // PlotHistogramsLocationDecayLog("PC3_PIP", "D12RP", cell_component, 25);
    // PlotHistogramsLocationDecayLog("PC3_PIP", "D12CP", cell_component, 25);
    // PlotHistogramsLocationDecayLog("PC3_PIP", "D5RP", cell_component, 25);
    // PlotHistogramsLocationDecayLog("PC3_PIP", "D5CP", cell_component, 25);

    // PlotHistogramsLocationDecayLog("PC3_Flu", "D12RP", cell_component, 25);
    // PlotHistogramsLocationDecayLog("PC3_Flu", "D12CP", cell_component, 25);
    // PlotHistogramsLocationDecayLog("PC3_Flu", "D5RP", cell_component, 25);
    // PlotHistogramsLocationDecayLog("PC3_Flu", "D5CP", cell_component, 25);
}
