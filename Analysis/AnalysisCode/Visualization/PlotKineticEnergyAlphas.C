#include <TH1.h>
#include <TCanvas.h>
#include <vector>

void MakePlots(std::string cellLine, std::string cellGeometry, std::string cellComponent, int activity)
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

    int diameter;

    if(cellGeometry=="D12RP"||cellGeometry=="D12CP")
    {
        diameter = 12;
    }
    if(cellGeometry=="D5RP"||cellGeometry=="D5CP")
    {
        diameter = 5;
    }

    std::string cellLine_Name;

    if(cellLine=="C4_2")
    {
        cellLine_Name = "C4-2";
    }
    if(cellLine=="PC3_PIP")
    {
        cellLine_Name = "PC3-PIP";
    }
    if(cellLine=="PC3_Flu")
    {
        cellLine_Name = "PC3-Flu";
    }

    std::string region;

    if(cellComponent=="Nucleus")
    {
        region = "nucleus";
    }
    if(cellComponent=="Membrane")
    {
        region = "membrane";
    }
    if(cellComponent=="Cytoplasm")
    {
        region = "cytoplasm";
    }
    if(cellComponent=="TotalCell")
    {
        region = "cell";
    }

    std::string nucleiDist;
    if(cellGeometry=="D12RP"||cellGeometry=="D5RP")
    {
        nucleiDist = "Uniformly";
    }
    if(cellGeometry=="D12CP"||cellGeometry=="D5CP")
    {
        nucleiDist = "Centrally";
    }

    if(cellLine=="PC3_Flu")
    {

        //-------------------------------
        // Loading activity histograms for scaling
        std::string fileName_activity = "../Output_" + cellGeometry + "/Output_"+ cellLine + "_" + std::to_string(activity) + "kBq.root";
        auto inputFile_activity = std::unique_ptr<TFile>(TFile::Open(fileName_activity.c_str()));

        std::string histogramName_KineticEnergy_FromSolution_activity = "i0_hKineticEnergyAlpha_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_" + cellComponent + "_FromSolution";


        //-----------------------------
        TH1D* histogram_KineticEnergy_FromSolution_activity = 0;
        inputFile_activity->GetObject(histogramName_KineticEnergy_FromSolution_activity.c_str(), histogram_KineticEnergy_FromSolution_activity);
        histogram_KineticEnergy_FromSolution_activity->SetDirectory(0);


        double intSolution_activity = histogram_KineticEnergy_FromSolution_activity->Integral();

        //------------------------------------------------
        // Loading 150kBq histograms
        std::string fileName_150kBq = "../Output_" + cellGeometry + "/Output_"+ cellLine + "_150kBq.root";
        auto inputFile_150kBq = std::unique_ptr<TFile>(TFile::Open(fileName_150kBq.c_str()));


        std::string histogramName_KineticEnergy_FromSolution = "i0_hKineticEnergyAlpha_212Pb_" + cellLine + "_150kBq_" + cellComponent + "_FromSolution";

        //-----------------------------
        TH1D* histogram_KineticEnergy_FromSolution = 0;

        if(cellComponent=="TotalCell")
        {
            inputFile_activity->GetObject(histogramName_KineticEnergy_FromSolution_activity.c_str(), histogram_KineticEnergy_FromSolution);
            histogram_KineticEnergy_FromSolution->SetDirectory(0);
        }
        else
        {
            inputFile_150kBq->GetObject(histogramName_KineticEnergy_FromSolution.c_str(), histogram_KineticEnergy_FromSolution);
            histogram_KineticEnergy_FromSolution->SetDirectory(0);
        }



        //-----------------------------
        double intSolution_150kBq = histogram_KineticEnergy_FromSolution->Integral();


        double scaleSolution = intSolution_activity/intSolution_150kBq;


        histogram_KineticEnergy_FromSolution->Scale(scaleSolution);


        gStyle->SetOptStat(0);

        //------------------------------
        // Create the canvas
        std::string canvasName = "c1_" + cellLine + "_" + std::to_string(activity) + "kBq";
        auto c1 = new TCanvas(canvasName.c_str(), "Kinetic Energy of Impingning Alpha Particle Location of 212-Pb Decay", 600,400);


        histogram_KineticEnergy_FromSolution->SetLineColor(colours[2]);
        // histogram_KineticEnergy_FromSolution->SetLineWidth(2);
        histogram_KineticEnergy_FromSolution->SetFillColorAlpha(colours[2], 0.5);;

        int reBin = 4000;

        histogram_KineticEnergy_FromSolution->Rebin(reBin);


        std::string yAxisName = "Fraction of impinging alpha particles / 40 keV bin";

        double maxY = histogram_KineticEnergy_FromSolution->GetMaximum() + 0.1*histogram_KineticEnergy_FromSolution->GetMaximum();


        std::string xAxisName = "Kinetic energy of alpha particle hitting " + region + " [MeV]";
        std::string title = cellLine_Name + ", " + nucleiDist + " Dist. Nuclei, d_{nuc} = " + Form("%d", diameter) + " #mum, Activity = " + Form("%d", activity)   + "kBq / 1mL";


        histogram_KineticEnergy_FromSolution->GetXaxis()->SetRangeUser(0.,9.);

        histogram_KineticEnergy_FromSolution->GetYaxis()->SetRangeUser(0.,maxY);


        histogram_KineticEnergy_FromSolution->SetTitle(title.c_str());

        histogram_KineticEnergy_FromSolution->GetXaxis()->CenterTitle(true);
        histogram_KineticEnergy_FromSolution->GetYaxis()->CenterTitle(true);
        histogram_KineticEnergy_FromSolution->GetXaxis()->SetTitle(xAxisName.c_str());
        histogram_KineticEnergy_FromSolution->GetYaxis()->SetTitle(yAxisName.c_str());



        histogram_KineticEnergy_FromSolution->Draw("HIST");


        auto legend = new TLegend(0.1,0.75,0.4,0.9);
        legend->SetHeader("Location of ^{212}Pb Decay","C"); // option "C" allows to center the header
        TLegendEntry *header = (TLegendEntry*)legend->GetListOfPrimitives()->First();
        header->SetTextAlign(22);
        header->SetTextSize(.04);

        legend->AddEntry(histogram_KineticEnergy_FromSolution,"Solution")->SetTextSize(0.04);
        legend->Draw();

        c1->Update();

        std::string output = "/Users/idarosenqvist/Desktop/Academics/MasterThesis/Thesis/figures/Results/" + cellGeometry + "/" + cellLine + "/KineticEnergyAlphas_" + cellLine + "_" + cellComponent + "_" + std::to_string(activity) + "kBq.pdf";
        c1->SaveAs(output.c_str());
    }
    else
    {
        //-------------------------------
        // Loading activity histograms for scaling
        std::string fileName_activity = "../Output_" + cellGeometry + "/Output_"+ cellLine + "_" + std::to_string(activity) + "kBq.root";
        auto inputFile_activity = std::unique_ptr<TFile>(TFile::Open(fileName_activity.c_str()));

        std::string histogramName_KineticEnergy_FromSolution_activity = "i0_hKineticEnergyAlpha_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_" + cellComponent + "_FromSolution";
        std::string histogramName_KineticEnergy_FromMembrane_activity = "i0_hKineticEnergyAlpha_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_" + cellComponent + "_FromMembrane";
        std::string histogramName_KineticEnergy_FromCytoplasm_activity = "i0_hKineticEnergyAlpha_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_" + cellComponent + "_FromCytoplasm";


        //-----------------------------
        TH1D* histogram_KineticEnergy_FromSolution_activity = 0;
        inputFile_activity->GetObject(histogramName_KineticEnergy_FromSolution_activity.c_str(), histogram_KineticEnergy_FromSolution_activity);
        histogram_KineticEnergy_FromSolution_activity->SetDirectory(0);

        TH1D* histogram_KineticEnergy_FromMembrane_activity = 0;
        inputFile_activity->GetObject(histogramName_KineticEnergy_FromMembrane_activity.c_str(), histogram_KineticEnergy_FromMembrane_activity);
        histogram_KineticEnergy_FromMembrane_activity->SetDirectory(0);

        TH1D* histogram_KineticEnergy_FromCytoplasm_activity = 0;
        inputFile_activity->GetObject(histogramName_KineticEnergy_FromCytoplasm_activity.c_str(), histogram_KineticEnergy_FromCytoplasm_activity);
        histogram_KineticEnergy_FromCytoplasm_activity->SetDirectory(0);


        double intSolution_activity = histogram_KineticEnergy_FromSolution_activity->Integral();
        double intMembrane_activity = histogram_KineticEnergy_FromMembrane_activity->Integral();
        double intCytoplasm_activity = histogram_KineticEnergy_FromCytoplasm_activity->Integral();

        //------------------------------------------------
        // Loading 150kBq histograms
        std::string fileName_150kBq = "../Output_" + cellGeometry + "/Output_"+ cellLine + "_150kBq.root";
        auto inputFile_150kBq = std::unique_ptr<TFile>(TFile::Open(fileName_150kBq.c_str()));


        std::string histogramName_KineticEnergy_FromSolution = "i0_hKineticEnergyAlpha_212Pb_" + cellLine + "_150kBq_" + cellComponent + "_FromSolution";
        std::string histogramName_KineticEnergy_FromMembrane = "i0_hKineticEnergyAlpha_212Pb_" + cellLine + "_150kBq_" + cellComponent + "_FromMembrane";
        std::string histogramName_KineticEnergy_FromCytoplasm = "i0_hKineticEnergyAlpha_212Pb_" + cellLine + "_150kBq_" + cellComponent + "_FromCytoplasm";

        TH1D* histogram_KineticEnergy_FromSolution = 0;
        TH1D* histogram_KineticEnergy_FromMembrane = 0;
        TH1D* histogram_KineticEnergy_FromCytoplasm = 0;

        if(cellComponent=="TotalCell")
        {
            inputFile_activity->GetObject(histogramName_KineticEnergy_FromSolution_activity.c_str(), histogram_KineticEnergy_FromSolution);
            histogram_KineticEnergy_FromSolution->SetDirectory(0);

            inputFile_activity->GetObject(histogramName_KineticEnergy_FromMembrane_activity.c_str(), histogram_KineticEnergy_FromMembrane);
            histogram_KineticEnergy_FromMembrane->SetDirectory(0);

            inputFile_activity->GetObject(histogramName_KineticEnergy_FromCytoplasm_activity.c_str(), histogram_KineticEnergy_FromCytoplasm);
            histogram_KineticEnergy_FromCytoplasm->SetDirectory(0);
        }
        else
        {
            inputFile_150kBq->GetObject(histogramName_KineticEnergy_FromSolution.c_str(), histogram_KineticEnergy_FromSolution);
            histogram_KineticEnergy_FromSolution->SetDirectory(0);

            inputFile_150kBq->GetObject(histogramName_KineticEnergy_FromMembrane.c_str(), histogram_KineticEnergy_FromMembrane);
            histogram_KineticEnergy_FromMembrane->SetDirectory(0);

            inputFile_150kBq->GetObject(histogramName_KineticEnergy_FromCytoplasm.c_str(), histogram_KineticEnergy_FromCytoplasm);
            histogram_KineticEnergy_FromCytoplasm->SetDirectory(0);
        }

        //-----------------------------
        double intSolution_150kBq = histogram_KineticEnergy_FromSolution->Integral();
        double intMembrane_150kBq = histogram_KineticEnergy_FromMembrane->Integral();
        double intCytoplasm_150kBq = histogram_KineticEnergy_FromCytoplasm->Integral();


        double scaleSolution = intSolution_activity/intSolution_150kBq;
        double scaleMembrane = intMembrane_activity/intMembrane_150kBq;
        double scaleCytoplasm = intCytoplasm_activity/intCytoplasm_150kBq;


        histogram_KineticEnergy_FromSolution->Scale(scaleSolution);
        histogram_KineticEnergy_FromMembrane->Scale(scaleMembrane);
        histogram_KineticEnergy_FromCytoplasm->Scale(scaleCytoplasm);

        gStyle->SetOptStat(0);

        //------------------------------
        // Create the canvas
        std::string canvasName = "c1_" + cellLine + "_" + std::to_string(activity) + "kBq";
        auto c1 = new TCanvas(canvasName.c_str(), "Kinetic Energy of Impingning Alpha Particle Location of 212-Pb Decay", 600,400);


        histogram_KineticEnergy_FromSolution->SetLineColor(colours[2]);
        // histogram_KineticEnergy_FromSolution->SetLineWidth(2);
        histogram_KineticEnergy_FromSolution->SetFillColorAlpha(colours[2], 0.5);

        histogram_KineticEnergy_FromMembrane->SetLineColor(colours[3]);
        // histogram_KineticEnergy_FromMembrane->SetLineWidth(2);
        histogram_KineticEnergy_FromMembrane->SetFillColorAlpha(colours[3], 0.5);

        histogram_KineticEnergy_FromCytoplasm->SetLineColor(colours[4]);
        // histogram_KineticEnergy_FromCytoplasm->SetLineWidth(2);
        histogram_KineticEnergy_FromCytoplasm->SetFillColorAlpha(colours[4], 0.5);

        int reBin = 2000;

        histogram_KineticEnergy_FromSolution->Rebin(reBin);
        histogram_KineticEnergy_FromMembrane->Rebin(reBin);
        histogram_KineticEnergy_FromCytoplasm->Rebin(reBin);

        std::string yAxisName = "Fraction of impinging alpha particles / 20 keV bin";

        double maxY_Mem = histogram_KineticEnergy_FromMembrane->GetMaximum() + 0.1*histogram_KineticEnergy_FromMembrane->GetMaximum();
        double maxY_Sol = histogram_KineticEnergy_FromSolution->GetMaximum() + 0.1*histogram_KineticEnergy_FromSolution->GetMaximum();
        double maxY_Cyt = histogram_KineticEnergy_FromCytoplasm->GetMaximum() + 0.1*histogram_KineticEnergy_FromCytoplasm->GetMaximum();

        std::vector<double> maximums = {maxY_Mem,maxY_Cyt,maxY_Sol};

        double maxY = 0.;
        for(auto& entry: maximums)
        {
            if(entry>maxY)
            {
                maxY = entry;
            }
        }

        std::string xAxisName = "Kinetic energy of alpha particle hitting " + region + " [MeV]";
        std::string title = cellLine_Name + ", " + nucleiDist + " Dist. Nuclei, d_{nuc} = " + Form("%d", diameter) + " #mum, Activity = " + Form("%d", activity)   + "kBq / 1mL";


        histogram_KineticEnergy_FromSolution->GetXaxis()->SetRangeUser(0.,9.);
        histogram_KineticEnergy_FromMembrane->GetXaxis()->SetRangeUser(0.,9.);
        histogram_KineticEnergy_FromCytoplasm->GetXaxis()->SetRangeUser(0.,9.);

        histogram_KineticEnergy_FromSolution->GetYaxis()->SetRangeUser(0.,maxY);
        histogram_KineticEnergy_FromMembrane->GetYaxis()->SetRangeUser(0.,maxY);
        histogram_KineticEnergy_FromCytoplasm->GetYaxis()->SetRangeUser(0.,maxY);


        histogram_KineticEnergy_FromSolution->SetTitle(title.c_str());
        histogram_KineticEnergy_FromMembrane->SetTitle(title.c_str());
        histogram_KineticEnergy_FromCytoplasm->SetTitle(title.c_str());

        histogram_KineticEnergy_FromSolution->GetXaxis()->CenterTitle(true);
        histogram_KineticEnergy_FromSolution->GetYaxis()->CenterTitle(true);
        histogram_KineticEnergy_FromSolution->GetXaxis()->SetTitle(xAxisName.c_str());
        histogram_KineticEnergy_FromSolution->GetYaxis()->SetTitle(yAxisName.c_str());

        histogram_KineticEnergy_FromMembrane->GetXaxis()->CenterTitle(true);
        histogram_KineticEnergy_FromMembrane->GetYaxis()->CenterTitle(true);
        histogram_KineticEnergy_FromMembrane->GetXaxis()->SetTitle(xAxisName.c_str());
        histogram_KineticEnergy_FromMembrane->GetYaxis()->SetTitle(yAxisName.c_str());

        histogram_KineticEnergy_FromCytoplasm->GetXaxis()->CenterTitle(true);
        histogram_KineticEnergy_FromCytoplasm->GetYaxis()->CenterTitle(true);
        histogram_KineticEnergy_FromCytoplasm->GetXaxis()->SetTitle(xAxisName.c_str());
        histogram_KineticEnergy_FromCytoplasm->GetYaxis()->SetTitle(yAxisName.c_str());

        histogram_KineticEnergy_FromSolution->Draw("HIST");
        histogram_KineticEnergy_FromMembrane->Draw("HIST SAME");
        histogram_KineticEnergy_FromCytoplasm->Draw("HIST SAME");


        auto legend = new TLegend(0.1,0.7,0.4,0.9);
        legend->SetHeader("Location of ^{212}Pb Decay","C"); // option "C" allows to center the header
        TLegendEntry *header = (TLegendEntry*)legend->GetListOfPrimitives()->First();
        header->SetTextAlign(22);
        header->SetTextSize(.04);

        legend->AddEntry(histogram_KineticEnergy_FromSolution,"Solution")->SetTextSize(0.04);
        legend->AddEntry(histogram_KineticEnergy_FromMembrane,"Membrane")->SetTextSize(0.04);
        legend->AddEntry(histogram_KineticEnergy_FromCytoplasm,"Cytoplasm")->SetTextSize(0.04);
        legend->Draw();

        std::string output = "/Users/idarosenqvist/Desktop/Academics/MasterThesis/Thesis/figures/Results/" + cellGeometry + "/" + cellLine + "/KineticEnergyAlphas_" + cellLine + "_" + cellComponent + "_" + std::to_string(activity) + "kBq.pdf";

        c1->Update();
        c1->SaveAs(output.c_str());

    }
}


void MakePlotsLog(std::string cellLine, std::string cellGeometry, std::string cellComponent, int activity)
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

    int diameter;

    if(cellGeometry=="D12RP"||cellGeometry=="D12CP")
    {
        diameter = 12;
    }
    if(cellGeometry=="D5RP"||cellGeometry=="D5CP")
    {
        diameter = 5;
    }

    std::string cellLine_Name;

    if(cellLine=="C4_2")
    {
        cellLine_Name = "C4-2";
    }
    if(cellLine=="PC3_PIP")
    {
        cellLine_Name = "PC3-PIP";
    }
    if(cellLine=="PC3_Flu")
    {
        cellLine_Name = "PC3-Flu";
    }

    std::string region;

    if(cellComponent=="Nucleus")
    {
        region = "nucleus";
    }
    if(cellComponent=="Membrane")
    {
        region = "membrane";
    }
    if(cellComponent=="Cytoplasm")
    {
        region = "cytoplasm";
    }
    if(cellComponent=="TotalCell")
    {
        region = "cell";
    }

    std::string nucleiDist;
    if(cellGeometry=="D12RP"||cellGeometry=="D5RP")
    {
        nucleiDist = "Uniformly";
    }
    if(cellGeometry=="D12CP"||cellGeometry=="D5CP")
    {
        nucleiDist = "Centrally";
    }

    if(cellLine=="PC3_Flu")
    {
        //-------------------------------
        // Loading activity histograms for scaling
        std::string fileName_activity = "../Output_" + cellGeometry + "/Output_"+ cellLine + "_" + std::to_string(activity) + "kBq.root";
        auto inputFile_activity = std::unique_ptr<TFile>(TFile::Open(fileName_activity.c_str()));

        std::string histogramName_KineticEnergy_FromSolution_activity = "i0_hKineticEnergyAlpha_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_" + cellComponent + "_FromSolution";


        //-----------------------------
        TH1D* histogram_KineticEnergy_FromSolution_activity = 0;
        inputFile_activity->GetObject(histogramName_KineticEnergy_FromSolution_activity.c_str(), histogram_KineticEnergy_FromSolution_activity);
        histogram_KineticEnergy_FromSolution_activity->SetDirectory(0);


        double intSolution_activity = histogram_KineticEnergy_FromSolution_activity->Integral();

        //------------------------------------------------
        // Loading 150kBq histograms
        std::string fileName_150kBq = "../Output_" + cellGeometry + "/Output_"+ cellLine + "_150kBq.root";
        auto inputFile_150kBq = std::unique_ptr<TFile>(TFile::Open(fileName_150kBq.c_str()));


        std::string histogramName_KineticEnergy_FromSolution = "i0_hKineticEnergyAlpha_212Pb_" + cellLine + "_150kBq_" + cellComponent + "_FromSolution";

        //-----------------------------
        TH1D* histogram_KineticEnergy_FromSolution = 0;

        if(cellComponent=="TotalCell")
        {
            inputFile_activity->GetObject(histogramName_KineticEnergy_FromSolution_activity.c_str(), histogram_KineticEnergy_FromSolution);
            histogram_KineticEnergy_FromSolution->SetDirectory(0);
        }
        else
        {
            inputFile_150kBq->GetObject(histogramName_KineticEnergy_FromSolution.c_str(), histogram_KineticEnergy_FromSolution);
            histogram_KineticEnergy_FromSolution->SetDirectory(0);
        }



        //-----------------------------
        double intSolution_150kBq = histogram_KineticEnergy_FromSolution->Integral();


        double scaleSolution = intSolution_activity/intSolution_150kBq;


        histogram_KineticEnergy_FromSolution->Scale(scaleSolution);


        gStyle->SetOptStat(0);

        //------------------------------
        // Create the canvas
        std::string canvasName = "c1_" + cellLine + "_" + std::to_string(activity) + "kBq";
        auto c1 = new TCanvas(canvasName.c_str(), "Kinetic Energy of Impingning Alpha Particle Location of 212-Pb Decay", 600,400);


        histogram_KineticEnergy_FromSolution->SetLineColor(colours[2]);
        // histogram_KineticEnergy_FromSolution->SetLineWidth(2);
        histogram_KineticEnergy_FromSolution->SetFillColorAlpha(colours[2], 0.5);

        int reBin = 4000;

        histogram_KineticEnergy_FromSolution->Rebin(reBin);


        std::string yAxisName = "Fraction of impinging alpha particles / 40 keV bin";

        double maxY = histogram_KineticEnergy_FromSolution->GetMaximum() + 0.1*histogram_KineticEnergy_FromSolution->GetMaximum();


        std::string xAxisName = "Kinetic energy of alpha particle hitting " + region + " [MeV]";
        std::string title = cellLine_Name + ", " + nucleiDist + " Dist. Nuclei, d_{nuc} = " + Form("%d", diameter) + " #mum, Activity = " + Form("%d", activity)   + "kBq / 1mL";


        histogram_KineticEnergy_FromSolution->GetXaxis()->SetRangeUser(0.,9.);

        histogram_KineticEnergy_FromSolution->GetYaxis()->SetRangeUser(1.e-7,1.);


        histogram_KineticEnergy_FromSolution->SetTitle(title.c_str());

        histogram_KineticEnergy_FromSolution->GetXaxis()->CenterTitle(true);
        histogram_KineticEnergy_FromSolution->GetYaxis()->CenterTitle(true);
        histogram_KineticEnergy_FromSolution->GetXaxis()->SetTitle(xAxisName.c_str());
        histogram_KineticEnergy_FromSolution->GetYaxis()->SetTitle(yAxisName.c_str());



        histogram_KineticEnergy_FromSolution->Draw("HIST");


        auto legend = new TLegend(0.1,0.75,0.4,0.9);
        legend->SetHeader("Location of ^{212}Pb Decay","C"); // option "C" allows to center the header
        TLegendEntry *header = (TLegendEntry*)legend->GetListOfPrimitives()->First();
        header->SetTextAlign(22);
        header->SetTextSize(.04);

        legend->AddEntry(histogram_KineticEnergy_FromSolution,"Solution")->SetTextSize(0.04);
        legend->Draw();

        c1->SetLogy();
        c1->Update();

        std::string output = "/Users/idarosenqvist/Desktop/Academics/MasterThesis/Thesis/figures/Results/" + cellGeometry + "/" + cellLine + "/KineticEnergyAlphas_" + cellLine + "_" + cellComponent + "_" + std::to_string(activity) + "kBq_logy.pdf";
        c1->SaveAs(output.c_str());
    }
    else
    {
        //-------------------------------
        // Loading activity histograms for scaling
        std::string fileName_activity = "../Output_" + cellGeometry + "/Output_"+ cellLine + "_" + std::to_string(activity) + "kBq.root";
        auto inputFile_activity = std::unique_ptr<TFile>(TFile::Open(fileName_activity.c_str()));

        std::string histogramName_KineticEnergy_FromSolution_activity = "i0_hKineticEnergyAlpha_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_" + cellComponent + "_FromSolution";
        std::string histogramName_KineticEnergy_FromMembrane_activity = "i0_hKineticEnergyAlpha_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_" + cellComponent + "_FromMembrane";
        std::string histogramName_KineticEnergy_FromCytoplasm_activity = "i0_hKineticEnergyAlpha_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_" + cellComponent + "_FromCytoplasm";


        //-----------------------------
        TH1D* histogram_KineticEnergy_FromSolution_activity = 0;
        inputFile_activity->GetObject(histogramName_KineticEnergy_FromSolution_activity.c_str(), histogram_KineticEnergy_FromSolution_activity);
        histogram_KineticEnergy_FromSolution_activity->SetDirectory(0);

        TH1D* histogram_KineticEnergy_FromMembrane_activity = 0;
        inputFile_activity->GetObject(histogramName_KineticEnergy_FromMembrane_activity.c_str(), histogram_KineticEnergy_FromMembrane_activity);
        histogram_KineticEnergy_FromMembrane_activity->SetDirectory(0);

        TH1D* histogram_KineticEnergy_FromCytoplasm_activity = 0;
        inputFile_activity->GetObject(histogramName_KineticEnergy_FromCytoplasm_activity.c_str(), histogram_KineticEnergy_FromCytoplasm_activity);
        histogram_KineticEnergy_FromCytoplasm_activity->SetDirectory(0);


        double intSolution_activity = histogram_KineticEnergy_FromSolution_activity->Integral();
        double intMembrane_activity = histogram_KineticEnergy_FromMembrane_activity->Integral();
        double intCytoplasm_activity = histogram_KineticEnergy_FromCytoplasm_activity->Integral();

        //------------------------------------------------
        // Loading 150kBq histograms
        std::string fileName_150kBq = "../Output_" + cellGeometry + "/Output_"+ cellLine + "_150kBq.root";
        auto inputFile_150kBq = std::unique_ptr<TFile>(TFile::Open(fileName_150kBq.c_str()));


        std::string histogramName_KineticEnergy_FromSolution = "i0_hKineticEnergyAlpha_212Pb_" + cellLine + "_150kBq_" + cellComponent + "_FromSolution";
        std::string histogramName_KineticEnergy_FromMembrane = "i0_hKineticEnergyAlpha_212Pb_" + cellLine + "_150kBq_" + cellComponent + "_FromMembrane";
        std::string histogramName_KineticEnergy_FromCytoplasm = "i0_hKineticEnergyAlpha_212Pb_" + cellLine + "_150kBq_" + cellComponent + "_FromCytoplasm";

        TH1D* histogram_KineticEnergy_FromSolution = 0;
        TH1D* histogram_KineticEnergy_FromMembrane = 0;
        TH1D* histogram_KineticEnergy_FromCytoplasm = 0;

        if(cellComponent=="TotalCell")
        {
            inputFile_activity->GetObject(histogramName_KineticEnergy_FromSolution_activity.c_str(), histogram_KineticEnergy_FromSolution);
            histogram_KineticEnergy_FromSolution->SetDirectory(0);

            inputFile_activity->GetObject(histogramName_KineticEnergy_FromMembrane_activity.c_str(), histogram_KineticEnergy_FromMembrane);
            histogram_KineticEnergy_FromMembrane->SetDirectory(0);

            inputFile_activity->GetObject(histogramName_KineticEnergy_FromCytoplasm_activity.c_str(), histogram_KineticEnergy_FromCytoplasm);
            histogram_KineticEnergy_FromCytoplasm->SetDirectory(0);
        }
        else
        {
            inputFile_150kBq->GetObject(histogramName_KineticEnergy_FromSolution.c_str(), histogram_KineticEnergy_FromSolution);
            histogram_KineticEnergy_FromSolution->SetDirectory(0);

            inputFile_150kBq->GetObject(histogramName_KineticEnergy_FromMembrane.c_str(), histogram_KineticEnergy_FromMembrane);
            histogram_KineticEnergy_FromMembrane->SetDirectory(0);

            inputFile_150kBq->GetObject(histogramName_KineticEnergy_FromCytoplasm.c_str(), histogram_KineticEnergy_FromCytoplasm);
            histogram_KineticEnergy_FromCytoplasm->SetDirectory(0);
        }


        //-----------------------------
        double intSolution_150kBq = histogram_KineticEnergy_FromSolution->Integral();
        double intMembrane_150kBq = histogram_KineticEnergy_FromMembrane->Integral();
        double intCytoplasm_150kBq = histogram_KineticEnergy_FromCytoplasm->Integral();


        double scaleSolution = intSolution_activity/intSolution_150kBq;
        double scaleMembrane = intMembrane_activity/intMembrane_150kBq;
        double scaleCytoplasm = intCytoplasm_activity/intCytoplasm_150kBq;


        histogram_KineticEnergy_FromSolution->Scale(scaleSolution);
        histogram_KineticEnergy_FromMembrane->Scale(scaleMembrane);
        histogram_KineticEnergy_FromCytoplasm->Scale(scaleCytoplasm);

        gStyle->SetOptStat(0);

        //------------------------------
        // Create the canvas
        std::string canvasName = "c1_" + cellLine + "_" + std::to_string(activity) + "kBq";
        auto c1 = new TCanvas(canvasName.c_str(), "Kinetic Energy of Impingning Alpha Particle Location of 212-Pb Decay", 600, 400);


        histogram_KineticEnergy_FromSolution->SetLineColor(colours[2]);
        // histogram_KineticEnergy_FromSolution->SetLineWidth(2);
        histogram_KineticEnergy_FromSolution->SetFillColorAlpha(colours[2], 0.5);

        histogram_KineticEnergy_FromMembrane->SetLineColor(colours[3]);
        // histogram_KineticEnergy_FromMembrane->SetLineWidth(2);
        histogram_KineticEnergy_FromMembrane->SetFillColorAlpha(colours[3], 0.5);

        histogram_KineticEnergy_FromCytoplasm->SetLineColor(colours[4]);
        // histogram_KineticEnergy_FromCytoplasm->SetLineWidth(2);
        histogram_KineticEnergy_FromCytoplasm->SetFillColorAlpha(colours[4], 0.5);

        int reBin = 2000;

        histogram_KineticEnergy_FromSolution->Rebin(reBin);
        histogram_KineticEnergy_FromMembrane->Rebin(reBin);
        histogram_KineticEnergy_FromCytoplasm->Rebin(reBin);

        std::string yAxisName = "Fraction of impinging alpha particles / 20 keV bin";

        double maxY_Mem = histogram_KineticEnergy_FromMembrane->GetMaximum() + 0.1*histogram_KineticEnergy_FromMembrane->GetMaximum();
        double maxY_Sol = histogram_KineticEnergy_FromSolution->GetMaximum() + 0.1*histogram_KineticEnergy_FromSolution->GetMaximum();
        double maxY_Cyt = histogram_KineticEnergy_FromCytoplasm->GetMaximum() + 0.1*histogram_KineticEnergy_FromCytoplasm->GetMaximum();

        std::vector<double> maximums = {maxY_Mem,maxY_Cyt,maxY_Sol};

        double maxY = 0.;
        for(auto& entry: maximums)
        {
            if(entry>maxY)
            {
                maxY = entry;
            }
        }

        std::string xAxisName = "Kinetic energy of alpha particle hitting " + region + " [MeV]";
        std::string title = cellLine_Name + ", " + nucleiDist + " Dist. Nuclei, d_{nuc} = " + Form("%d", diameter) + " #mum, Activity = " + Form("%d", activity)   + "kBq / 1mL";


        histogram_KineticEnergy_FromSolution->GetXaxis()->SetRangeUser(0.,9.);
        histogram_KineticEnergy_FromMembrane->GetXaxis()->SetRangeUser(0.,9.);
        histogram_KineticEnergy_FromCytoplasm->GetXaxis()->SetRangeUser(0.,9.);

        histogram_KineticEnergy_FromSolution->GetYaxis()->SetRangeUser(1.e-7,1.);
        histogram_KineticEnergy_FromMembrane->GetYaxis()->SetRangeUser(1.e-7,1.);
        histogram_KineticEnergy_FromCytoplasm->GetYaxis()->SetRangeUser(1.e-7,1.);


        histogram_KineticEnergy_FromSolution->SetTitle(title.c_str());
        histogram_KineticEnergy_FromMembrane->SetTitle(title.c_str());
        histogram_KineticEnergy_FromCytoplasm->SetTitle(title.c_str());

        histogram_KineticEnergy_FromSolution->GetXaxis()->CenterTitle(true);
        histogram_KineticEnergy_FromSolution->GetYaxis()->CenterTitle(true);
        histogram_KineticEnergy_FromSolution->GetXaxis()->SetTitle(xAxisName.c_str());
        histogram_KineticEnergy_FromSolution->GetYaxis()->SetTitle(yAxisName.c_str());

        histogram_KineticEnergy_FromMembrane->GetXaxis()->CenterTitle(true);
        histogram_KineticEnergy_FromMembrane->GetYaxis()->CenterTitle(true);
        histogram_KineticEnergy_FromMembrane->GetXaxis()->SetTitle(xAxisName.c_str());
        histogram_KineticEnergy_FromMembrane->GetYaxis()->SetTitle(yAxisName.c_str());

        histogram_KineticEnergy_FromCytoplasm->GetXaxis()->CenterTitle(true);
        histogram_KineticEnergy_FromCytoplasm->GetYaxis()->CenterTitle(true);
        histogram_KineticEnergy_FromCytoplasm->GetXaxis()->SetTitle(xAxisName.c_str());
        histogram_KineticEnergy_FromCytoplasm->GetYaxis()->SetTitle(yAxisName.c_str());

        histogram_KineticEnergy_FromSolution->Draw("HIST");
        histogram_KineticEnergy_FromMembrane->Draw("HIST SAME");
        histogram_KineticEnergy_FromCytoplasm->Draw("HIST SAME");


        auto legend = new TLegend(0.1,0.7,0.4,0.9);
        legend->SetHeader("Location of ^{212}Pb Decay","C"); // option "C" allows to center the header
        TLegendEntry *header = (TLegendEntry*)legend->GetListOfPrimitives()->First();
        header->SetTextAlign(22);
        header->SetTextSize(.04);

        legend->AddEntry(histogram_KineticEnergy_FromSolution,"Solution")->SetTextSize(0.04);
        legend->AddEntry(histogram_KineticEnergy_FromMembrane,"Membrane")->SetTextSize(0.04);
        legend->AddEntry(histogram_KineticEnergy_FromCytoplasm,"Cytoplasm")->SetTextSize(0.04);
        legend->Draw();

        std::string output = "/Users/idarosenqvist/Desktop/Academics/MasterThesis/Thesis/figures/Results/" + cellGeometry + "/" + cellLine + "/KineticEnergyAlphas_" + cellLine + "_" + cellComponent + "_" + std::to_string(activity) + "kBq_logy.pdf";

        c1->SetLogy();
        c1->Update();
        c1->SaveAs(output.c_str());

    }
}

void PlotKineticEnergyAlphas()
{
    // std::string cellLine = "C4_2";
    // std::string cellLine = "PC3_PIP";
    std::string cellLine = "PC3_Flu";

    MakePlots(cellLine, "D12RP", "Nucleus", 25);
    MakePlots(cellLine, "D12CP", "Nucleus", 25);
    MakePlots(cellLine, "D5CP", "Nucleus", 25);
    MakePlots(cellLine, "D5RP", "Nucleus", 25);

    MakePlotsLog(cellLine, "D12RP", "Nucleus", 25);
    MakePlotsLog(cellLine, "D12CP", "Nucleus", 25);
    MakePlotsLog(cellLine, "D5CP", "Nucleus", 25);
    MakePlotsLog(cellLine, "D5RP", "Nucleus", 25);

    MakePlots(cellLine, "D12RP", "TotalCell", 25);
    MakePlots(cellLine, "D12CP", "TotalCell", 25);
    MakePlots(cellLine, "D5CP", "TotalCell", 25);
    MakePlots(cellLine, "D5RP", "TotalCell", 25);

    MakePlotsLog(cellLine, "D12RP", "TotalCell", 25);
    MakePlotsLog(cellLine, "D12CP", "TotalCell", 25);
    MakePlotsLog(cellLine, "D5CP", "TotalCell", 25);
    MakePlotsLog(cellLine, "D5RP", "TotalCell", 25);
}