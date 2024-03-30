
void MakePlots(std::string cellLine, std::string cellGeometry, std::string cellComponent, int activity)
{
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


    std::vector<int> colours;
    colours.push_back(kCyan+1);
    colours.push_back(kGreen+2);
    colours.push_back(kOrange-3);
    colours.push_back(kViolet);
    colours.push_back(kRed);
    colours.push_back(kBlue);
    colours.push_back(kOrange-6);
    colours.push_back(kRed);


    //------------------------------------------------
    std::string fileName = "../Output_" + cellGeometry + "/Output_"+ cellLine + "_" + std::to_string(activity) + "kBq.root";
    auto inputFile = std::unique_ptr<TFile>(TFile::Open(fileName.c_str()));


    std::string histogramName_KineticEnergy_FromSolution = "i0_hKineticEnergyAlpha_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_" + cellComponent + "_FromSolution";
    std::string histogramName_KineticEnergy_FromMembrane = "i0_hKineticEnergyAlpha_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_" + cellComponent + "_FromMembrane";
    std::string histogramName_KineticEnergy_FromCytoplasm = "i0_hKineticEnergyAlpha_212Pb_" + cellLine + "_" + std::to_string(activity) + "kBq_" + cellComponent + "_FromCytoplasm";

    //-----------------------------
    TH1D* histogram_KineticEnergy_FromSolution = 0;
    inputFile->GetObject(histogramName_KineticEnergy_FromSolution.c_str(), histogram_KineticEnergy_FromSolution);
    histogram_KineticEnergy_FromSolution->SetDirectory(0);

    TH1D* histogram_KineticEnergy_FromMembrane = 0;
    inputFile->GetObject(histogramName_KineticEnergy_FromMembrane.c_str(), histogram_KineticEnergy_FromMembrane);
    histogram_KineticEnergy_FromMembrane->SetDirectory(0);

    TH1D* histogram_KineticEnergy_FromCytoplasm = 0;
    inputFile->GetObject(histogramName_KineticEnergy_FromCytoplasm.c_str(), histogram_KineticEnergy_FromCytoplasm);
    histogram_KineticEnergy_FromCytoplasm->SetDirectory(0);

    gStyle->SetOptStat(0);

    //------------------------------
    // Create the canvas
    std::string canvasName = "c1_" + cellLine + "_" + std::to_string(activity) + "kBq";
    auto c1 = new TCanvas(canvasName.c_str(), "Kinetic Energy of Impingning Alpha Particle Location of 212-Pb Decay", 600, 400);


    histogram_KineticEnergy_FromSolution->SetLineColor(colours[0]);
    histogram_KineticEnergy_FromSolution->SetLineWidth(2);
    histogram_KineticEnergy_FromSolution->SetFillColorAlpha(colours[0], 0.3);

    histogram_KineticEnergy_FromMembrane->SetLineColor(colours[1]);
    histogram_KineticEnergy_FromMembrane->SetLineWidth(2);
    histogram_KineticEnergy_FromMembrane->SetFillColorAlpha(colours[1], 0.3);

    histogram_KineticEnergy_FromCytoplasm->SetLineColor(colours[3]);
    histogram_KineticEnergy_FromCytoplasm->SetLineWidth(2);
    histogram_KineticEnergy_FromCytoplasm->SetFillColorAlpha(colours[3], 0.3);

    int reBin;

    double xMin;
    double xMax;
    double yMin;
    double yMax;

    std::string yAxisName;

    if(cellLine=="C4_2")
    {
        if(cellGeometry=="D12RP"||cellGeometry=="D12CP")
        {
            xMin = 1.e-10;
            xMax = 9.0;
            yMin = 1.e-4;
            yMax = 0.09;

            reBin = 2000;

            yAxisName = "Fraction of impinging alpha particles / 20 keV bin";
        }
        if(cellGeometry=="D5RP"||cellGeometry=="D5CP")
        {
            xMin = 1.e-10;
            xMax = 9.0;
            yMin = 1.e-4;
            yMax = 0.12;

            reBin = 4000;

            yAxisName = "Fraction of impinging alpha particles / 40 keV bin";
        }
    }
    if(cellLine=="PC3_PIP")
    {
        if(cellGeometry=="D12RP"||cellGeometry=="D12CP")
        {
            xMin = 1.e-10;
            xMax = 9.0;
            yMin = 1.e-4;
            yMax = 0.09;

            reBin = 2000;

            yAxisName = "Fraction of impinging alpha particles / 20 keV bin";
        }
        if(cellGeometry=="D5RP"||cellGeometry=="D5CP")
        {
            xMin = 1.e-10;
            xMax = 9.0;
            yMin = 1.e-4;
            yMax = 0.11;

            reBin = 4000;

            yAxisName = "Fraction of impinging alpha particles / 40 keV bin";
        }
    }

    std::string xAxisName = "Kinetic energy of alpha particle hitting " + region + " [MeV]";
    std::string title = cellLine_Name + ", " + nucleiDist + " Dist. Nuclei, d_{nuc} = " + Form("%d", diameter) + " #mum, Activity = " + Form("%d", activity)   + "kBq / 1mL";

    histogram_KineticEnergy_FromSolution->Rebin(reBin);
    histogram_KineticEnergy_FromMembrane->Rebin(reBin);
    histogram_KineticEnergy_FromCytoplasm->Rebin(reBin);

    histogram_KineticEnergy_FromSolution->GetXaxis()->SetRangeUser(xMin,xMax);
    histogram_KineticEnergy_FromMembrane->GetXaxis()->SetRangeUser(xMin,xMax);
    histogram_KineticEnergy_FromCytoplasm->GetXaxis()->SetRangeUser(xMin,xMax);

    histogram_KineticEnergy_FromSolution->GetYaxis()->SetRangeUser(yMin,yMax);
    histogram_KineticEnergy_FromMembrane->GetYaxis()->SetRangeUser(yMin,yMax);
    histogram_KineticEnergy_FromCytoplasm->GetYaxis()->SetRangeUser(yMin,yMax);


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


    auto legend = new TLegend(0.67,0.7,0.9,0.9);
    legend->SetHeader("Location of ^{212}Pb Decay","C"); // option "C" allows to center the header
    TLegendEntry *header = (TLegendEntry*)legend->GetListOfPrimitives()->First();
    header->SetTextAlign(22);
    header->SetTextSize(.032);

    legend->AddEntry(histogram_KineticEnergy_FromSolution,"Solution")->SetTextSize(0.03);
    legend->AddEntry(histogram_KineticEnergy_FromMembrane,"Membrane")->SetTextSize(0.03);
    legend->AddEntry(histogram_KineticEnergy_FromCytoplasm,"Cytoplasm")->SetTextSize(0.03);
    legend->Draw();

    c1->Update();

    std::string output = "Figures_" + cellGeometry + "/KineticEnergyAlphas_" + cellLine + "_" + std::to_string(activity) + "kBq.pdf";
    c1->SaveAs(output.c_str());

    //----------------------------
    // Log scale

    if(cellLine=="C4_2")
    {
        if(cellGeometry=="D12RP"||cellGeometry=="D12CP")
        {
            xMin = 1.e-10;
            xMax = 9.0;
            yMin = 1.e-7;
            yMax = 10.0;
        }
        if(cellGeometry=="D5RP"||cellGeometry=="D5CP")
        {
            xMin = 1.e-10;
            xMax = 9.0;
            yMin = 1.e-7;
            yMax = 10.0;
        }
    }
    if(cellLine=="PC3_PIP")
    {
        if(cellGeometry=="D12RP"||cellGeometry=="D12CP")
        {
            xMin = 1.e-10;
            xMax = 9.0;
            yMin = 1.e-7;
            yMax = 10.0;
        }
        if(cellGeometry=="D5RP"||cellGeometry=="D5CP")
        {
            xMin = 1.e-10;
            xMax = 9.0;
            yMin = 1.e-7;
            yMax = 10.0;
        }
    }

    histogram_KineticEnergy_FromSolution->GetXaxis()->SetRangeUser(xMin,xMax);
    histogram_KineticEnergy_FromMembrane->GetXaxis()->SetRangeUser(xMin,xMax);
    histogram_KineticEnergy_FromCytoplasm->GetXaxis()->SetRangeUser(xMin,xMax);

    histogram_KineticEnergy_FromSolution->GetYaxis()->SetRangeUser(yMin,yMax);
    histogram_KineticEnergy_FromMembrane->GetYaxis()->SetRangeUser(yMin,yMax);
    histogram_KineticEnergy_FromCytoplasm->GetYaxis()->SetRangeUser(yMin,yMax);

    c1->SetLogy();
    c1->Update();

    output = "Figures_" + cellGeometry + "/KineticEnergyAlphas_" + cellLine + "_" + std::to_string(activity) + "kBq_logy.pdf";
    c1->SaveAs(output.c_str());
}

void PlotKineticEnergyAlphas()
{

    // MakePlots("C4_2", "D12RP", "Nucleus", 150);
    MakePlots("C4_2", "D5RP", "Nucleus", 150);

    // MakePlots("PC3_PIP", "D5RP", "Nucleus", 150);
    // MakePlots("PC3_PIP", "D12RP", "Nucleus", 150);
}