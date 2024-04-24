#include "../include/SurvivalFit.hpp"
#include "../include/DoseAnalysis.hpp"


DoseAnalysis::DoseAnalysis(SurvivalFit survivalFitInstance)
{
    hDose_Activity_Vec = survivalFitInstance.Get_hDose_Activity_Vec();

    hDose_HitsCellComponent_Vec = survivalFitInstance.Get_hDose_HitsCellComponent_Vec();

    hHitMultiplicity_vec = survivalFitInstance.Get_hHitMultiplicity_vec();

    regionName = survivalFitInstance.Get_RegionName();

    parametersFit = survivalFitInstance.Get_ParametersAndUncertainties_Vec();

    //-----------------------------
        // Define number of bins and max number of hits
    nBinsHits = std::get<1>(hDose_HitsCellComponent_Vec[1])->GetNbinsY();
    nBinsDose = std::get<1>(hDose_HitsCellComponent_Vec[1])->GetNbinsX();
    maxDose = std::get<1>(hDose_HitsCellComponent_Vec[2])->GetXaxis()->GetXmax();

    scalingFactorHistograms = ((double)std::get<1>(hDose_Activity_Vec[1])->GetEntries())/((double)std::get<1>(hDose_Activity_Vec[1])->Integral());
}

void DoseAnalysis::MakeMeanDose_PerNHits_Graphs()
{
    auto FillMeanDose_PerNHits_Graph_OneActivity = [&](double activity, TH2D* hDose_hitsAlpha_CellComponent_OneActivity)
    {
        //--------------------------
        // Making histogram for mean dose per hit
        std::string grMeanDose_PerNHits_OneActivity_Name = "grMeanDose_" + regionName + "_" + std::to_string((int)activity) + "kBq";
        TGraphAsymmErrors* grMeanDose_PerNHits_OneActivity = new TGraphAsymmErrors();
        grMeanDose_PerNHits_OneActivity->SetName(grMeanDose_PerNHits_OneActivity_Name.c_str());
        grMeanDose_PerNHits_OneActivity->GetXaxis()->SetTitle("N Number of Hits to Cell Component by Alpha Particle");
        grMeanDose_PerNHits_OneActivity->GetYaxis()->SetTitle("Mean Dose [Gy]");

        int graphPoint = 0;

        for(int i=0; i<hDose_hitsAlpha_CellComponent_OneActivity->GetNbinsY(); i++)
        {
            //---------------------------
            // Extract dose histogram for i number of hits
            std::string nameProjectedHist = "doseDelivered_NHits_" + std::to_string(i) + "_" + std::to_string((int)activity) + "kBq";
            TH1D* projectedDose_ForNHits = hDose_hitsAlpha_CellComponent_OneActivity->ProjectionX(nameProjectedHist.c_str(),i+1,i+1);

            double intDose_ThisN = projectedDose_ForNHits->Integral();

            if(intDose_ThisN>0)
            {
                int n = 3;
                double x_q[3];
                double y_q[3] = {0.16, 0.50, 0.84};

                projectedDose_ForNHits->GetQuantiles(n, x_q, y_q);

                double medianDose = x_q[1];
                double doseAtQuantileMinus = x_q[0];
                double doseAtQuantilePlus = x_q[2];

                grMeanDose_PerNHits_OneActivity->SetPoint(graphPoint, ((double)i),medianDose);
                grMeanDose_PerNHits_OneActivity->SetPointEYhigh(graphPoint,doseAtQuantilePlus-medianDose);
                grMeanDose_PerNHits_OneActivity->SetPointEYlow(graphPoint,medianDose-doseAtQuantileMinus);

                graphPoint++;
            }
            projectedDose_ForNHits->SetDirectory(0);
        }

        return grMeanDose_PerNHits_OneActivity;
    };

    for(auto & entry : hDose_HitsCellComponent_Vec)
    {
        double activity = std::get<0>(entry);
        TH2D* hDose_hitsAlpha_CellComponent_ThisActivity = std::get<1>(entry);

        std::tuple<double,TGraphAsymmErrors*> tuple = std::make_tuple(activity,FillMeanDose_PerNHits_Graph_OneActivity(activity,hDose_hitsAlpha_CellComponent_ThisActivity));
        grMeanDose_PerNHits_Vec.push_back(tuple);
    }
}

// hDose_HitsCellComponent_Average_Vec

void DoseAnalysis::MakeDose_PerNHits_Average_Histograms()
{
    for(int i=0; i<nBinsHits; i++)
    {
        std::string hDose_ThisNHits_OverlappedAllActivities_Name = "hDose_NHits_" + std::to_string(i) + "_Overlapped_AllActivities";
        TH1D* hDose_ThisNHits_OverlappedAllActivities = new TH1D(hDose_ThisNHits_OverlappedAllActivities_Name.c_str(),"Dose Delivered for one NHit, All Activities Overlapped", nBinsDose, 0., maxDose);
        hDose_ThisNHits_OverlappedAllActivities->GetXaxis()->SetTitle("Dose Delivered [Gy]");
        hDose_ThisNHits_OverlappedAllActivities->GetXaxis()->SetTitle("Fraction of Cells Hit with Dose");

        int numberHistogramsAdded = 0;

        for(auto & entry : hDose_HitsCellComponent_Vec)
        {
            double activity = std::get<0>(entry);
            TH2D* hDose_hitsAlpha_CellComponent_ThisActivity = std::get<1>(entry);

            std::string hDose_ThisNHits_Name = "hDose_ThisNHits";
            TH1D* hDose_ThisNHits = hDose_hitsAlpha_CellComponent_ThisActivity->ProjectionX(hDose_ThisNHits_Name.c_str(),i+1,i+1);

            if(hDose_ThisNHits->Integral()>0.)
            {
                hDose_ThisNHits_OverlappedAllActivities->Add(hDose_ThisNHits);

                numberHistogramsAdded++;
            }
            hDose_ThisNHits->SetDirectory(0);
        }

        hDose_ThisNHits_OverlappedAllActivities->Scale(1./(numberHistogramsAdded));

        hDose_HitsCellComponent_Average_Vec.push_back(std::make_tuple(i,hDose_ThisNHits_OverlappedAllActivities));
        hDose_ThisNHits_OverlappedAllActivities->SetDirectory(0);
    }
}

void DoseAnalysis::MakeMeanDose_PerNHits_Average_Graph()
{
    grMeanDose_PerNHits_Average->SetName("grMeanDose_PerNHits_Average");
    grMeanDose_PerNHits_Average->GetXaxis()->SetTitle("N Number of Hits By Alpha Particle");
    grMeanDose_PerNHits_Average->GetYaxis()->SetTitle("Mean Dose [Gy]");

    int graphPoint = 0;

    for(auto & entry : hDose_HitsCellComponent_Average_Vec)
    {
        int NHits = std::get<0>(entry);
        TH1D* hDose_ThisNHits = std::get<1>(entry);

        if(hDose_ThisNHits->Integral()>0.)
        {
            int n = 3;
            double x_q[3];
            double y_q[3] = {0.16, 0.50, 0.84};

            hDose_ThisNHits->GetQuantiles(n, x_q, y_q);

            double medianDose = x_q[1];
            double doseAtQuantileMinus = x_q[0];
            double doseAtQuantilePlus = x_q[2];

            grMeanDose_PerNHits_Average->SetPoint(graphPoint, ((double)NHits), medianDose);
            grMeanDose_PerNHits_Average->SetPointEYhigh(graphPoint, doseAtQuantilePlus-medianDose);
            grMeanDose_PerNHits_Average->SetPointEYlow(graphPoint, medianDose-doseAtQuantileMinus);

            graphPoint++;
        }
    }
}


void DoseAnalysis::MakeDoseAnalysis()
{
    MakeMeanDose_PerNHits_Graphs();
    MakeDose_PerNHits_Average_Histograms();
    MakeMeanDose_PerNHits_Average_Graph();
}

void DoseAnalysis::WriteToFile(TFile* file)
{
    file->cd();
    for(auto& entry : grMeanDose_PerNHits_Vec)
    {
        (std::get<1>(entry))->Write();
    }

    grMeanDose_PerNHits_Average->Write();
}