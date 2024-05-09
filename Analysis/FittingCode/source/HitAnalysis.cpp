#include "../include/SurvivalFit.hpp"
#include "../include/HitAnalysis.hpp"


HitAnalysis::HitAnalysis(SurvivalFit survivalFitInstance)
{
    hDose_Activity_Vec = survivalFitInstance.Get_hDose_Activity_Vec();

    hDose_HitsCellComponent_Vec = survivalFitInstance.Get_hDose_HitsCellComponent_Vec();

    hHitMultiplicity_vec = survivalFitInstance.Get_hHitMultiplicity_vec();

    regionName = survivalFitInstance.Get_RegionName();

    parametersFit = survivalFitInstance.Get_ParametersAndUncertainties_Vec();

    //-----------------------------
        // Define number of bins and max number of hits
    nBins = std::get<1>(hDose_HitsCellComponent_Vec[1])->GetNbinsY();
    maxDose = std::get<1>(hDose_HitsCellComponent_Vec[2])->GetXaxis()->GetXmax();

    scalingFactorHistograms = ((double)std::get<1>(hDose_Activity_Vec[1])->GetEntries())/((double)std::get<1>(hDose_Activity_Vec[1])->Integral());
}

void HitAnalysis::MakeHitMulitplicityGraphs()
{
    auto MakeHistogramIntoGraph = [&](double activity, TH1D* hHitMultiplicity_OneActivity)
    {
        std::string grHitMultiplicity_OneActivity_Name = "grHitMulitplicity_" + std::to_string((int)activity) + "kBq";
        TGraph* grHitMultiplicity_OneActivity = new TGraphErrors();
        grHitMultiplicity_OneActivity->SetName(grHitMultiplicity_OneActivity_Name.c_str());
        grHitMultiplicity_OneActivity->GetXaxis()->SetTitle("N Number of Hits By Alpha Particle");
        grHitMultiplicity_OneActivity->GetYaxis()->SetTitle("Fraction of Cells Hit N Number of Times");

        int graphPoint = 0;
        //----------------------
        for(int i=0; i<hHitMultiplicity_OneActivity->GetNbinsX(); i++)
        {
            double fraction = hHitMultiplicity_OneActivity->GetBinContent(i+1);

            if(fraction>0.)
            {
                grHitMultiplicity_OneActivity->SetPoint(graphPoint, ((double)i), fraction);

                graphPoint++;
            }
        }

        return grHitMultiplicity_OneActivity;
    };

    for(auto & entry : hHitMultiplicity_vec)
    {
        double activity = std::get<0>(entry);
        TH1D* hHitMuliplicity_ThisActivity = std::get<1>(entry);

        std::tuple<double, TGraph*> tuple = std::make_tuple(activity,MakeHistogramIntoGraph(activity,hHitMuliplicity_ThisActivity));
        grHitMultiplicity_vec.push_back(tuple);
    }
}

void HitAnalysis::MakeHitMultiplicity_SurvivalFraction_Histograms()
{
    double alpha = std::get<0>(parametersFit[0]);
    double dAlpha = std::get<1>(parametersFit[0]);

    //--------------------------------------------
    auto Fill_hHitMultiplicity_SurvivalFraction_OneActivity = [&](double activity, TH2D* hDose_hitsAlpha_CellComponent_OneActivity, TH1D* hDose_mGy_OneActivity, TH1D* hHitMultiplicity_OneActivity)
    {

        //--------------------------
        // Making histogram for cell survival fraction for hit multiplicity
        std::string hHitMultiplicity_SurvivalFraction_ThisActivity_Name = "hHitMultiplicity_" + regionName + "_SurvivalFraction_" + std::to_string((int)activity) + "kBq";
        TH1D* hHitMultiplicity_SurvivalFraction_ThisActivity = new TH1D(hHitMultiplicity_SurvivalFraction_ThisActivity_Name.c_str(), "Cell Survival Fraction per N Number of Hits by Alpha Particle", nBins, 0.,nBins);
        hHitMultiplicity_SurvivalFraction_ThisActivity->GetXaxis()->SetTitle("N Number of Hits to Cell Component by Alpha Particle");
        hHitMultiplicity_SurvivalFraction_ThisActivity->GetYaxis()->SetTitle("Fraction of Total Cells Survived");


        //---------------------------------
        // Function to calculate the error in the survival fraction of one bin at a certain dose
        auto CalculateUncertainty_HitMultiplicity_SurvivalFraction_OneBin = [&](double doseInBin, double fractionHitDoseBin, double fractionSurvivedDoseBin)
        {

            // double a = std::pow((1./fractionHitDoseBin)*dFractionHitDoseBin,2.);
            // double b = std::pow(doseInBin*dAlpha,2.);

            // return fractionSurvivedDoseBin*std::sqrt(a);
            // return doseInBin*dAlpha;

            double dFractionHitDoseBin = (1./scalingFactorHistograms)*std::sqrt(scalingFactorHistograms*fractionHitDoseBin);

            double a = std::pow((1./fractionHitDoseBin)*dFractionHitDoseBin,2.);
            double b = std::pow(doseInBin*dAlpha,2.);

            return fractionSurvivedDoseBin*std::sqrt(a+b);
        };

        //---------------------------
        // Vector to store fraction of cells survived at each number of hits
        // <nHits, survivalFractionNHits, uncertainty survivalFractionNHits>
        std::vector<std::tuple<int,double,double>> hitMultiplicity_SurvivalFraction_Vec;

        //----------------------------
        // Loop over number of hits (y-axis)
        for(int n=0; n<hDose_hitsAlpha_CellComponent_OneActivity->GetNbinsY(); n++)
        {

            // Those not hit always survive
            double survivalFractionThisHitNumber = 0.;

            // Vector to store the uncertainty in the survival fraction at every bin
            std::vector<double> dSurvivaFractionThisDose_Vec;

            double fractionHitNTimes = hHitMultiplicity_OneActivity->GetBinContent(n+1);

            //--------------------------
            // Loop over x axis (dose deposited)
            for(int j=0; j<hDose_hitsAlpha_CellComponent_OneActivity->GetNbinsX(); j++)
            {
                double doseDep = (hDose_hitsAlpha_CellComponent_OneActivity->GetXaxis())->GetBinCenter(j+1);
                double fractionHitThisDose = hDose_hitsAlpha_CellComponent_OneActivity->GetBinContent(j+1,n+1);

                // Only calculate survival if some fraction has dose deposited
                if(fractionHitThisDose>0.)
                {
                    double survivaFractionThisDose = fractionHitThisDose*TMath::Exp(-alpha*doseDep);

                    double dSurvivaFractionThisDose = CalculateUncertainty_HitMultiplicity_SurvivalFraction_OneBin(doseDep, fractionHitThisDose, survivaFractionThisDose);
                    dSurvivaFractionThisDose_Vec.push_back(dSurvivaFractionThisDose);

                    // Adding survival fraction to total survival
                    survivalFractionThisHitNumber += survivaFractionThisDose;
                }
            }

            //------------------------------
            // Calculating uncertainty in survival fraction for hit number
            double dSurvivalFractionThisHitNumber_squared = 0.;
            for(auto& entry : dSurvivaFractionThisDose_Vec)
            {
                dSurvivalFractionThisHitNumber_squared += std::pow(entry,2.);
            }
            double dSurvivalFractionThisHitNumber = std::sqrt(dSurvivalFractionThisHitNumber_squared);
            if(fractionHitNTimes>0.)
            {
                hitMultiplicity_SurvivalFraction_Vec.push_back(std::make_tuple(n,survivalFractionThisHitNumber, dSurvivalFractionThisHitNumber));
            }
        }

        //-------------------------
        // Adding survival fraction to histogram
        for(auto& entry : hitMultiplicity_SurvivalFraction_Vec)
        {
            int nHits = std::get<0>(entry);
            double survivalFraction = std::get<1>(entry);
            double dSurvivalFraction = std::get<2>(entry);

            hHitMultiplicity_SurvivalFraction_ThisActivity->SetBinContent(nHits+1,survivalFraction);
            hHitMultiplicity_SurvivalFraction_ThisActivity->SetBinError(nHits+1,dSurvivalFraction);
        }

        return hHitMultiplicity_SurvivalFraction_ThisActivity;
    };

    for(int i=0; i<hDose_Activity_Vec.size(); i++)
    {
        TH2D* hDose_HitsCellComponent_ThisActivity = std::get<1>(hDose_HitsCellComponent_Vec[i]);
        TH1D* hDose_Activity_ThisActivity = std::get<2>(hDose_Activity_Vec[i]);
        TH1D* hHitMultiplicity_ThisActivity = std::get<1>(hHitMultiplicity_vec[i]);

        double activity = std::get<0>(hDose_Activity_Vec[i]);

        std::tuple<double,TH1D*> tuple = std::make_tuple(activity, Fill_hHitMultiplicity_SurvivalFraction_OneActivity(activity, hDose_HitsCellComponent_ThisActivity, hDose_Activity_ThisActivity,hHitMultiplicity_ThisActivity));
        hHitMultiplicity_SurvivalFraction_vec.push_back(tuple);
    }

}


void HitAnalysis::MakeHitMultiplicity_SurvivalFraction_Graphs()
{
    auto MakeHistogramIntoGraph = [&](double activity, TH1D* hHitMultiplicity_SurvivalFraction_OneActivity)
    {
        std::string grHitMultiplicity_SurvivalFraction_OneActivity_Name = "grHitMulitplicity_SurvivalFraction_" + std::to_string((int)activity) + "kBq";
        TGraphErrors* grHitMultiplicity_SurvivalFraction_OneActivity = new TGraphErrors();
        grHitMultiplicity_SurvivalFraction_OneActivity->SetName(grHitMultiplicity_SurvivalFraction_OneActivity_Name.c_str());
        grHitMultiplicity_SurvivalFraction_OneActivity->GetXaxis()->SetTitle("N Number of Hits By Alpha Particle");
        grHitMultiplicity_SurvivalFraction_OneActivity->GetYaxis()->SetTitle("Fraction of Total Cells Surviving N Hits");

        int graphPoint = 0;

        //----------------------
        for(int i=0; i<hHitMultiplicity_SurvivalFraction_OneActivity->GetNbinsX(); i++)
        {
            double fraction = hHitMultiplicity_SurvivalFraction_OneActivity->GetBinContent(i+1);
            double dFraction = hHitMultiplicity_SurvivalFraction_OneActivity->GetBinError(i+1);

            if(fraction>0.)
            {
                grHitMultiplicity_SurvivalFraction_OneActivity->SetPoint(graphPoint, ((double)i), fraction);
                grHitMultiplicity_SurvivalFraction_OneActivity->SetPointError(graphPoint, 0.0, dFraction);

                graphPoint++;
            }
        }

        return grHitMultiplicity_SurvivalFraction_OneActivity;
    };

    for(auto & entry : hHitMultiplicity_SurvivalFraction_vec)
    {
        double activity = std::get<0>(entry);
        TH1D* hHitMuliplicity_SurvivalFraction_ThisActivity = std::get<1>(entry);

        std::tuple<double, TGraphErrors*> tuple = std::make_tuple(activity,MakeHistogramIntoGraph(activity,hHitMuliplicity_SurvivalFraction_ThisActivity));
        grHitMultiplicity_vec.push_back(tuple);
    }
}

void HitAnalysis::MakeHitMultiplicity_PercentKilled_Histograms()
{


    //-----------------------------
    auto Fill_hHitMultiplicity_PercentKilled_OneActivity = [&](double activity, TH1D* hHitMultiplicity_SurvivalFraction_OneActivity, TH1D* hHitMultiplicity_OneActivity)
    {

        //--------------------------
        // Making histogram for percent of cells killed for hit multiplicity
        std::string hHitMultiplicity_PercentKilled_Name = "hHitMultiplicity_" + regionName + "_PercentKilled_" + std::to_string((int)activity) + "kBq";
        TH1D* hHitMultiplicity_PercentKilled = new TH1D(hHitMultiplicity_PercentKilled_Name.c_str(), "Cell Death per N Number of Hits by Alpha Particle", nBins, 0.,nBins);
        hHitMultiplicity_PercentKilled->GetXaxis()->SetTitle("N Number of Hits to Cell Component by Alpha Particle");
        hHitMultiplicity_PercentKilled->GetYaxis()->SetTitle("Percentage of Cells Killed");

        //----------------------------
        // auto CalculateUncertainty_PercentKilled_NHits = [&](double fractionHitNHits, double dFractionHitNHits, double fractionDeadNHits, double percentDeathNHits, double dFractionSurvivedNHits)
        auto CalculateUncertainty_PercentKilled_NHits = [&](double fractionHitNHits, double fractionHitNHits_survived, double dFractionSurvivedNHits)
        {
            return (100./fractionHitNHits)*dFractionSurvivedNHits;
            // return (100./fractionHitNHits)*fractionHitNHits_survived*dFractionSurvivedNHits;
            // return 0.01*(100.*(fractionHitNHits-fractionHitNHits_survived)/fractionHitNHits);

        };

        //-------------------------
        // Looping over number of hits to cell component
        for(int i=0; i<hHitMultiplicity_SurvivalFraction_OneActivity->GetNbinsX(); i++)
        {
            double fractionHitThisNHits = hHitMultiplicity_OneActivity->GetBinContent(i+1);
            double dFractionHitThisNHits = (1./scalingFactorHistograms)*std::sqrt(scalingFactorHistograms*fractionHitThisNHits);

            if(fractionHitThisNHits>0.)
            {
                double fractionSurvivedThisNHits = hHitMultiplicity_SurvivalFraction_OneActivity->GetBinContent(i+1);
                double dFractionSurvivedThisNHits = hHitMultiplicity_SurvivalFraction_OneActivity->GetBinError(i+1);

                double fractionDeadThisNHit = fractionHitThisNHits - dFractionSurvivedThisNHits;

                double percentKilled = 100.*(fractionHitThisNHits-fractionSurvivedThisNHits)/fractionHitThisNHits;


                double dPercentKilled = CalculateUncertainty_PercentKilled_NHits(fractionHitThisNHits, fractionSurvivedThisNHits, dFractionSurvivedThisNHits);


                hHitMultiplicity_PercentKilled->SetBinContent(i+1, percentKilled);
                hHitMultiplicity_PercentKilled->SetBinError(i+1, dPercentKilled);
            }
        }

        return hHitMultiplicity_PercentKilled;
    };

     for(int i=0; i<hHitMultiplicity_vec.size(); i++)
    {
        double activity = std::get<0>(hHitMultiplicity_vec[i]);
        TH1D* hHitMultiplicity_ThisActivity = std::get<1>(hHitMultiplicity_vec[i]);
        TH1D* hHitMultiplicity_SurvivalFraction_ThisActivity = std::get<1>(hHitMultiplicity_SurvivalFraction_vec[i]);

        std::tuple<double,TH1D*> tuple = std::make_tuple(activity, Fill_hHitMultiplicity_PercentKilled_OneActivity(activity, hHitMultiplicity_SurvivalFraction_ThisActivity, hHitMultiplicity_ThisActivity));
        hHitMultiplicity_PercentKilled_Vec.push_back(tuple);
    }
}


void HitAnalysis::MakeHitMultiplicity_PercentKilled_Graphs()
{
    auto MakeHistogramIntoGraph = [&](double activity, TH1D* hHitMultiplicity_PercentKilled_OneActivity)
    {
        std::string grHitMultiplicity_PercentKilled_OneActivity_Name = "grHitMulitplicity_PercentKilled_" + std::to_string((int)activity) + "kBq";
        TGraphErrors* grHitMultiplicity_PercentKilled_OneActivity = new TGraphErrors();
        grHitMultiplicity_PercentKilled_OneActivity->SetName(grHitMultiplicity_PercentKilled_OneActivity_Name.c_str());
        grHitMultiplicity_PercentKilled_OneActivity->GetXaxis()->SetTitle("N Number of Hits By Alpha Particle");
        grHitMultiplicity_PercentKilled_OneActivity->GetYaxis()->SetTitle("Probability Death for N Number of Hits");

        int graphPoint = 0;
        //----------------------
        for(int i=0; i<hHitMultiplicity_PercentKilled_OneActivity->GetNbinsX(); i++)
        {
            double percent = hHitMultiplicity_PercentKilled_OneActivity->GetBinContent(i+1);
            double dPercent = hHitMultiplicity_PercentKilled_OneActivity->GetBinError(i+1);

            if(percent>0.)
            {
                grHitMultiplicity_PercentKilled_OneActivity->SetPoint(graphPoint, ((double)i), percent);
                grHitMultiplicity_PercentKilled_OneActivity->SetPointError(graphPoint, 0.0, dPercent);

                graphPoint++;
            }
        }

        return grHitMultiplicity_PercentKilled_OneActivity;
    };

    for(auto & entry : hHitMultiplicity_PercentKilled_Vec)
    {
        double activity = std::get<0>(entry);
        TH1D* hHitMuliplicity_PercentKilled_ThisActivity = std::get<1>(entry);

        std::tuple<double, TGraphErrors*> tuple = std::make_tuple(activity,MakeHistogramIntoGraph(activity,hHitMuliplicity_PercentKilled_ThisActivity));
        grHitMultiplicity_PercentKilled_Vec.push_back(tuple);
    }
}

void HitAnalysis::MakeGraph_FractionHit_PerNHits(int MaxNHits)
{
    auto Fill_grFractionHit_PerNHits = [&](int NHits)
    {
        std::string grFractionHit_PerOneHit_Name = "grFractionHit_" + std::to_string(NHits) + "_NumberHitsToCellComponent";
        TGraphErrors* grFractionHit_PerOneHit = new TGraphErrors();
        grFractionHit_PerOneHit->SetName(grFractionHit_PerOneHit_Name.c_str());
        grFractionHit_PerOneHit->GetXaxis()->SetTitle("Activity [kBq/mL");
        grFractionHit_PerOneHit->GetXaxis()->SetTitle("Fraction of cells in sample hit");

        std::vector<std::tuple<double,double>> fractionHit_Activity_ThisNHit_Vec;

        for(int i=0; i<hHitMultiplicity_vec.size(); i++)
        {
            double activity = std::get<0>(hHitMultiplicity_vec[i]);
            TH1D* hHitMultiplicity_OneActivity = std::get<1>(hHitMultiplicity_vec[i]);

            if(hHitMultiplicity_OneActivity->GetNbinsX()>0)
            {
                double fractionHit = hHitMultiplicity_OneActivity->GetBinContent(NHits+1);

                if(fractionHit>0.)
                {
                    fractionHit_Activity_ThisNHit_Vec.push_back(std::make_tuple(activity,fractionHit));
                }
            }
        }
        for(int i=0; i<fractionHit_Activity_ThisNHit_Vec.size(); i++)
        {
            double activity = std::get<0>(fractionHit_Activity_ThisNHit_Vec[i]);
            double fractionHit_NHits = std::get<1>(fractionHit_Activity_ThisNHit_Vec[i]);

            grFractionHit_PerOneHit->SetPoint(i, activity, fractionHit_NHits);
        }

        return grFractionHit_PerOneHit;
    };

    for(int i=0; i<MaxNHits; i++)
    {
        TGraphErrors* gr = Fill_grFractionHit_PerNHits(i);
        if(gr->GetN()>0)
        {
            std::tuple<int,TGraphErrors*> tuple = std::make_tuple(i,gr);
            grFractionHit_PerOneHit_Vec.push_back(tuple);
        }
    }
}


void HitAnalysis::MakeGraph_ProbabilityDeath_PerNHits(int MaxNHits)
{
    //---------------------------
    auto Fill_grPercentKilled_PerOneNHit = [&](int NHits)
    {
        //----------------------------
        std::string grPercentKilled_PerOneNHit_Name = "grProbabilityDeath_" + std::to_string(NHits) + "_NumberHitsToCellComponent";
        TGraphErrors* grPercentKilled_PerOneNHit = new TGraphErrors();
        grPercentKilled_PerOneNHit->SetName(grPercentKilled_PerOneNHit_Name.c_str());
        grPercentKilled_PerOneNHit->GetXaxis()->SetTitle("Activity [kBq/mL]");
        grPercentKilled_PerOneNHit->GetYaxis()->SetTitle("Probability of Cell Death");

        std::vector<std::tuple<double,double,double>> percentKilled_Activity_ThisNHit_Vec;

        //--------------------------------
        for(int i=0; i<hHitMultiplicity_PercentKilled_Vec.size(); i++)
        {
            double activity = std::get<0>(hHitMultiplicity_PercentKilled_Vec[i]);
            TH1D* hHitMultiplicity_PercentKilled_OneActivity = std::get<1>(hHitMultiplicity_PercentKilled_Vec[i]);

            if(hHitMultiplicity_PercentKilled_OneActivity->GetNbinsX()>0)
            {
                double percentKilled_NHits = hHitMultiplicity_PercentKilled_OneActivity->GetBinContent(NHits+1);
                double dPercentKilled_NHits = hHitMultiplicity_PercentKilled_OneActivity->GetBinError(NHits+1);

                if(percentKilled_NHits>0.){percentKilled_Activity_ThisNHit_Vec.push_back(std::make_tuple(activity,percentKilled_NHits,dPercentKilled_NHits));}
            }
        }

        for(int i=0; i<percentKilled_Activity_ThisNHit_Vec.size(); i++)
        {
            double activity = std::get<0>(percentKilled_Activity_ThisNHit_Vec[i]);
            double percentKilled_NHits = std::get<1>(percentKilled_Activity_ThisNHit_Vec[i]);
            double dPercentKilled_NHits = std::get<2>(percentKilled_Activity_ThisNHit_Vec[i]);

            grPercentKilled_PerOneNHit->SetPoint(i, activity, percentKilled_NHits);
            grPercentKilled_PerOneNHit->SetPointError(i, 0., dPercentKilled_NHits);
        }

        return grPercentKilled_PerOneNHit;
    };

    for(int i=0; i<MaxNHits; i++)
    {
        TGraphErrors* gr = Fill_grPercentKilled_PerOneNHit(i);
        if(gr->GetN()>0)
        {
            std::tuple<int,TGraphErrors*> tuple = std::make_tuple(i,gr);
            grPercentKilled_PerOneHit_Vec.push_back(tuple);
        }
    }
}


void HitAnalysis::MakeGraph_ProbabilityDeath_UWA_ForNHits()
{
    //----------------------------
    grProbabilityDeath_ForNHits->SetName("grProbabilityDeath_UWA_ForNHits");
    grProbabilityDeath_ForNHits->GetXaxis()->SetTitle("Number of Hits by Alpha Particle");
    grProbabilityDeath_ForNHits->GetYaxis()->SetTitle("Probabilitry of Death at N Number of Hits");

    int NPoints = 0;
    for(auto & entry : grPercentKilled_PerOneHit_Vec)
    {
        int NHits = std::get<0>(entry);
        TGraphErrors* grProbDeath = std::get<1>(entry);

        //--------------------
        double sum_w_i = 0.;
        double sum_mu_i_w = 0.;

        for(int n=0; n<grProbDeath->GetN(); n++)
        {
            double probabilityDeath = grProbDeath->GetPointY(n);
            double dProbabilityDeath = grProbDeath->GetErrorY(n);

            double w_i = 1./std::pow(dProbabilityDeath,2.);
            double mu_i_w = w_i*probabilityDeath;

            sum_w_i += w_i;
            sum_mu_i_w += mu_i_w;
        }

        double UWA_percentDeath_perN = sum_mu_i_w/sum_w_i;
        double dUWA_percentDeath_perN = 1./std::sqrt(sum_w_i);

        grProbabilityDeath_ForNHits->SetPoint(NPoints, ((double)NHits), UWA_percentDeath_perN);
        grProbabilityDeath_ForNHits->SetPointError(NPoints, 0., dUWA_percentDeath_perN);

        NPoints++;
    }
}

void HitAnalysis::MakeHitAnalysis(int MaxNHits)
{
    MakeHitMulitplicityGraphs();
    MakeHitMultiplicity_SurvivalFraction_Histograms();
    MakeHitMultiplicity_SurvivalFraction_Graphs();
    MakeHitMultiplicity_PercentKilled_Histograms();
    MakeHitMultiplicity_PercentKilled_Graphs();
    MakeGraph_FractionHit_PerNHits(MaxNHits);
    MakeGraph_ProbabilityDeath_PerNHits(MaxNHits);
    MakeGraph_ProbabilityDeath_UWA_ForNHits();
}

void HitAnalysis::WriteToFile(TFile* file)
{
    file->cd();
    for(auto& entry : grHitMultiplicity_vec)
    {
        (std::get<1>(entry))->Write();
    }

    for(auto& entry : grHitMultiplicity_SurvivalFraction_vec)
    {
        (std::get<1>(entry))->Write();
    }

    for(auto& entry : grHitMultiplicity_PercentKilled_Vec)
    {
        (std::get<1>(entry))->Write();
    }

    for(auto& entry : grPercentKilled_PerOneHit_Vec)
    {
        (std::get<1>(entry))->Write();
    }

    for(auto & entry : grFractionHit_PerOneHit_Vec)
    {
        (std::get<1>(entry))->Write();
    }

    grProbabilityDeath_ForNHits->Write();
}