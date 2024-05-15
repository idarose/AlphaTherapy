#include "../include/SurvivalFit.hpp"
#include "../include/CellSurvival.hpp"
#include <iostream>
#include <fstream>

SurvivalFit::SurvivalFit() {
}

//----------------------------
void SurvivalFit::LoadHistogramsFromAnalysis(CellSurvival cellSurvivalInstance, std::string regionName)
{

    //--------------------------
    // Get cell survival data
    cellSurvivalData = cellSurvivalInstance.GetCellSurvivalData();

    //-----------------------------
    // General filepath to output
    std::string generalFilePath = "../AnalysisCode/Output_" + cellSurvivalInstance.GetCellGeometryType() + "/Output_" + cellSurvivalInstance.GetCellLine() + "_";

    // Specific filepath
    std::string filePath;

    //---------------------------
    // Making zero activity case by cloning the first histogram and resetting it
    double firstActivity = std::get<0>(cellSurvivalData[0]);
    filePath = generalFilePath + std::to_string(((int)firstActivity)) + "kBq.root";

    //-------------------------
    // // Making specific filepath
    // double firstActivity = std::get<0>(cellSurvivalData[0]);
    // filePath = generalFilePath + std::to_string(((int)firstActivity)) + "kBq_case_3.root";

    // Read file containing histogram
    TFile* inputFile = new TFile(filePath.c_str(), "READ");

    //-------------------------------------
    // Extracting dose deposition histogram for cell component with mGY binning
    TH1D* hDose_CellComponent_mGy = nullptr;
    std::string histogramName_mGy = "i0_hDose_212Pb_" + cellSurvivalInstance.GetCellLine() + "_" + std::to_string(((int)firstActivity)) +"kBq_" + regionName + "_mGyBinning";
    std::string histogramNewName = "i0_hDose_212Pb_" + cellSurvivalInstance.GetCellLine() + "_" + std::to_string(0) +"kBq_" + regionName + "_mGyBinning";


    inputFile->GetObject(histogramName_mGy.c_str(), hDose_CellComponent_mGy);
    hDose_CellComponent_mGy->SetDirectory(0);
    hDose_CellComponent_mGy->SetName(histogramNewName.c_str());
    hDose_CellComponent_mGy->Reset();


    //-------------------------------------
    // Extracting dose deposition histogram for nucleus with uGyBinning
    TH1D* hDose_CellComponent_uGy = nullptr;
    std::string histogramName_uGy = "i0_hDose_212Pb_" + cellSurvivalInstance.GetCellLine() + "_" + std::to_string(((int)firstActivity)) +"kBq_" + regionName + "_uGyBinning";
    histogramNewName = "i0_hDose_212Pb_" + cellSurvivalInstance.GetCellLine() + "_" + std::to_string(0) +"kBq_" + regionName + "_uGyBinning";

    inputFile->GetObject(histogramName_uGy.c_str(), hDose_CellComponent_uGy);
    hDose_CellComponent_uGy->SetDirectory(0);
    hDose_CellComponent_uGy->SetName(histogramNewName.c_str());
    hDose_CellComponent_uGy->Reset();


    //-------------------------------------
    // Extracting 2D histogram for dose in cell component per number of alpha-particle hits
    TH2D* hDose_hitsAlpha_CellComponent = nullptr;
    std::string histogramName_Dose_HitsAlpha_CellComponent = "i0_hDose_212Pb_" + cellSurvivalInstance.GetCellLine() + "_" + std::to_string(((int)firstActivity)) +"kBq_HitsAlpha_" + regionName;
    histogramNewName = "i0_hDose_212Pb_" + cellSurvivalInstance.GetCellLine() + "_" + std::to_string(0) +"kBq_HitsAlpha_" + regionName;


    inputFile->GetObject(histogramName_Dose_HitsAlpha_CellComponent.c_str(), hDose_hitsAlpha_CellComponent);
    hDose_hitsAlpha_CellComponent->SetDirectory(0);
    hDose_hitsAlpha_CellComponent->SetName(histogramNewName.c_str());
    hDose_hitsAlpha_CellComponent->Reset();


    //-------------------------------------
    // Extracting hit multiplicity histogram for cell component
    TH1D* hHitMultiplicity_CellComponent = nullptr;
    std::string histogramName_hitMultiplicity = "i0_hFractionHitsAlpha_212Pb_" + cellSurvivalInstance.GetCellLine() + "_" + std::to_string(((int)firstActivity)) +"kBq_" + regionName;
    histogramNewName = "i0_hFractionHitsAlpha_212Pb_" + cellSurvivalInstance.GetCellLine() + "_" + std::to_string(0) +"kBq_" + regionName;


    inputFile->GetObject(histogramName_hitMultiplicity.c_str(), hHitMultiplicity_CellComponent);
    hHitMultiplicity_CellComponent->SetDirectory(0);
    hHitMultiplicity_CellComponent->SetName(histogramNewName.c_str());
    hHitMultiplicity_CellComponent->Reset();


    //-------------------------------------
    inputFile->Close();


    hDose_Activity_Vec.push_back(std::make_tuple(0.,hDose_CellComponent_uGy,hDose_CellComponent_mGy));
    hHitMultiplicity_vec.push_back(std::make_tuple(0,hHitMultiplicity_CellComponent));
    hDose_HitsCellComponent_Vec.push_back(std::make_tuple(0.,hDose_hitsAlpha_CellComponent));


    //---------------------------------------
    std::vector<double> activities = {1.,3.,5.,10.,25.,50.,75.,100.,150.};

    for(int i=0; i<activities.size(); i++)
    {
        // Extracting activity
        double activity = activities[i];

        // Making specific filepath
        filePath = generalFilePath + std::to_string(((int)activity)) + "kBq.root";

        // if(cellLine=="C4_2"||cellLine=="PC3_Flu")
        // {
        //     filePath = generalFilePath + std::to_string(((int)activity)) + "kBq_case_1.root";
        // }
        // else
        // {
        //     filePath = generalFilePath + std::to_string(((int)activity)) + "kBq_case_3.root";
        // }

        // Read file containing histogram
        TFile* inputFile = new TFile(filePath.c_str(), "READ");


        //-------------------------------------
        // Extracting dose deposition histogram for nucleus with keV binning
        TH1D* hDose_CellComponent_mGy = nullptr;
        std::string histogramName_mGy = "i0_hDose_212Pb_" + cellSurvivalInstance.GetCellLine() + "_" + std::to_string(((int)activity)) +"kBq_" + regionName + "_mGyBinning";

        inputFile->GetObject(histogramName_mGy.c_str(), hDose_CellComponent_mGy);
        hDose_CellComponent_mGy->SetDirectory(0);


        //-------------------------------------
        // Extracting dose deposition histogram for nucleus with eV binning
        TH1D* hDose_CellComponent_uGy = nullptr;
        std::string histogramName_uGy = "i0_hDose_212Pb_" + cellSurvivalInstance.GetCellLine() + "_" + std::to_string(((int)activity)) +"kBq_" + regionName + "_uGyBinning";

        inputFile->GetObject(histogramName_uGy.c_str(), hDose_CellComponent_uGy);
        hDose_CellComponent_uGy->SetDirectory(0);

        //-------------------------------------
        // Extracting 2D histogram for dose in cell component per number of alpha-particle hits
        TH2D* hDose_hitsAlpha_CellComponent = nullptr;
        std::string histogramName_Dose_HitsAlpha_CellComponent = "i0_hDose_212Pb_" + cellSurvivalInstance.GetCellLine() + "_" + std::to_string(((int)activity)) +"kBq_HitsAlpha_" + regionName;

        inputFile->GetObject(histogramName_Dose_HitsAlpha_CellComponent.c_str(), hDose_hitsAlpha_CellComponent);
        hDose_hitsAlpha_CellComponent->SetDirectory(0);

        //-------------------------------------
        // Extracting hit multiplicity histogram for cell component
        TH1D* hHitMultiplicity_CellComponent = nullptr;
        std::string histogramName_hitMultiplicity = "i0_hFractionHitsAlpha_212Pb_" + cellSurvivalInstance.GetCellLine() + "_" + std::to_string(((int)activity)) +"kBq_" + regionName;

        inputFile->GetObject(histogramName_hitMultiplicity.c_str(), hHitMultiplicity_CellComponent);
        hHitMultiplicity_CellComponent->SetDirectory(0);


        //-------------------------------------
        inputFile->Close();


        hDose_Activity_Vec.push_back(std::make_tuple(activity,hDose_CellComponent_uGy,hDose_CellComponent_mGy));
        hHitMultiplicity_vec.push_back(std::make_tuple(activity,hHitMultiplicity_CellComponent));
        hDose_HitsCellComponent_Vec.push_back(std::make_tuple(activity,hDose_hitsAlpha_CellComponent));
    }

    for(auto & entry : hDose_Activity_Vec)
    {
        std::cout << "Activity : " << std::get<2>(entry)->Integral() << std::endl;
    }
}



//----------------------------
void SurvivalFit::FitCellSurvival(CellSurvival cellSurvivalInstance, std::string modelName_in, int region)
{
    //--------------------------------
    // Function fitting the cell survival data provided with a model, either:
    // LM : Linear model, beta set to zero automatically
    // LQ : Linear-quadratic model
    // The region is
    // Returns a vector of tuples containing the fit parameters and their uncertainties

    modelName = modelName_in;

    cellLine = cellSurvivalInstance.GetCellLine();

    cellGeometryType = cellSurvivalInstance.GetCellGeometryType();

    //------------------------
    // Defininng volume used for fit
    if(region==1)
    {
        // Fitting for membrane
        regionName = "Membrane";
    }
    else if(region==2)
    {
        //Fitting for cytoplasm
        regionName = "Cytoplasm";
    }
    else if(region==3)
    {
        // Fitting for nucleus
        regionName = "Nucleus";
    }
    else if(region==4)
    {
        regionName = "TotalCell";
    }

    //----------------------------
    // Load histograms
    LoadHistogramsFromAnalysis(cellSurvivalInstance,regionName);


    //--------------------------
    // Getting and graphing cell survival data
    std::vector<std::tuple<double,double,double>> data_cellSurvival;
    data_cellSurvival = cellSurvivalInstance.GetCellSurvivalData();


    //-----------------------------------

    gr_clonogenicSurvival = new TGraphErrors();
    gr_clonogenicSurvival->SetName("gr_clonogenicSurvival");

    auto GraphData = [&](TGraphErrors* gr_clonogenicSurvival_CellLine, std::vector<std::tuple<double, double, double>> data_cellSurvival_CellLine)
    {
        int graphPointN_data = 0;
        for(auto & entry : data_cellSurvival_CellLine)
        {
            double activity = std::get<0>(entry);
            double cellSurvival = std::get<1>(entry);
            double cellSurvivalUncertainty = std::get<2>(entry);

            gr_clonogenicSurvival_CellLine->SetPoint(graphPointN_data, activity, cellSurvival);
            gr_clonogenicSurvival_CellLine->SetPointError(graphPointN_data, 0.0, cellSurvivalUncertainty);

            graphPointN_data++;
        }
    };

    GraphData(gr_clonogenicSurvival, data_cellSurvival);

    //------------------------------
    // Linear quadratic model
    auto CalculateCellSurvivalFraction_LQModel = [&](TH1D *hDose_uGy, TH1D *hDose_mGy, double *par)
    {
        double alpha = par[0];
        double beta = par[1];

        //------------------------------------------------
        double fractionOfComponentsHit = hDose_mGy->Integral();
        double fractionOfComponentsMissed = (1.0 - fractionOfComponentsHit);

        //------------------------------------------------
        //      The fraction of missed cells obviously all survive so that is immediately added.
        double fractionOfTotalCellsHitSurviving = fractionOfComponentsMissed;

        // Looping over uGy binned histogram
        for(int i=0; i<hDose_uGy->GetNbinsX();i++)
        {
            double dose = hDose_uGy->GetBinCenter(i+1);
            double fractionOfTotalCellsHit = hDose_uGy->GetBinContent(i+1);

            //----------------------------
            double cellSurvivalfraction = TMath::Exp(-(alpha*dose + beta*TMath::Power(dose, 2.0)));
            double fractionOfTotalCellsHit_survivingFraction = cellSurvivalfraction*fractionOfTotalCellsHit;

            fractionOfTotalCellsHitSurviving += fractionOfTotalCellsHit_survivingFraction;
        }

        // Looping over mGy binned histogram
        for(int i=1000; i<hDose_mGy->GetNbinsX(); i++)
        {
            double dose = hDose_mGy->GetBinCenter(i+1);
            double fractionOfTotalCellsHit = hDose_mGy->GetBinContent(i+1);

            //----------------------------
            double cellSurvivalfraction = TMath::Exp(-(alpha*dose + beta*TMath::Power(dose, 2.0)));
            double fractionOfTotalCellsHit_survivingFraction = cellSurvivalfraction*fractionOfTotalCellsHit;

            fractionOfTotalCellsHitSurviving += fractionOfTotalCellsHit_survivingFraction;
        }


        return fractionOfTotalCellsHitSurviving;
    };

    //------------------------------
    // Linear model
    auto CalculateCellSurvivalFraction_LModel = [&](TH1D *hDose_uGy, TH1D *hDose_mGy, double *par)
    {
        double alpha = par[0];

        //------------------------------------------------
        double fractionOfComponentsHit = hDose_mGy->Integral();
        double fractionOfComponentsMissed = (1.0 - fractionOfComponentsHit);

        //------------------------------------------------
        //      The fraction of missed cells obviously all survive so that is immediately added.
        double fractionOfTotalCellsHitSurviving = fractionOfComponentsMissed;

        // Looping over uGy binned histogram
        for(int i=0; i<hDose_uGy->GetNbinsX();i++)
        {
            double dose = hDose_uGy->GetBinCenter(i+1);
            double fractionOfTotalCellsHit = hDose_uGy->GetBinContent(i+1);

            //----------------------------
            double cellSurvivalfraction = TMath::Exp(-(alpha*dose));
            double fractionOfTotalCellsHit_survivingFraction = cellSurvivalfraction*fractionOfTotalCellsHit;


            fractionOfTotalCellsHitSurviving += fractionOfTotalCellsHit_survivingFraction;
        }

        // Looping over mGy binned histogram
        for(int i=1000; i<hDose_mGy->GetNbinsX(); i++)
        {
            double dose = hDose_mGy->GetBinCenter(i+1);
            double fractionOfTotalCellsHit = hDose_mGy->GetBinContent(i+1);

            //----------------------------
            double cellSurvivalfraction = TMath::Exp(-(alpha*dose));
            double fractionOfTotalCellsHit_survivingFraction = cellSurvivalfraction*fractionOfTotalCellsHit;

            fractionOfTotalCellsHitSurviving += fractionOfTotalCellsHit_survivingFraction;
        }

        //------------------------------------------------
        return fractionOfTotalCellsHitSurviving;
    };

    //------------------------------------------------------------------------------------------------
    //      This generates a new graph corresponding to the predicted survivability.
    //      This graph is calculated everytime the parameters (par) are changed.
    auto GenerateGraph_CellSurvivalFraction = [&](std::vector<std::tuple<double,TH1D*,TH1D*>> vec, double *par)
    {
        TGraphErrors gr;

        int graphPointN = 0;

        for(auto & entry : vec)
        {
            double activity = std::get<0>(entry);
            auto histDose_uGy = std::get<1>(entry);
            auto histDose_mGy = std::get<2>(entry);
            double cellSurvival;

            if(modelName=="LQ")
            {
                cellSurvival = CalculateCellSurvivalFraction_LQModel(histDose_uGy, histDose_mGy, par);
            }

            if(modelName=="LM")
            {
                cellSurvival = CalculateCellSurvivalFraction_LModel(histDose_uGy, histDose_mGy, par);
            }

            gr.SetPoint(graphPointN, activity, cellSurvival);

            graphPointN++;
        }

        return gr;
    };

     //------------------------------------------------------------------------------------------------
    int nParameters;

    if(modelName=="LQ")
    {
        nParameters = 2;
    }
    else if(modelName=="LM")
    {
        nParameters = 1;
    }

    double savedParameters[nParameters];

    //------------------------------------------------------------------------------------------------
    f_cellSurvivalVsDose = new TF1("f_cellSurvivalVsDose",
        [&](double*x, double *p)
        {
            double activity_Bq = x[0];

            //------------------------------------------------
            //      This block of code is to test whether the parameters have changed.
            //      A new graph of survivability is generated only if the parameters are found to have changed.
            //      This makes the code much faster, which is useful if your calculation for cell survivability becomes more computationally heavy
            bool foundChangeInParameters = false;

            for(int i=0; i<nParameters; i++)
            {
                if(savedParameters[i]!=p[i])
                {
                    foundChangeInParameters = true;
                    savedParameters[i] = p[i];
                }
            }

            if(foundChangeInParameters)
            {
                gr_cellSurvivability_vs_activitykBqPerMl = GenerateGraph_CellSurvivalFraction(hDose_Activity_Vec, p);
            }

            //------------------------------------------------
            double cellSurvival = gr_cellSurvivability_vs_activitykBqPerMl.Eval(x[0], 0, "");

            return cellSurvival;

        }, 0.0, 150.0, nParameters);


    //------------------------
    f_cellSurvivalVsDose->SetNpx(10000);

    if(modelName=="LM")
    {
        f_cellSurvivalVsDose->SetParLimits(0, 0.0, 1.0e+05);

        f_cellSurvivalVsDose->SetParameter(0, 1.0e+00);
    }
    if(modelName=="LQ")
    {
        f_cellSurvivalVsDose->SetParLimits(0, 0.0, 1.0e+05);
        f_cellSurvivalVsDose->SetParLimits(1, 0.0, 1.0e+05);

        f_cellSurvivalVsDose->SetParameter(0, 1.0e+00);
        f_cellSurvivalVsDose->SetParameter(1, 1.0e+00);
    }


    //-----------------------------
    std::cout << " Cell Line " + cellSurvivalInstance.GetCellLine() << " Model : " << modelName << " Region : " << regionName << std::endl;
    gr_clonogenicSurvival->Fit(f_cellSurvivalVsDose, "", "", 0.0, 150.0);


    double chi_sq = f_cellSurvivalVsDose->GetChisquare();

    int deg_freedom = f_cellSurvivalVsDose->GetNDF();

    std::cout << "Chi squared: " << chi_sq << std::endl;
    std::cout << "Deg Freedom: " << deg_freedom << std::endl;
    std::cout << "Reduced Chi squared: " << chi_sq/deg_freedom << std::endl;


    //-----------------------
    // Extracting fit parameters

    double alpha = f_cellSurvivalVsDose->GetParameter(0);
    double dAlpha = f_cellSurvivalVsDose->GetParError(0);

    // std::cout << alpha << std::endl;


    parametersAndUncertainties_Vec.push_back(std::make_tuple(alpha,dAlpha));

    double beta = 0.;
    double dBeta = 0.;

    if(modelName=="LQ")
    {
        beta = f_cellSurvivalVsDose->GetParameter(1);
        dBeta = f_cellSurvivalVsDose->GetParError(1);
    }


    //---------------------------
    auto CalculateUncertainty_HitMultiplicity_SurvivalFraction_OneBin = [&](double doseInBin, double fractionHitDoseBin, double fractionSurvivedDoseBin, double scalingFactor)
    {
        double dFractionHitDoseBin = (1./scalingFactor)*std::sqrt(scalingFactor*fractionHitDoseBin);

        double a = std::pow((1./fractionHitDoseBin)*dFractionHitDoseBin,2.);
        double b = std::pow(doseInBin*dAlpha,2.);
        double c = std::pow(std::pow(doseInBin,2.)*dBeta,2.);

        // std::cout << " a : " << a << " b : " << b << std::endl;
        // std::cout << scalingFactor*fractionHitDoseBin << std::endl;

        return fractionSurvivedDoseBin*std::sqrt(a+b+c);
    };

    //------------------------------
    // Linear quadratic model
    auto CalculateCellSurvivalFraction_final = [&](TH1D *hDose_uGy, TH1D *hDose_mGy)
    {
        std::vector<double> uncertaintyInSurvOneDose_Vec;

        //------------------------------------------------
        double fractionOfComponentsHit = hDose_mGy->Integral();
        double fractionOfComponentsMissed = (1.0 - fractionOfComponentsHit);

        //------------------------------------------------
        //      The fraction of missed cells obviously all survive so that is immediately added.
        double fractionOfTotalCellsHitSurviving = fractionOfComponentsMissed;

        double scalingFactorHistogram = ((double)hDose_uGy->GetEntries())/((double)hDose_uGy->Integral());

        // Looping over uGy binned histogram
        for(int i=0; i<hDose_uGy->GetNbinsX();i++)
        {
            double dose = hDose_uGy->GetBinCenter(i+1);
            double fractionOfTotalCellsHit = hDose_uGy->GetBinContent(i+1);

            if(fractionOfTotalCellsHit>0.)
            {
                double cellSurvivalfraction = TMath::Exp(-alpha*dose - beta*std::pow(dose,2.));
                double fractionOfTotalCellsHit_survivingFraction = cellSurvivalfraction*fractionOfTotalCellsHit;

                double d_fractionOfTotalCellsHit_survivingFraction = CalculateUncertainty_HitMultiplicity_SurvivalFraction_OneBin(dose,fractionOfTotalCellsHit,fractionOfTotalCellsHit_survivingFraction, scalingFactorHistogram);

                uncertaintyInSurvOneDose_Vec.push_back(d_fractionOfTotalCellsHit_survivingFraction);

                fractionOfTotalCellsHitSurviving += fractionOfTotalCellsHit_survivingFraction;
            }
        }

        scalingFactorHistogram = ((double)hDose_mGy->GetEntries())/((double)hDose_mGy->Integral());

        // Looping over mGy binned histogram
        for(int i=1000; i<hDose_mGy->GetNbinsX(); i++)
        {
            double dose = hDose_mGy->GetBinCenter(i+1);
            double fractionOfTotalCellsHit = hDose_mGy->GetBinContent(i+1);

            if(fractionOfTotalCellsHit>0.)
            {
                //----------------------------
                double cellSurvivalfraction = TMath::Exp(-alpha*dose - beta*std::pow(dose,2.));
                double fractionOfTotalCellsHit_survivingFraction = cellSurvivalfraction*fractionOfTotalCellsHit;

                double d_fractionOfTotalCellsHit_survivingFraction = CalculateUncertainty_HitMultiplicity_SurvivalFraction_OneBin(dose,fractionOfTotalCellsHit,fractionOfTotalCellsHit_survivingFraction, scalingFactorHistogram);

                uncertaintyInSurvOneDose_Vec.push_back(d_fractionOfTotalCellsHit_survivingFraction);

                fractionOfTotalCellsHitSurviving += fractionOfTotalCellsHit_survivingFraction;
            }
        }

        double d_fractionOfTotalCellsHitSurviving = 0.;
        for(auto & entry : uncertaintyInSurvOneDose_Vec)
        {
            d_fractionOfTotalCellsHitSurviving += std::pow(entry,2.);

        }
        d_fractionOfTotalCellsHitSurviving = std::pow(d_fractionOfTotalCellsHitSurviving,1./2.);

        //------------------------------------------------
        return std::make_tuple(fractionOfTotalCellsHitSurviving,d_fractionOfTotalCellsHitSurviving);
    };

    grFitPoints->SetName("grFitPoints");
    grFitPoints->SetTitle("Cell Survival Fraction");
    grFitPoints->GetXaxis()->SetTitle("Activity [kBq/mL]");
    grFitPoints->GetYaxis()->SetTitle("Fraction of cells in sample surviving");

    int graphPoint = 0;
    for(int i=0; i<hDose_Activity_Vec.size(); i++)
    {
        double activity = std::get<0>(hDose_Activity_Vec[i]);

        if(activity>0.)
        {
            // double surv = f_cellSurvivalVsDose->Eval(activity);
            TH1D* hDose_uGy_ThisActivity = std::get<1>(hDose_Activity_Vec[i]);
            TH1D* hDose_mGy_ThisActivity = std::get<2>(hDose_Activity_Vec[i]);
            double surv = std::get<0>(CalculateCellSurvivalFraction_final(hDose_uGy_ThisActivity,hDose_mGy_ThisActivity));
            double dSurv = std::get<1>(CalculateCellSurvivalFraction_final(hDose_uGy_ThisActivity,hDose_mGy_ThisActivity));

            // std::cout << activity << ": " << dSurv << std::endl;

            grFitPoints->SetPoint(graphPoint, activity, surv);
            grFitPoints->SetPointError(graphPoint, 0., dSurv);

            graphPoint++;
        }
    }

}



//----------------------------
void SurvivalFit::WriteToFile(TFile *file)
{
    // std::string outputName = "Output_AnalyseCellSurvival_" + cellLine + "_" + regionName + "_" + modelName + ".root";
    // TFile *outputFile = new TFile(outputName.c_str(), "RECREATE");

    file->cd();

    std::string titleGraph = "Cell Survival Using Energy Deposition in " + regionName;
    gr_clonogenicSurvival->SetTitle(titleGraph.c_str());
    gr_clonogenicSurvival->GetXaxis()->SetTitle("Activity [kBq/mL]");
    gr_clonogenicSurvival->GetYaxis()->SetTitle("Survival Fraction");
    gr_clonogenicSurvival->Write();


    f_cellSurvivalVsDose->Write();
    gr_cellSurvivability_vs_activitykBqPerMl.Write();

    grFitPoints->Write();



    // outputFile->Write();
    // outputFile->Close();

}