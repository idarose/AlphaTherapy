#include "../include/EnergyDepositionHistograms.hpp"
#include "../include/DecayDynamics.hpp"
#include "../include/CellHit.hpp"

//------------------–----------
EnergyDepositionHistograms::EnergyDepositionHistograms(double EMax_in, int iteration)
{
    EMax = EMax_in;
    it = iteration;
}


//------------------–----------
void EnergyDepositionHistograms::GenerateEmptyHistograms(DecayDynamics decayDynamicsInstance)
{
    int NBins_eVBinning = std::pow(10.,6.);
    double EMin_eVBinning = 0.;
    double EMax_eVBinning = 1.;

    int NBins_keVBinning = 3.*std::pow(10.,4.);
    double EMin_keVBinning = 0.;
    double EMax_keVBinning = 300.;

    int NBins_Hits = 700;
    double MaxHits = 700.;
    double MinHits = 0.;

    int NBins_mGyBinning = 8.*std::pow(10.,4.);
    double MaxDose_mGyBinning = 800.;
    double MinDose_mGyBinning = 0.;

    int NBins_uGyBinning = std::pow(10.,6.);
    double MaxDose_uGyBinning = 10.;
    double MinDose_uGyBinning = 0.;



    //-----------------------------
    // Energy deposition histograms

    //-----------------------------
    // Making names for histograms
    std::string generalHistogramName = "i" + std::to_string(it) + "_hEnergyDeps_212Pb_" + decayDynamicsInstance.GetCellLine() + "_" + std::to_string(decayDynamicsInstance.GetActivity()) + "kBq_";

    std::string histogramNameNucleus_eVBinning = generalHistogramName + "Nucleus_eVBinning";
    std::string histogramNameMembrane_eVBinning = generalHistogramName + "Membrane_eVBinning";
    std::string histogramNameCytoplasm_eVBinning = generalHistogramName + "Cytoplasm_eVBinning";
    std::string histogramNameTotalCell_eVBinning = generalHistogramName + "TotalCell_eVBinning";

    std::string histogramNameNucleus_keVBinning = generalHistogramName + "Nucleus_keVBinning";
    std::string histogramNameMembrane_keVBinning = generalHistogramName + "Membrane_keVBinning";
    std::string histogramNameCytoplasm_keVBinning = generalHistogramName + "Cytoplasm_keVBinning";
    std::string histogramNameTotalCell_keVBinning = generalHistogramName + "TotalCell_keVBinning";

    std::string histogramNameNucleus_FromSolution = generalHistogramName + "Nucleus_FromSolution";
    std::string histogramNameNucleus_FromMembrane = generalHistogramName + "Nucleus_FromMembrane";
    std::string histogramNameNucleus_FromCytoplasm = generalHistogramName + "Nucleus_FromCytoplasm";

    std::string histogramNameMembrane_FromSolution = generalHistogramName + "Membrane_FromSolution";
    std::string histogramNameMembrane_FromMembrane = generalHistogramName + "Membrane_FromMembrane";
    std::string histogramNameMembrane_FromCytoplasm = generalHistogramName + "Membrane_FromCytoplasm";

    std::string histogramNameCytoplasm_FromSolution = generalHistogramName + "Cytoplasm_FromSolution";
    std::string histogramNameCytoplasm_FromMembrane = generalHistogramName + "Cytoplasm_FromMembrane";
    std::string histogramNameCytoplasm_FromCytoplasm = generalHistogramName + "Cytoplasm_FromCytoplasm";

    std::string histogramNameTotalCell_FromSolution = generalHistogramName + "TotalCell_FromSolution";
    std::string histogramNameTotalCell_FromMembrane = generalHistogramName + "TotalCell_FromMembrane";
    std::string histogramNameTotalCell_FromCytoplasm = generalHistogramName + "TotalCell_FromCytoplasm";



    //-------------------------------
    // Making empty histograms

    hEnergyDepsNucleus_eVBinning = new TH1D(histogramNameNucleus_eVBinning.c_str(), "Energy Deposition in Cell Nucleus", NBins_eVBinning, EMin_eVBinning, EMax_eVBinning);
    hEnergyDepsMembrane_eVBinning = new TH1D(histogramNameMembrane_eVBinning.c_str(), "Energy Deposition in Cell Membrane", NBins_eVBinning, EMin_eVBinning, EMax_eVBinning);
    hEnergyDepsCytoplasm_eVBinning = new TH1D(histogramNameCytoplasm_eVBinning.c_str(), "Energy Deposition in Cell Cytoplasm", NBins_eVBinning, EMin_eVBinning, EMax_eVBinning);
    hEnergyDepsTotalCell_eVBinning = new TH1D(histogramNameTotalCell_eVBinning.c_str(), "Energy Deposition in Cell", NBins_eVBinning, EMin_eVBinning, EMax_eVBinning);

    hEnergyDepsNucleus_keVBinning = new TH1D(histogramNameNucleus_keVBinning.c_str(), "Energy Deposition in Cell Nucleus", NBins_keVBinning, EMin_keVBinning, EMax_keVBinning);
    hEnergyDepsMembrane_keVBinning = new TH1D(histogramNameMembrane_keVBinning.c_str(), "Energy Deposition in Cell Membrane ", NBins_keVBinning, EMin_keVBinning, EMax_keVBinning);
    hEnergyDepsCytoplasm_keVBinning = new TH1D(histogramNameCytoplasm_keVBinning.c_str(), "Energy Deposition in Cell Cytoplasm", NBins_keVBinning, EMin_keVBinning, EMax_keVBinning);
    hEnergyDepsTotalCell_keVBinning = new TH1D(histogramNameTotalCell_keVBinning.c_str(), "Energy Deposition in Cell", NBins_keVBinning, EMin_keVBinning, EMax_keVBinning);

    hEnergyDepsNucleus_FromSolution = new TH1D(histogramNameNucleus_FromSolution.c_str(), "Energy Deposition in Cell Nucleus (Decay originated in Solution) ", NBins_keVBinning, EMin_keVBinning, EMax_keVBinning);
    hEnergyDepsNucleus_FromMembrane = new TH1D(histogramNameNucleus_FromMembrane.c_str(), "Energy Deposition in Cell Nucleus (Decay originated in Membrane) ", NBins_keVBinning, EMin_keVBinning, EMax_keVBinning);
    hEnergyDepsNucleus_FromCytoplasm = new TH1D(histogramNameNucleus_FromCytoplasm.c_str(), "Energy Deposition in Cell Nucleus (Decay originated in Cytoplasm) ", NBins_keVBinning, EMin_keVBinning, EMax_keVBinning);

    hEnergyDepsMembrane_FromSolution = new TH1D(histogramNameMembrane_FromSolution.c_str(), "Energy Deposition in Cell Membrane (Decay originated in Solution) ", NBins_keVBinning, EMin_keVBinning, EMax_keVBinning);
    hEnergyDepsMembrane_FromMembrane = new TH1D(histogramNameMembrane_FromMembrane.c_str(), "Energy Deposition in Cell Membrane (Decay originated in Membrane) ", NBins_keVBinning, EMin_keVBinning, EMax_keVBinning);
    hEnergyDepsMembrane_FromCytoplasm = new TH1D(histogramNameMembrane_FromCytoplasm.c_str(), "Energy Deposition in Cell Membrane (Decay originated in Cytoplasm) ", NBins_keVBinning, EMin_keVBinning, EMax_keVBinning);

    hEnergyDepsCytoplasm_FromSolution = new TH1D(histogramNameCytoplasm_FromSolution.c_str(), "Energy Deposition in Cell Cytoplasm (Decay originated in Solution) ", NBins_keVBinning, EMin_keVBinning, EMax_keVBinning);
    hEnergyDepsCytoplasm_FromMembrane = new TH1D(histogramNameCytoplasm_FromMembrane.c_str(), "Energy Deposition in Cell Cytoplasm (Decay originated in Membrane) ", NBins_keVBinning, EMin_keVBinning, EMax_keVBinning);
    hEnergyDepsCytoplasm_FromCytoplasm = new TH1D(histogramNameCytoplasm_FromCytoplasm.c_str(), "Energy Deposition in Cell Cytoplasm (Decay originated in Cytoplasm) ", NBins_keVBinning, EMin_keVBinning, EMax_keVBinning);

    hEnergyDepsTotalCell_FromSolution = new TH1D(histogramNameTotalCell_FromSolution.c_str(), "Energy Deposition in Cell (Decay originated in Solution) ", NBins_keVBinning, EMin_keVBinning, EMax_keVBinning);
    hEnergyDepsTotalCell_FromMembrane = new TH1D(histogramNameTotalCell_FromMembrane.c_str(), "Energy Deposition in Cell (Decay originated in Membrane) ", NBins_keVBinning, EMin_keVBinning, EMax_keVBinning);
    hEnergyDepsTotalCell_FromCytoplasm = new TH1D(histogramNameTotalCell_FromCytoplasm.c_str(), "Energy Deposition in Cell (Decay originated in Cytoplasm) ", NBins_keVBinning, EMin_keVBinning, EMax_keVBinning);


    // naming axis
    hEnergyDepsMembrane_eVBinning->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    hEnergyDepsMembrane_eVBinning->GetYaxis()->SetTitle("Fraction of Cells Hit With this Energy");

    hEnergyDepsCytoplasm_eVBinning->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    hEnergyDepsCytoplasm_eVBinning->GetYaxis()->SetTitle("Fraction of Cells Hit With this Energy");

    hEnergyDepsNucleus_eVBinning->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    hEnergyDepsNucleus_eVBinning->GetYaxis()->SetTitle("Fraction of Cells Hit With this Energy");

    hEnergyDepsTotalCell_eVBinning->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    hEnergyDepsTotalCell_eVBinning->GetYaxis()->SetTitle("Fraction of Cells Hit With this Energy");

    hEnergyDepsMembrane_keVBinning->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    hEnergyDepsMembrane_keVBinning->GetYaxis()->SetTitle("Fraction of Cells Hit With this Energy");

    hEnergyDepsCytoplasm_keVBinning->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    hEnergyDepsCytoplasm_eVBinning->GetYaxis()->SetTitle("Fraction of Cells Hit With this Energy");

    hEnergyDepsNucleus_keVBinning->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    hEnergyDepsNucleus_eVBinning->GetYaxis()->SetTitle("Fraction of Cells Hit With this Energy");

    hEnergyDepsTotalCell_keVBinning->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    hEnergyDepsTotalCell_keVBinning->GetYaxis()->SetTitle("Fraction of Cells Hit With this Energy");

    hEnergyDepsNucleus_FromSolution->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    hEnergyDepsNucleus_FromSolution->GetYaxis()->SetTitle("Fraction of Cells Hit With this Energy");

    hEnergyDepsNucleus_FromMembrane->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    hEnergyDepsNucleus_FromMembrane->GetYaxis()->SetTitle("Fraction of Cells Hit With this Energy");

    hEnergyDepsNucleus_FromCytoplasm->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    hEnergyDepsNucleus_FromCytoplasm->GetYaxis()->SetTitle("Fraction of Cells Hit With this Energy");


    hEnergyDepsMembrane_FromSolution->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    hEnergyDepsMembrane_FromSolution->GetYaxis()->SetTitle("Fraction of Cells Hit With this Energy");

    hEnergyDepsMembrane_FromMembrane->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    hEnergyDepsMembrane_FromMembrane->GetYaxis()->SetTitle("Fraction of Cells Hit With this Energy");

    hEnergyDepsMembrane_FromCytoplasm->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    hEnergyDepsMembrane_FromCytoplasm->GetYaxis()->SetTitle("Fraction of Cells Hit With this Energy");


    hEnergyDepsCytoplasm_FromSolution->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    hEnergyDepsCytoplasm_FromSolution->GetYaxis()->SetTitle("Fraction of Cells Hit With this Energy");

    hEnergyDepsCytoplasm_FromMembrane->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    hEnergyDepsCytoplasm_FromMembrane->GetYaxis()->SetTitle("Fraction of Cells Hit With this Energy");

    hEnergyDepsCytoplasm_FromCytoplasm->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    hEnergyDepsCytoplasm_FromCytoplasm->GetYaxis()->SetTitle("Fraction of Cells Hit With this Energy");


    hEnergyDepsTotalCell_FromSolution->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    hEnergyDepsTotalCell_FromSolution->GetYaxis()->SetTitle("Fraction of Cells Hit With this Energy");

    hEnergyDepsTotalCell_FromMembrane->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    hEnergyDepsTotalCell_FromMembrane->GetYaxis()->SetTitle("Fraction of Cells Hit With this Energy");

    hEnergyDepsTotalCell_FromCytoplasm->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    hEnergyDepsTotalCell_FromCytoplasm->GetYaxis()->SetTitle("Fraction of Cells Hit With this Energy");


    //----------------------------
    // Hit vs energy 2D histograms

    std::string histogramNameTotalCell_HitsAlpha =  generalHistogramName + "HitsAlpha_TotalCell";
    std::string histogramNameNucleus_HitsAlpha = generalHistogramName + "HitsAlpha_Nucleus";
    std::string histogramNameMembrane_HitsAlpha =  generalHistogramName + "HitsAlpha_Membrane";
    std::string histogramNameCytoplasm_HitsAlpha =  generalHistogramName + "HitsAlpha_Cytoplasm";


    // Making empty histograms
    hEnergyDepsTotalCell_HitsAlpha = new TH2D(histogramNameTotalCell_HitsAlpha.c_str(), "Energy Depsition in Cells per N Number of Hits to Cell by Alphas", NBins_keVBinning, EMin_keVBinning, EMax_keVBinning, NBins_Hits, MinHits, MaxHits);
    hEnergyDepsNucleus_HitsAlpha = new TH2D(histogramNameNucleus_HitsAlpha.c_str(), "Energy Depsition in Cell Nuclei per N Number of Hits to Cell Nucleus by Alphas", NBins_keVBinning, EMin_keVBinning, EMax_keVBinning, NBins_Hits, MinHits, MaxHits);
    hEnergyDepsMembrane_HitsAlpha = new TH2D(histogramNameMembrane_HitsAlpha.c_str(), "Energy Depsition in Cell Membranes per N Number of Hits to Cell Membrane by Alphas", NBins_keVBinning, EMin_keVBinning, EMax_keVBinning, NBins_Hits, MinHits, MaxHits);
    hEnergyDepsCytoplasm_HitsAlpha = new TH2D(histogramNameCytoplasm_HitsAlpha.c_str(), "Energy Depsition in Cell Cytoplasms per N Number of Hits to Cell Cytoplasm by Alphas", NBins_keVBinning, EMin_keVBinning, EMax_keVBinning, NBins_Hits, MinHits, MaxHits);


    // Naming axis
    hEnergyDepsTotalCell_HitsAlpha->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    hEnergyDepsTotalCell_HitsAlpha->GetYaxis()->SetTitle("N Number of Hits by Alpha Particle");

    hEnergyDepsNucleus_HitsAlpha->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    hEnergyDepsNucleus_HitsAlpha->GetYaxis()->SetTitle("N Number of Hits by Alpha Particle");

    hEnergyDepsMembrane_HitsAlpha->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    hEnergyDepsMembrane_HitsAlpha->GetYaxis()->SetTitle("N Number of Hits by Alpha Particle");

    hEnergyDepsCytoplasm_HitsAlpha->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    hEnergyDepsCytoplasm_HitsAlpha->GetYaxis()->SetTitle("N Number of Hits by Alpha Particle");



    //--------------------------
    // Hit multiplicity histograms

    std::string histogramNameFractionHitsAlpha_TotalCell = "i" + std::to_string(it) + "_hFractionHitsAlpha_212Pb_" + decayDynamicsInstance.GetCellLine() + "_" + std::to_string(decayDynamicsInstance.GetActivity()) + "kBq_TotalCell";
    std::string histogramNameFractionHitsAlpha_Nucleus = "i" + std::to_string(it) + "_hFractionHitsAlpha_212Pb_" + decayDynamicsInstance.GetCellLine() + "_" + std::to_string(decayDynamicsInstance.GetActivity()) + "kBq_Nucleus";
    std::string histogramNameFractionHitsAlpha_Membrane = "i" + std::to_string(it) + "_hFractionHitsAlpha_212Pb_" + decayDynamicsInstance.GetCellLine() + "_" + std::to_string(decayDynamicsInstance.GetActivity()) + "kBq_Membrane";
    std::string histogramNameFractionHitsAlpha_Cytoplasm = "i" + std::to_string(it) + "_hFractionHitsAlpha_212Pb_" + decayDynamicsInstance.GetCellLine() + "_" + std::to_string(decayDynamicsInstance.GetActivity()) + "kBq_Cytoplasm";

    // Making empty histograms
    hFractionHitsAlpha_TotalCell = new TH1D(histogramNameFractionHitsAlpha_TotalCell.c_str(), "Fraction of Cells Hit N Number of Times by an Alpha Particle", NBins_Hits, MinHits, MaxHits);
    hFractionHitsAlpha_Nucleus = new TH1D(histogramNameFractionHitsAlpha_Nucleus.c_str(), "Fraction of Cell Nuclei Hit N Number Number of Times by an Alpha Particle", NBins_Hits, MinHits, MaxHits);
    hFractionHitsAlpha_Membrane = new TH1D(histogramNameFractionHitsAlpha_Membrane.c_str(), "Fraction of Cell Membranes Hit N Number of Times by an Alpha Particle", NBins_Hits, MinHits, MaxHits);
    hFractionHitsAlpha_Cytoplasm = new TH1D(histogramNameFractionHitsAlpha_Cytoplasm.c_str(), "Fraction of Cell Cytoplasms Hit N Number of Times by an Alpha Particle", NBins_Hits, MinHits, MaxHits);

    // Naming axis
    hFractionHitsAlpha_TotalCell->GetXaxis()->SetTitle("N Number of hits by Alpha Particle");
    hFractionHitsAlpha_TotalCell->GetYaxis()->SetTitle("Fraction of cells hit");

    hFractionHitsAlpha_Nucleus->GetXaxis()->SetTitle("N Number of hits by Alpha Particle");
    hFractionHitsAlpha_Nucleus->GetYaxis()->SetTitle("Fraction of cells hit");

    hFractionHitsAlpha_Membrane->GetXaxis()->SetTitle("N Number of hits by Alpha Particle");
    hFractionHitsAlpha_Membrane->GetYaxis()->SetTitle("Fraction of cells hit");

    hFractionHitsAlpha_Cytoplasm->GetXaxis()->SetTitle("N Number of hits by Alpha Particle");
    hFractionHitsAlpha_Cytoplasm->GetYaxis()->SetTitle("Fraction of cells hit");


    //-----------------------------
    // Dose histograms

    generalHistogramName = "i" + std::to_string(it) + "_hDose_212Pb_" + decayDynamicsInstance.GetCellLine() + "_" + std::to_string(decayDynamicsInstance.GetActivity()) + "kBq_";

    std::string histogramNameNucleus_uGyBinning = generalHistogramName + "Nucleus_uGyBinning";
    std::string histogramNameMembrane_uGyBinning = generalHistogramName + "Membrane_uGyBinning";
    std::string histogramNameCytoplasm_uGyBinning = generalHistogramName + "Cytoplasm_uGyBinning";
    std::string histogramNameTotalCell_uGyBinning = generalHistogramName + "TotalCell_uGyBinning";

    std::string histogramNameNucleus_mGyBinning = generalHistogramName + "Nucleus_mGyBinning";
    std::string histogramNameMembrane_mGyBinning = generalHistogramName + "Membrane_mGyBinning";
    std::string histogramNameCytoplasm_mGyBinning = generalHistogramName + "Cytoplasm_mGyBinning";
    std::string histogramNameTotalCell_mGyBinning = generalHistogramName + "TotalCell_mGyBinning";

    histogramNameNucleus_FromSolution = generalHistogramName + "Nucleus_FromSolution";
    histogramNameNucleus_FromMembrane = generalHistogramName + "Nucleus_FromMembrane";
    histogramNameNucleus_FromCytoplasm = generalHistogramName + "Nucleus_FromCytoplasm";

    histogramNameMembrane_FromSolution = generalHistogramName + "Membrane_FromSolution";
    histogramNameMembrane_FromMembrane = generalHistogramName + "Membrane_FromMembrane";
    histogramNameMembrane_FromCytoplasm = generalHistogramName + "Membrane_FromCytoplasm";

    histogramNameCytoplasm_FromSolution = generalHistogramName + "Cytoplasm_FromSolution";
    histogramNameCytoplasm_FromMembrane = generalHistogramName + "Cytoplasm_FromMembrane";
    histogramNameCytoplasm_FromCytoplasm = generalHistogramName + "Cytoplasm_FromCytoplasm";

    histogramNameTotalCell_FromSolution = generalHistogramName + "TotalCell_FromSolution";
    histogramNameTotalCell_FromMembrane = generalHistogramName + "TotalCell_FromMembrane";
    histogramNameTotalCell_FromCytoplasm = generalHistogramName + "TotalCell_FromCytoplasm";

    //-------------------------------
    // Making empty histograms
    hDoseNucleus_uGyBinning = new TH1D(histogramNameNucleus_uGyBinning.c_str(), "Dose Delived in Cell Nucleus", NBins_uGyBinning, MinDose_uGyBinning, MaxDose_uGyBinning);
    hDoseMembrane_uGyBinning = new TH1D(histogramNameMembrane_uGyBinning.c_str(), "Dose Delivered in Cell Membrane", NBins_uGyBinning, MinDose_uGyBinning, MaxDose_uGyBinning);
    hDoseCytoplasm_uGyBinning = new TH1D(histogramNameCytoplasm_uGyBinning.c_str(), "Dose Delivered in Cell Cytoplasm", NBins_uGyBinning, MinDose_uGyBinning, MaxDose_uGyBinning);
    hDoseTotalCell_uGyBinning = new TH1D(histogramNameTotalCell_uGyBinning.c_str(), "Dose Delivered in Cell", NBins_uGyBinning, MinDose_uGyBinning, MaxDose_uGyBinning);

    hDoseNucleus_mGyBinning = new TH1D(histogramNameNucleus_mGyBinning.c_str(), "Dose Delivered in Cell Nucleus", NBins_mGyBinning, MinDose_mGyBinning, MaxDose_mGyBinning);
    hDoseMembrane_mGyBinning = new TH1D(histogramNameMembrane_mGyBinning.c_str(), "Dose Delivered in Cell Membrane", NBins_mGyBinning, MinDose_mGyBinning, MaxDose_mGyBinning);
    hDoseCytoplasm_mGyBinning = new TH1D(histogramNameCytoplasm_mGyBinning.c_str(), "Dose Delivered in Cell Cytoplasm", NBins_mGyBinning, MinDose_mGyBinning, MaxDose_mGyBinning);
    hDoseTotalCell_mGyBinning = new TH1D(histogramNameTotalCell_mGyBinning.c_str(), "Dose Delivered in Cell", NBins_mGyBinning, MinDose_mGyBinning, MaxDose_mGyBinning);

    hDoseNucleus_FromSolution = new TH1D(histogramNameNucleus_FromSolution.c_str(), "Dose Delivered in Cell Nucleus (Decay originated in Solution) ", NBins_mGyBinning, MinDose_mGyBinning, MaxDose_mGyBinning);
    hDoseNucleus_FromMembrane = new TH1D(histogramNameNucleus_FromMembrane.c_str(), "Dose Delivered in Cell Nucleus (Decay originated in Membrane) ", NBins_mGyBinning, MinDose_mGyBinning, MaxDose_mGyBinning);
    hDoseNucleus_FromCytoplasm = new TH1D(histogramNameNucleus_FromCytoplasm.c_str(), "Dose Delivered in Cell Nucleus (Decay originated in Cytoplasm) ", NBins_mGyBinning, MinDose_mGyBinning, MaxDose_mGyBinning);

    hDoseMembrane_FromSolution = new TH1D(histogramNameMembrane_FromSolution.c_str(), "Dose Delivered in Cell Membrane (Decay originated in Solution) ", NBins_mGyBinning, MinDose_mGyBinning, MaxDose_mGyBinning);
    hDoseMembrane_FromMembrane = new TH1D(histogramNameMembrane_FromMembrane.c_str(), "Dose Delivered in Cell Membrane (Decay originated in Membrane) ", NBins_mGyBinning, MinDose_mGyBinning, MaxDose_mGyBinning);
    hDoseMembrane_FromCytoplasm = new TH1D(histogramNameMembrane_FromCytoplasm.c_str(), "Dose Delivered in Cell Membrane (Decay originated in Cytoplasm) ", NBins_mGyBinning, MinDose_mGyBinning, MaxDose_mGyBinning);

    hDoseCytoplasm_FromSolution = new TH1D(histogramNameCytoplasm_FromSolution.c_str(), "Dose Delivered in Cell Cytoplasm (Decay originated in Solution) ", NBins_mGyBinning, MinDose_mGyBinning, MaxDose_mGyBinning);
    hDoseCytoplasm_FromMembrane = new TH1D(histogramNameCytoplasm_FromMembrane.c_str(), "Dose Delivered in Cell Cytoplasm (Decay originated in Membrane) ", NBins_mGyBinning, MinDose_mGyBinning, MaxDose_mGyBinning);
    hDoseCytoplasm_FromCytoplasm = new TH1D(histogramNameCytoplasm_FromCytoplasm.c_str(), "Dose Delivered in Cell Cytoplasm (Decay originated in Cytoplasm) ", NBins_mGyBinning, MinDose_mGyBinning, MaxDose_mGyBinning);

    hDoseTotalCell_FromSolution = new TH1D(histogramNameTotalCell_FromSolution.c_str(), "Dose Delivered in Cell (Decay originated in Solution) ", NBins_mGyBinning, MinDose_mGyBinning, MaxDose_mGyBinning);
    hDoseTotalCell_FromMembrane = new TH1D(histogramNameTotalCell_FromMembrane.c_str(), "Dose Delivered in Cell (Decay originated in Membrane) ", NBins_mGyBinning, MinDose_mGyBinning, MaxDose_mGyBinning);
    hDoseTotalCell_FromCytoplasm = new TH1D(histogramNameTotalCell_FromCytoplasm.c_str(), "Dose Delivered in Cell (Decay originated in Cytoplasm) ", NBins_mGyBinning, MinDose_mGyBinning, MaxDose_mGyBinning);

    // Naming axis
    hDoseMembrane_uGyBinning->GetXaxis()->SetTitle("Dose Delivered [Gy]");
    hDoseMembrane_uGyBinning->GetYaxis()->SetTitle("Fraction of Cells Hit with this Dose");

    hDoseCytoplasm_uGyBinning->GetXaxis()->SetTitle("Dose Delivered [Gy]");
    hDoseCytoplasm_uGyBinning->GetYaxis()->SetTitle("Fraction of Cells Hit with This Dose");

    hDoseNucleus_uGyBinning->GetXaxis()->SetTitle("Dose Delivered [Gy]");
    hDoseNucleus_uGyBinning->GetYaxis()->SetTitle("Fraction of Cells Hit with This Dose");

    hDoseTotalCell_uGyBinning->GetXaxis()->SetTitle("Dose Delivered [Gy]");
    hDoseTotalCell_uGyBinning->GetYaxis()->SetTitle("Fraction of Cells Hit with This Dose");

    hDoseMembrane_mGyBinning->GetXaxis()->SetTitle("Dose Delivered [Gy]");
    hDoseMembrane_mGyBinning->GetYaxis()->SetTitle("Fraction of Cells Hit with This Dose");

    hDoseCytoplasm_mGyBinning->GetXaxis()->SetTitle("Dose Delivered [Gy]");
    hDoseCytoplasm_uGyBinning->GetYaxis()->SetTitle("Fraction of Cells Hit with This Dose");

    hDoseNucleus_mGyBinning->GetXaxis()->SetTitle("Dose Delivered [Gy]");
    hDoseNucleus_uGyBinning->GetYaxis()->SetTitle("Fraction of Cells Hit with This Dose");

    hDoseTotalCell_mGyBinning->GetXaxis()->SetTitle("Dose Delivered [Gy]");
    hDoseTotalCell_mGyBinning->GetYaxis()->SetTitle("Fraction of Cells Hit with This Dose");

    hDoseNucleus_FromSolution->GetXaxis()->SetTitle("Dose Delivered [Gy]");
    hDoseNucleus_FromSolution->GetYaxis()->SetTitle("Fraction of Cells Hit with This Dose");

    hDoseNucleus_FromMembrane->GetXaxis()->SetTitle("Dose Delivered [Gy]");
    hDoseNucleus_FromMembrane->GetYaxis()->SetTitle("Fraction of Cells Hit with This Dose");

    hDoseNucleus_FromCytoplasm->GetXaxis()->SetTitle("Dose Delivered [Gy]");
    hDoseNucleus_FromCytoplasm->GetYaxis()->SetTitle("Fraction of Cells Hit with This Dose");


    hDoseMembrane_FromSolution->GetXaxis()->SetTitle("Dose Delivered [Gy]");
    hDoseMembrane_FromSolution->GetYaxis()->SetTitle("Fraction of Cells Hit with This Dose");

    hDoseMembrane_FromMembrane->GetXaxis()->SetTitle("Dose Delivered [Gy]");
    hDoseMembrane_FromMembrane->GetYaxis()->SetTitle("Fraction of Cells Hit with This Dose");

    hDoseMembrane_FromCytoplasm->GetXaxis()->SetTitle("Dose Delivered [Gy]");
    hDoseMembrane_FromCytoplasm->GetYaxis()->SetTitle("Fraction of Cells Hit with This Dose");


    hDoseCytoplasm_FromSolution->GetXaxis()->SetTitle("Dose Delivered [Gy]");
    hDoseCytoplasm_FromSolution->GetYaxis()->SetTitle("Fraction of Cells Hit with This Dose");

    hDoseCytoplasm_FromMembrane->GetXaxis()->SetTitle("Dose Delivered [Gy]");
    hDoseCytoplasm_FromMembrane->GetYaxis()->SetTitle("Fraction of Cells Hit with This Dose");

    hDoseCytoplasm_FromCytoplasm->GetXaxis()->SetTitle("Dose Delivered [Gy]");
    hDoseCytoplasm_FromCytoplasm->GetYaxis()->SetTitle("Fraction of Cells Hit with This Dose");


    hDoseTotalCell_FromSolution->GetXaxis()->SetTitle("Dose Delivered [Gy]");
    hDoseTotalCell_FromSolution->GetYaxis()->SetTitle("Fraction of Cells Hit with This Dose");

    hDoseTotalCell_FromMembrane->GetXaxis()->SetTitle("Dose Delivered [Gy]");
    hDoseTotalCell_FromMembrane->GetYaxis()->SetTitle("Fraction of Cells Hit with This Dose");

    hDoseTotalCell_FromCytoplasm->GetXaxis()->SetTitle("Dose Delivered [Gy]");
    hDoseTotalCell_FromCytoplasm->GetYaxis()->SetTitle("Fraction of Cells Hit with This Dose");

    //----------------------------
    // Hit vs dose 2D histograms

    histogramNameTotalCell_HitsAlpha =  generalHistogramName + "HitsAlpha_TotalCell";
    histogramNameNucleus_HitsAlpha = generalHistogramName + "HitsAlpha_Nucleus";
    histogramNameMembrane_HitsAlpha =  generalHistogramName + "HitsAlpha_Membrane";
    histogramNameCytoplasm_HitsAlpha =  generalHistogramName + "HitsAlpha_Cytoplasm";


    // making empty histograms
    hDoseTotalCell_HitsAlpha = new TH2D(histogramNameTotalCell_HitsAlpha.c_str(), "Dose Delivered in Cell per N Number of hits by alphas", NBins_mGyBinning, MinDose_mGyBinning, MaxDose_mGyBinning, NBins_Hits, MinHits, MaxHits);
    hDoseNucleus_HitsAlpha = new TH2D(histogramNameNucleus_HitsAlpha.c_str(), "Dose Delivered in Cell Nucleus per N Number of hits by alphas", NBins_mGyBinning, MinDose_mGyBinning, MaxDose_mGyBinning, NBins_Hits, MinHits, MaxHits);
    hDoseMembrane_HitsAlpha = new TH2D(histogramNameMembrane_HitsAlpha.c_str(), "Dose Delivered in Cell Membrane per N Number of hits by alphas", NBins_mGyBinning, MinDose_mGyBinning, MaxDose_mGyBinning, NBins_Hits, MinHits, MaxHits);
    hDoseCytoplasm_HitsAlpha = new TH2D(histogramNameCytoplasm_HitsAlpha.c_str(), "Dose Delivered in Cell Cytoplasm per N Number of hits by alphas", NBins_mGyBinning, MinDose_mGyBinning, MaxDose_mGyBinning, NBins_Hits, MinHits, MaxHits);


    // Naming axis
    hDoseTotalCell_HitsAlpha->GetXaxis()->SetTitle("Dose Delivered[Gy]");
    hDoseTotalCell_HitsAlpha->GetYaxis()->SetTitle("N Number of Hits by Alpha Particle");

    hDoseNucleus_HitsAlpha->GetXaxis()->SetTitle("Dose Delivered[Gy]");
    hDoseNucleus_HitsAlpha->GetYaxis()->SetTitle("N Number of Hits by Alpha Particle");

    hDoseMembrane_HitsAlpha->GetXaxis()->SetTitle("Dose Delivered[Gy]");
    hDoseMembrane_HitsAlpha->GetYaxis()->SetTitle("N Number of Hits by Alpha Particle");

    hDoseCytoplasm_HitsAlpha->GetXaxis()->SetTitle("Dose Delivered[Gy]");
    hDoseCytoplasm_HitsAlpha->GetYaxis()->SetTitle("N Number of Hits by Alpha Particle");
}

//------------------–----------
void EnergyDepositionHistograms::AddCellHitsToHistograms(CellHit cellHit)
{
    double massMembrane = cellHit.GetMassMembrane();
    double massNucleus = cellHit.GetMassNucleus();
    double massCytoplasm = cellHit.GetMassCytoplasm();
    double massCell = cellHit.GetMassCell();

    double convert_MeV_joules = 1.60218e-13;

    //----------------------------
    // Fill histograms

    //--------------------------------
    // Filling for energy depositions
    if(cellHit.GetEnergyDepositionMembrane()>0.0)
    {
        hEnergyDepsMembrane_eVBinning->Fill(cellHit.GetEnergyDepositionMembrane());
        hEnergyDepsMembrane_keVBinning->Fill(cellHit.GetEnergyDepositionMembrane());

        hDoseMembrane_uGyBinning->Fill(cellHit.GetEnergyDepositionMembrane()*convert_MeV_joules/massMembrane);
        hDoseMembrane_mGyBinning->Fill(cellHit.GetEnergyDepositionMembrane()*convert_MeV_joules/massMembrane);
    }
    if(cellHit.GetEnergyDepositionCytoplasm()>0.0)
    {
        hEnergyDepsCytoplasm_eVBinning->Fill(cellHit.GetEnergyDepositionCytoplasm());
        hEnergyDepsCytoplasm_keVBinning->Fill(cellHit.GetEnergyDepositionCytoplasm());

        hDoseCytoplasm_uGyBinning->Fill(cellHit.GetEnergyDepositionCytoplasm()*convert_MeV_joules/massCytoplasm);
        hDoseCytoplasm_mGyBinning->Fill(cellHit.GetEnergyDepositionCytoplasm()*convert_MeV_joules/massCytoplasm);
    }
    if(cellHit.GetEnergyDepositionNucleus()>0.0)
    {
        hEnergyDepsNucleus_eVBinning->Fill(cellHit.GetEnergyDepositionNucleus());
        hEnergyDepsNucleus_keVBinning->Fill(cellHit.GetEnergyDepositionNucleus());

        hDoseNucleus_uGyBinning->Fill(cellHit.GetEnergyDepositionNucleus()*convert_MeV_joules/massNucleus);
        hDoseNucleus_mGyBinning->Fill(cellHit.GetEnergyDepositionNucleus()*convert_MeV_joules/massNucleus);
    }
    if(cellHit.GetEnergyDepositionTotalCell()>0.0)
    {
        hEnergyDepsTotalCell_eVBinning->Fill(cellHit.GetEnergyDepositionTotalCell());
        hEnergyDepsTotalCell_keVBinning->Fill(cellHit.GetEnergyDepositionTotalCell());

        hDoseTotalCell_uGyBinning->Fill(cellHit.GetEnergyDepositionTotalCell()*convert_MeV_joules/massCell);
        hDoseTotalCell_mGyBinning->Fill(cellHit.GetEnergyDepositionTotalCell()*convert_MeV_joules/massCell);
    }

    //------------------------------------
    // Filling for decay originating in different components
    if(cellHit.GetEnergyDepositionNucleus_FromSolution()>0.0)
    {
        hEnergyDepsNucleus_FromSolution->Fill(cellHit.GetEnergyDepositionNucleus_FromSolution());

        hDoseNucleus_FromSolution->Fill(cellHit.GetEnergyDepositionNucleus_FromSolution()*convert_MeV_joules/massNucleus);
    }
    if(cellHit.GetEnergyDepositionNucleus_FromMembrane()>0.0)
    {
        hEnergyDepsNucleus_FromMembrane->Fill(cellHit.GetEnergyDepositionNucleus_FromMembrane());

        hDoseNucleus_FromMembrane->Fill(cellHit.GetEnergyDepositionNucleus_FromMembrane()*convert_MeV_joules/massNucleus);
    }
    if(cellHit.GetEnergyDepositionNucleus_FromCytoplasm()>0.0)
    {
        hEnergyDepsNucleus_FromCytoplasm->Fill(cellHit.GetEnergyDepositionNucleus_FromCytoplasm());

        hDoseNucleus_FromCytoplasm->Fill(cellHit.GetEnergyDepositionNucleus_FromCytoplasm()*convert_MeV_joules/massNucleus);
    }
    if(cellHit.GetEnergyDepositionMembrane_FromSolution()>0.0)
    {
        hEnergyDepsMembrane_FromSolution->Fill(cellHit.GetEnergyDepositionMembrane_FromSolution());

        hDoseMembrane_FromSolution->Fill(cellHit.GetEnergyDepositionMembrane_FromSolution()*convert_MeV_joules/massMembrane);
    }
    if(cellHit.GetEnergyDepositionMembrane_FromMembrane()>0.0)
    {
        hEnergyDepsMembrane_FromMembrane->Fill(cellHit.GetEnergyDepositionMembrane_FromMembrane());

        hDoseMembrane_FromMembrane->Fill(cellHit.GetEnergyDepositionMembrane_FromMembrane()*convert_MeV_joules/massMembrane);
    }
    if(cellHit.GetEnergyDepositionMembrane_FromCytoplasm()>0.0)
    {
        hEnergyDepsMembrane_FromCytoplasm->Fill(cellHit.GetEnergyDepositionMembrane_FromCytoplasm());

        hDoseMembrane_FromCytoplasm->Fill(cellHit.GetEnergyDepositionMembrane_FromCytoplasm()*convert_MeV_joules/massMembrane);
    }
    if(cellHit.GetEnergyDepositionCytoplasm_FromSolution()>0.0)
    {
        hEnergyDepsCytoplasm_FromSolution->Fill(cellHit.GetEnergyDepositionCytoplasm_FromSolution());

        hDoseCytoplasm_FromSolution->Fill(cellHit.GetEnergyDepositionCytoplasm_FromSolution()*convert_MeV_joules/massCytoplasm);
    }
    if(cellHit.GetEnergyDepositionCytoplasm_FromMembrane()>0.0)
    {
        hEnergyDepsCytoplasm_FromMembrane->Fill(cellHit.GetEnergyDepositionCytoplasm_FromMembrane());

        hDoseCytoplasm_FromMembrane->Fill(cellHit.GetEnergyDepositionCytoplasm_FromMembrane()*convert_MeV_joules/massCytoplasm);
    }
    if(cellHit.GetEnergyDepositionCytoplasm_FromCytoplasm()>0.0)
    {
        hEnergyDepsCytoplasm_FromCytoplasm->Fill(cellHit.GetEnergyDepositionCytoplasm_FromCytoplasm());

        hDoseCytoplasm_FromCytoplasm->Fill(cellHit.GetEnergyDepositionCytoplasm_FromCytoplasm()*convert_MeV_joules/massCytoplasm);
    }
    if(cellHit.GetEnergyDepositionTotalCell_FromSolution()>0.0)
    {
        hEnergyDepsTotalCell_FromSolution->Fill(cellHit.GetEnergyDepositionTotalCell_FromSolution());

        hDoseTotalCell_FromSolution->Fill(cellHit.GetEnergyDepositionTotalCell_FromSolution()*convert_MeV_joules/massCell);
    }
    if(cellHit.GetEnergyDepositionTotalCell_FromMembrane()>0.0)
    {
        hEnergyDepsTotalCell_FromMembrane->Fill(cellHit.GetEnergyDepositionTotalCell_FromMembrane());

        hDoseTotalCell_FromMembrane->Fill(cellHit.GetEnergyDepositionTotalCell_FromMembrane()*convert_MeV_joules/massCell);
    }
    if(cellHit.GetEnergyDepositionTotalCell_FromCytoplasm()>0.0)
    {
        hEnergyDepsTotalCell_FromCytoplasm->Fill(cellHit.GetEnergyDepositionTotalCell_FromCytoplasm());

        hDoseTotalCell_FromCytoplasm->Fill(cellHit.GetEnergyDepositionTotalCell_FromCytoplasm()*convert_MeV_joules/massCell);
    }


    //------------------------------
    // Filling 2D histograms for energy deposition and number of hits
    if(cellHit.GetEnergyDepositionTotalCell()>0.0)
    {
        hEnergyDepsTotalCell_HitsAlpha->Fill(cellHit.GetEnergyDepositionTotalCell(), cellHit.GetNumberHitsAlphas_TotalCell());

        hDoseTotalCell_HitsAlpha->Fill(cellHit.GetEnergyDepositionTotalCell()*convert_MeV_joules/massCell, cellHit.GetNumberHitsAlphas_TotalCell());
    }
    if(cellHit.GetEnergyDepositionNucleus()>0.0)
    {
        hEnergyDepsNucleus_HitsAlpha->Fill(cellHit.GetEnergyDepositionNucleus(), cellHit.GetNumberHitsAlphas_Nucleus());

        hDoseNucleus_HitsAlpha->Fill(cellHit.GetEnergyDepositionNucleus()*convert_MeV_joules/massNucleus, cellHit.GetNumberHitsAlphas_Nucleus());
    }
    if(cellHit.GetEnergyDepositionMembrane()>0.0)
    {
        hEnergyDepsMembrane_HitsAlpha->Fill(cellHit.GetEnergyDepositionMembrane(), cellHit.GetNumberHitsAlphas_Membrane());

        hDoseMembrane_HitsAlpha->Fill(cellHit.GetEnergyDepositionMembrane()*convert_MeV_joules/massMembrane, cellHit.GetNumberHitsAlphas_Membrane());
    }
    if(cellHit.GetEnergyDepositionCytoplasm()>0.0)
    {
        hEnergyDepsCytoplasm_HitsAlpha->Fill(cellHit.GetEnergyDepositionCytoplasm(), cellHit.GetNumberHitsAlphas_Cytoplasm());

        hDoseCytoplasm_HitsAlpha->Fill(cellHit.GetEnergyDepositionCytoplasm()*convert_MeV_joules/massCytoplasm, cellHit.GetNumberHitsAlphas_Cytoplasm());
    }


    //-------------------------------
    // Filling hit multiplicity histograms
    if(cellHit.GetEnergyDepositionTotalCell()>0.0){hFractionHitsAlpha_TotalCell->Fill(cellHit.GetNumberHitsAlphas_TotalCell());}
    if(cellHit.GetEnergyDepositionNucleus()>0.0){hFractionHitsAlpha_Nucleus->Fill(cellHit.GetNumberHitsAlphas_Nucleus());}
    if(cellHit.GetEnergyDepositionMembrane()>0.0){hFractionHitsAlpha_Membrane->Fill(cellHit.GetNumberHitsAlphas_Membrane());}
    if(cellHit.GetEnergyDepositionCytoplasm()>0.0){hFractionHitsAlpha_Cytoplasm->Fill(cellHit.GetNumberHitsAlphas_Cytoplasm());}

}

//------------------–----------
void EnergyDepositionHistograms::ScaleHistograms(double factor)
{
    //---------------------------
    // Scaling histograms

    //--------------------
    // Energy deposition histograms
    hEnergyDepsMembrane_eVBinning->Scale(factor);
    hEnergyDepsCytoplasm_eVBinning->Scale(factor);
    hEnergyDepsNucleus_eVBinning->Scale(factor);
    hEnergyDepsTotalCell_eVBinning->Scale(factor);

    hEnergyDepsMembrane_keVBinning->Scale(factor);
    hEnergyDepsCytoplasm_keVBinning->Scale(factor);
    hEnergyDepsNucleus_keVBinning->Scale(factor);
    hEnergyDepsTotalCell_keVBinning->Scale(factor);

    hEnergyDepsNucleus_FromSolution->Scale(factor);
    hEnergyDepsNucleus_FromMembrane->Scale(factor);
    hEnergyDepsNucleus_FromCytoplasm->Scale(factor);

    hEnergyDepsMembrane_FromSolution->Scale(factor);
    hEnergyDepsMembrane_FromMembrane->Scale(factor);
    hEnergyDepsMembrane_FromCytoplasm->Scale(factor);

    hEnergyDepsCytoplasm_FromSolution->Scale(factor);
    hEnergyDepsCytoplasm_FromMembrane->Scale(factor);
    hEnergyDepsCytoplasm_FromCytoplasm->Scale(factor);

    hEnergyDepsTotalCell_FromSolution->Scale(factor);
    hEnergyDepsTotalCell_FromMembrane->Scale(factor);
    hEnergyDepsTotalCell_FromCytoplasm->Scale(factor);

    //--------------------------------
    // Dose histograms

    hDoseMembrane_uGyBinning->Scale(factor);
    hDoseCytoplasm_uGyBinning->Scale(factor);
    hDoseNucleus_uGyBinning->Scale(factor);
    hDoseTotalCell_uGyBinning->Scale(factor);

    hDoseMembrane_mGyBinning->Scale(factor);
    hDoseCytoplasm_mGyBinning->Scale(factor);
    hDoseNucleus_mGyBinning->Scale(factor);
    hDoseTotalCell_mGyBinning->Scale(factor);

    hDoseNucleus_FromSolution->Scale(factor);
    hDoseNucleus_FromMembrane->Scale(factor);
    hDoseNucleus_FromCytoplasm->Scale(factor);

    hDoseMembrane_FromSolution->Scale(factor);
    hDoseMembrane_FromMembrane->Scale(factor);
    hDoseMembrane_FromCytoplasm->Scale(factor);

    hDoseCytoplasm_FromSolution->Scale(factor);
    hDoseCytoplasm_FromMembrane->Scale(factor);
    hDoseCytoplasm_FromCytoplasm->Scale(factor);

    hDoseTotalCell_FromSolution->Scale(factor);
    hDoseTotalCell_FromMembrane->Scale(factor);
    hDoseTotalCell_FromCytoplasm->Scale(factor);

    //--------------------------------
    // Number of times hit vs energy deposition 2D histograms
    hEnergyDepsTotalCell_HitsAlpha->Scale(factor);
    hEnergyDepsNucleus_HitsAlpha->Scale(factor);
    hEnergyDepsMembrane_HitsAlpha->Scale(factor);
    hEnergyDepsCytoplasm_HitsAlpha->Scale(factor);

    hDoseTotalCell_HitsAlpha->Scale(factor);
    hDoseNucleus_HitsAlpha->Scale(factor);
    hDoseMembrane_HitsAlpha->Scale(factor);
    hDoseCytoplasm_HitsAlpha->Scale(factor);

    //-------------------------------
    // Hit multiplicity histograms
    hFractionHitsAlpha_TotalCell->Scale(factor);
    hFractionHitsAlpha_Nucleus->Scale(factor);
    hFractionHitsAlpha_Membrane->Scale(factor);
    hFractionHitsAlpha_Cytoplasm->Scale(factor);
}

//------------------–----------
void EnergyDepositionHistograms::WriteHistogramsToFile()
{
    //------------------–----------
    // Writing histograms to file

    //------------------------
    // Energy deposition histograms
    hEnergyDepsMembrane_eVBinning->Write();
    hEnergyDepsCytoplasm_eVBinning->Write();
    hEnergyDepsNucleus_eVBinning->Write();
    hEnergyDepsTotalCell_eVBinning->Write();

    hEnergyDepsMembrane_keVBinning->Write();
    hEnergyDepsCytoplasm_keVBinning->Write();
    hEnergyDepsNucleus_keVBinning->Write();
    hEnergyDepsTotalCell_keVBinning->Write();

    hEnergyDepsNucleus_FromSolution->Write();
    hEnergyDepsNucleus_FromMembrane->Write();
    hEnergyDepsNucleus_FromCytoplasm->Write();

    hEnergyDepsMembrane_FromSolution->Write();
    hEnergyDepsMembrane_FromMembrane->Write();
    hEnergyDepsMembrane_FromCytoplasm->Write();

    hEnergyDepsCytoplasm_FromSolution->Write();
    hEnergyDepsCytoplasm_FromMembrane->Write();
    hEnergyDepsCytoplasm_FromCytoplasm->Write();

    hEnergyDepsTotalCell_FromSolution->Write();
    hEnergyDepsTotalCell_FromMembrane->Write();
    hEnergyDepsTotalCell_FromCytoplasm->Write();

    //-----------------------------
    // Dose Histograms

    hDoseMembrane_uGyBinning->Write();
    hDoseCytoplasm_uGyBinning->Write();
    hDoseNucleus_uGyBinning->Write();
    hDoseTotalCell_uGyBinning->Write();

    hDoseMembrane_mGyBinning->Write();
    hDoseCytoplasm_mGyBinning->Write();
    hDoseNucleus_mGyBinning->Write();
    hDoseTotalCell_mGyBinning->Write();

    hDoseNucleus_FromSolution->Write();
    hDoseNucleus_FromMembrane->Write();
    hDoseNucleus_FromCytoplasm->Write();

    hDoseMembrane_FromSolution->Write();
    hDoseMembrane_FromMembrane->Write();
    hDoseMembrane_FromCytoplasm->Write();

    hDoseCytoplasm_FromSolution->Write();
    hDoseCytoplasm_FromMembrane->Write();
    hDoseCytoplasm_FromCytoplasm->Write();

    hDoseTotalCell_FromSolution->Write();
    hDoseTotalCell_FromMembrane->Write();
    hDoseTotalCell_FromCytoplasm->Write();



    //--------------------------------
    // Number of times hit vs energy deposition 2D histograms
    hEnergyDepsTotalCell_HitsAlpha->Write();
    hEnergyDepsNucleus_HitsAlpha->Write();
    hEnergyDepsMembrane_HitsAlpha->Write();
    hEnergyDepsCytoplasm_HitsAlpha->Write();

    hDoseTotalCell_HitsAlpha->Write();
    hDoseNucleus_HitsAlpha->Write();
    hDoseMembrane_HitsAlpha->Write();
    hDoseCytoplasm_HitsAlpha->Write();


    //------------------–----------
    // Hit multiplicity histograms

    hFractionHitsAlpha_TotalCell->Write();
    hFractionHitsAlpha_Nucleus->Write();
    hFractionHitsAlpha_Membrane->Write();
    hFractionHitsAlpha_Cytoplasm->Write();
}

void EnergyDepositionHistograms::ResetHistograms()
{
        //------------------------
    // Energy deposition histograms
    hEnergyDepsMembrane_eVBinning->Reset();
    hEnergyDepsCytoplasm_eVBinning->Reset();
    hEnergyDepsNucleus_eVBinning->Reset();
    hEnergyDepsTotalCell_eVBinning->Reset();

    hEnergyDepsMembrane_keVBinning->Reset();
    hEnergyDepsCytoplasm_keVBinning->Reset();
    hEnergyDepsNucleus_keVBinning->Reset();
    hEnergyDepsTotalCell_keVBinning->Reset();

    hEnergyDepsNucleus_FromSolution->Reset();
    hEnergyDepsNucleus_FromMembrane->Reset();
    hEnergyDepsNucleus_FromCytoplasm->Reset();

    hEnergyDepsMembrane_FromSolution->Reset();
    hEnergyDepsMembrane_FromMembrane->Reset();
    hEnergyDepsMembrane_FromCytoplasm->Reset();

    hEnergyDepsCytoplasm_FromSolution->Reset();
    hEnergyDepsCytoplasm_FromMembrane->Reset();
    hEnergyDepsCytoplasm_FromCytoplasm->Reset();

    hEnergyDepsTotalCell_FromSolution->Reset();
    hEnergyDepsTotalCell_FromMembrane->Reset();
    hEnergyDepsTotalCell_FromCytoplasm->Reset();

    //-----------------------------
    // Dose Histograms

    hDoseMembrane_uGyBinning->Reset();
    hDoseCytoplasm_uGyBinning->Reset();
    hDoseNucleus_uGyBinning->Reset();
    hDoseTotalCell_uGyBinning->Reset();

    hDoseMembrane_mGyBinning->Reset();
    hDoseCytoplasm_mGyBinning->Reset();
    hDoseNucleus_mGyBinning->Reset();
    hDoseTotalCell_mGyBinning->Reset();

    hDoseNucleus_FromSolution->Reset();
    hDoseNucleus_FromMembrane->Reset();
    hDoseNucleus_FromCytoplasm->Reset();

    hDoseMembrane_FromSolution->Reset();
    hDoseMembrane_FromMembrane->Reset();
    hDoseMembrane_FromCytoplasm->Reset();

    hDoseCytoplasm_FromSolution->Reset();
    hDoseCytoplasm_FromMembrane->Reset();
    hDoseCytoplasm_FromCytoplasm->Reset();

    hDoseTotalCell_FromSolution->Reset();
    hDoseTotalCell_FromMembrane->Reset();
    hDoseTotalCell_FromCytoplasm->Reset();



    //--------------------------------
    // Number of times hit vs energy deposition 2D histograms
    hEnergyDepsTotalCell_HitsAlpha->Reset();
    hEnergyDepsNucleus_HitsAlpha->Reset();
    hEnergyDepsMembrane_HitsAlpha->Reset();
    hEnergyDepsCytoplasm_HitsAlpha->Reset();

    hDoseTotalCell_HitsAlpha->Reset();
    hDoseNucleus_HitsAlpha->Reset();
    hDoseMembrane_HitsAlpha->Reset();
    hDoseCytoplasm_HitsAlpha->Reset();


    //------------------–----------
    // Hit multiplicity histograms

    hFractionHitsAlpha_TotalCell->Reset();
    hFractionHitsAlpha_Nucleus->Reset();
    hFractionHitsAlpha_Membrane->Reset();
    hFractionHitsAlpha_Cytoplasm->Reset();
}