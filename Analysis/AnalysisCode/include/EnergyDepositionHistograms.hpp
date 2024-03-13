#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TMath.h>
#include <TChain.h>

#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>
#include <ctime>
#include <tuple>
#include <future>
#include <thread>

#ifndef ENERGYDEPOSITIONHISTOGRAMS_HPP
#define ENERGYDEPOSITIONHISTOGRAMS_HPP

class AddEnergyDepositionHistograms;
class DecayDynamics;
class CellHit;


//------------------–----------
class EnergyDepositionHistograms
{
    public:
        EnergyDepositionHistograms(double EMax_in, int iteration);
        // ~EnergyDepositionHistograms();

        void GenerateEmptyHistograms(DecayDynamics decayDynamicsInstance);
        void AddCellHitsToHistograms(CellHit cellHit);
        void ScaleHistograms(double numberCells_in, double numberIterations_in);
        void ScaleHistogramKineticEnergyAlphas_PerIteration(double numberHitsAlphasTotalCell_ThisIteration, double numberHitsAlphasNucleus_ThisIteration);
        void WriteHistogramsToFile();
        void ResetHistograms();

        friend class AddEnergyDepositionHistograms;

    private:
        double EMax;
        int it;

        //----------------------------
        // Energy deposition histograms

        // Histogram for total energy deposited in one nuclei per number of cells
        TH1D *hEnergyDepsNucleus_eVBinning;
        TH1D *hEnergyDepsNucleus_keVBinning;
        TH1D *hEnergyDepsNucleus_FromSolution;
        TH1D *hEnergyDepsNucleus_FromMembrane;
        TH1D *hEnergyDepsNucleus_FromCytoplasm;

        // Histogram for total energy deposited in one membrane per number of cells
        TH1D *hEnergyDepsMembrane_eVBinning;
        TH1D *hEnergyDepsMembrane_keVBinning;
        TH1D *hEnergyDepsMembrane_FromSolution;
        TH1D *hEnergyDepsMembrane_FromMembrane;
        TH1D *hEnergyDepsMembrane_FromCytoplasm;

        // Histogram for total energy deposited in one cytoplasm per number of cells
        TH1D *hEnergyDepsCytoplasm_eVBinning;
        TH1D *hEnergyDepsCytoplasm_keVBinning;
        TH1D *hEnergyDepsCytoplasm_FromSolution;
        TH1D *hEnergyDepsCytoplasm_FromMembrane;
        TH1D *hEnergyDepsCytoplasm_FromCytoplasm;

        // Histogram for total energy deposited in one cell per number of cells
        TH1D *hEnergyDepsTotalCell_eVBinning;
        TH1D *hEnergyDepsTotalCell_keVBinning;
        TH1D *hEnergyDepsTotalCell_FromSolution;
        TH1D *hEnergyDepsTotalCell_FromMembrane;
        TH1D *hEnergyDepsTotalCell_FromCytoplasm;

        // Histograms for number of hits by alpha particles
        TH2D *hEnergyDepsTotalCell_HitsAlpha;
        TH2D *hEnergyDepsNucleus_HitsAlpha;
        TH2D *hEnergyDepsMembrane_HitsAlpha;
        TH2D *hEnergyDepsCytoplasm_HitsAlpha;

        //--------------------------
        // Dose histograms

        // Histogram for total energy deposited in one nuclei per number of cells
        TH1D *hDoseNucleus_mGyBinning;
        TH1D *hDoseNucleus_uGyBinning;
        TH1D *hDoseNucleus_FromSolution;
        TH1D *hDoseNucleus_FromMembrane;
        TH1D *hDoseNucleus_FromCytoplasm;

        // Histogram for total energy deposited in one membrane per number of cells
        TH1D *hDoseMembrane_mGyBinning;
        TH1D *hDoseMembrane_uGyBinning;
        TH1D *hDoseMembrane_FromSolution;
        TH1D *hDoseMembrane_FromMembrane;
        TH1D *hDoseMembrane_FromCytoplasm;


        // Histogram for total energy deposited in one cytoplasm per number of cells
        TH1D *hDoseCytoplasm_mGyBinning;
        TH1D *hDoseCytoplasm_uGyBinning;
        TH1D *hDoseCytoplasm_FromSolution;
        TH1D *hDoseCytoplasm_FromMembrane;
        TH1D *hDoseCytoplasm_FromCytoplasm;

        // Histogram for total energy deposited in one cell per number of cells
        TH1D *hDoseTotalCell_mGyBinning;
        TH1D *hDoseTotalCell_uGyBinning;
        TH1D *hDoseTotalCell_FromSolution;
        TH1D *hDoseTotalCell_FromMembrane;
        TH1D *hDoseTotalCell_FromCytoplasm;

        // Histograms for dose per number of hits by alpha particles
        TH2D* hDoseTotalCell_HitsAlpha;
        TH2D* hDoseNucleus_HitsAlpha;
        TH2D* hDoseMembrane_HitsAlpha;
        TH2D* hDoseCytoplasm_HitsAlpha;

        //-------------------------
        // Alpha particle hits histograms

        //--------------------------
        // Hit multiplicity histograms
        TH1D *hFractionHitsAlpha_TotalCell;
        TH1D *hFractionHitsAlpha_Nucleus;
        TH1D *hFractionHitsAlpha_Membrane;
        TH1D *hFractionHitsAlpha_Cytoplasm;

        //-----------------------
        // Histrogram for kinetic energy of alpha particle before entering a volume
        TH1D* hKineticEnergyAlphaTotalCell_FromSolution;
        TH1D* hKineticEnergyAlphaTotalCell_FromMembrane;
        TH1D* hKineticEnergyAlphaTotalCell_FromCytoplasm;

        TH1D* hKineticEnergyAlphaNucleus_FromSolution;
        TH1D* hKineticEnergyAlphaNucleus_FromMembrane;
        TH1D* hKineticEnergyAlphaNucleus_FromCytoplasm;
};

#endif // ENERGYDEPOSITIONHISTOGRAMS_HPP