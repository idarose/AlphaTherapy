#include "../include/EnergyDepositionHistograms.hpp"
#include "../include/AddEnergyDepositionHistograms.hpp"

void AddEnergyDepositionHistograms::AddHistograms(EnergyDepositionHistograms& h1, const EnergyDepositionHistograms& h2)
{
         // Histogram for total energy deposited in one nuclei per number of cells
        h1.hEnergyDepsNucleus_eVBinning->Add(h2.hEnergyDepsNucleus_eVBinning);
        h1.hEnergyDepsNucleus_keVBinning->Add(h2.hEnergyDepsNucleus_keVBinning);
        h1.hEnergyDepsNucleus_FromSolution->Add(h2.hEnergyDepsNucleus_FromSolution);
        h1.hEnergyDepsNucleus_FromMembrane->Add(h2.hEnergyDepsNucleus_FromMembrane);
        h1.hEnergyDepsNucleus_FromCytoplasm->Add(h2.hEnergyDepsNucleus_FromCytoplasm);

        h1.hEnergyDepsNucleus_FractionFromSolution->Add(h2.hEnergyDepsNucleus_FractionFromSolution);
        h1.hEnergyDepsNucleus_FractionFromMembrane->Add(h2.hEnergyDepsNucleus_FractionFromMembrane);
        h1.hEnergyDepsNucleus_FractionFromCytoplasm->Add(h2.hEnergyDepsNucleus_FractionFromCytoplasm);

        // Histogram for total energy deposited in one membrane per number of cells
        h1.hEnergyDepsMembrane_eVBinning->Add(h2.hEnergyDepsMembrane_eVBinning);
        h1.hEnergyDepsMembrane_keVBinning->Add(h2.hEnergyDepsMembrane_keVBinning);
        h1.hEnergyDepsMembrane_FromSolution->Add(h2.hEnergyDepsMembrane_FromSolution);
        h1.hEnergyDepsMembrane_FromMembrane->Add(h2.hEnergyDepsMembrane_FromMembrane);
        h1.hEnergyDepsMembrane_FromCytoplasm->Add(h2.hEnergyDepsMembrane_FromCytoplasm);

        h1.hEnergyDepsMembrane_FractionFromSolution->Add(h2.hEnergyDepsMembrane_FractionFromSolution);
        h1.hEnergyDepsMembrane_FractionFromMembrane->Add(h2.hEnergyDepsMembrane_FractionFromMembrane);
        h1.hEnergyDepsMembrane_FractionFromCytoplasm->Add(h2.hEnergyDepsMembrane_FractionFromCytoplasm);


        // Histogram for total energy deposited in one cytoplasm per number of cells
        h1.hEnergyDepsCytoplasm_eVBinning->Add(h2.hEnergyDepsCytoplasm_eVBinning);
        h1.hEnergyDepsCytoplasm_keVBinning->Add(h2.hEnergyDepsCytoplasm_keVBinning);
        h1.hEnergyDepsCytoplasm_FromSolution->Add(h2.hEnergyDepsCytoplasm_FromSolution);
        h1.hEnergyDepsCytoplasm_FromMembrane->Add(h2.hEnergyDepsCytoplasm_FromMembrane);
        h1.hEnergyDepsCytoplasm_FromCytoplasm->Add(h2.hEnergyDepsCytoplasm_FromCytoplasm);

        h1.hEnergyDepsCytoplasm_FractionFromSolution->Add(h2.hEnergyDepsCytoplasm_FractionFromSolution);
        h1.hEnergyDepsCytoplasm_FractionFromMembrane->Add(h2.hEnergyDepsCytoplasm_FractionFromMembrane);
        h1.hEnergyDepsCytoplasm_FractionFromCytoplasm->Add(h2.hEnergyDepsCytoplasm_FractionFromCytoplasm);

        // Histogram for total energy deposited in one cell per number of cells
        h1.hEnergyDepsTotalCell_eVBinning->Add(h2.hEnergyDepsTotalCell_eVBinning);
        h1.hEnergyDepsTotalCell_keVBinning->Add(h2.hEnergyDepsTotalCell_keVBinning);
        h1.hEnergyDepsTotalCell_FromSolution->Add(h2.hEnergyDepsTotalCell_FromSolution);
        h1.hEnergyDepsTotalCell_FromMembrane->Add(h2.hEnergyDepsTotalCell_FromMembrane);
        h1.hEnergyDepsTotalCell_FromCytoplasm->Add(h2.hEnergyDepsTotalCell_FromCytoplasm);

        h1.hEnergyDepsTotalCell_FractionFromSolution->Add(h2.hEnergyDepsTotalCell_FractionFromSolution);
        h1.hEnergyDepsTotalCell_FractionFromMembrane->Add(h2.hEnergyDepsTotalCell_FractionFromMembrane);
        h1.hEnergyDepsTotalCell_FractionFromCytoplasm->Add(h2.hEnergyDepsTotalCell_FractionFromCytoplasm);

        // Histograms for number of hits by alpha particles
        h1.hEnergyDepsTotalCell_HitsAlpha->Add(h2.hEnergyDepsTotalCell_HitsAlpha);
        h1.hEnergyDepsNucleus_HitsAlpha->Add(h2.hEnergyDepsNucleus_HitsAlpha);
        h1.hEnergyDepsMembrane_HitsAlpha->Add(h2.hEnergyDepsMembrane_HitsAlpha);
        h1.hEnergyDepsCytoplasm_HitsAlpha->Add(h2.hEnergyDepsCytoplasm_HitsAlpha);

        // General hit histogram
        h1.hFractionHitsAlpha_TotalCell->Add(h2.hFractionHitsAlpha_TotalCell);
        h1.hFractionHitsAlpha_Nucleus->Add(h2.hFractionHitsAlpha_Nucleus);
        h1.hFractionHitsAlpha_Membrane->Add(h2.hFractionHitsAlpha_Membrane);
        h1.hFractionHitsAlpha_Cytoplasm->Add(h2.hFractionHitsAlpha_Cytoplasm);

        //--------------------------
        // Dose histograms

        // Histogram for total energy deposited in one nuclei per number of cells
        h1.hDoseNucleus_uGyBinning->Add(h2.hDoseNucleus_uGyBinning);
        h1.hDoseNucleus_mGyBinning->Add(h2.hDoseNucleus_mGyBinning);
        h1.hDoseNucleus_FromSolution->Add(h2.hDoseNucleus_FromSolution);
        h1.hDoseNucleus_FromMembrane->Add(h2.hDoseNucleus_FromMembrane);
        h1.hDoseNucleus_FromCytoplasm->Add(h2.hDoseNucleus_FromCytoplasm);

        h1.hDoseNucleus_FractionFromSolution->Add(h2.hDoseNucleus_FractionFromSolution);
        h1.hDoseNucleus_FractionFromMembrane->Add(h2.hDoseNucleus_FractionFromMembrane);
        h1.hDoseNucleus_FractionFromCytoplasm->Add(h2.hDoseNucleus_FractionFromCytoplasm);

        // Histogram for total energy deposited in one membrane per number of cells
        h1.hDoseMembrane_uGyBinning->Add(h2.hDoseMembrane_uGyBinning);
        h1.hDoseMembrane_mGyBinning->Add(h2.hDoseMembrane_mGyBinning);
        h1.hDoseMembrane_FromSolution->Add(h2.hDoseMembrane_FromSolution);
        h1.hDoseMembrane_FromMembrane->Add(h2.hDoseMembrane_FromMembrane);
        h1.hDoseMembrane_FromCytoplasm->Add(h2.hDoseMembrane_FromCytoplasm);

        h1.hDoseMembrane_FractionFromSolution->Add(h2.hDoseMembrane_FractionFromSolution);
        h1.hDoseMembrane_FractionFromMembrane->Add(h2.hDoseMembrane_FractionFromMembrane);
        h1.hDoseMembrane_FractionFromCytoplasm->Add(h2.hDoseMembrane_FractionFromCytoplasm);


        // Histogram for total energy deposited in one cytoplasm per number of cells
        h1.hDoseCytoplasm_uGyBinning->Add(h2.hDoseCytoplasm_uGyBinning);
        h1.hDoseCytoplasm_mGyBinning->Add(h2.hDoseCytoplasm_mGyBinning);
        h1.hDoseCytoplasm_FromSolution->Add(h2.hDoseCytoplasm_FromSolution);
        h1.hDoseCytoplasm_FromMembrane->Add(h2.hDoseCytoplasm_FromMembrane);
        h1.hDoseCytoplasm_FromCytoplasm->Add(h2.hDoseCytoplasm_FromCytoplasm);

        h1.hDoseCytoplasm_FractionFromSolution->Add(h2.hDoseCytoplasm_FractionFromSolution);
        h1.hDoseCytoplasm_FractionFromMembrane->Add(h2.hDoseCytoplasm_FractionFromMembrane);
        h1.hDoseCytoplasm_FractionFromCytoplasm->Add(h2.hDoseCytoplasm_FractionFromCytoplasm);

        // Histogram for total energy deposited in one cell per number of cells
        h1.hDoseTotalCell_uGyBinning->Add(h2.hDoseTotalCell_uGyBinning);
        h1.hDoseTotalCell_mGyBinning->Add(h2.hDoseTotalCell_mGyBinning);
        h1.hDoseTotalCell_FromSolution->Add(h2.hDoseTotalCell_FromSolution);
        h1.hDoseTotalCell_FromMembrane->Add(h2.hDoseTotalCell_FromMembrane);
        h1.hDoseTotalCell_FromCytoplasm->Add(h2.hDoseTotalCell_FromCytoplasm);

        h1.hDoseTotalCell_FractionFromSolution->Add(h2.hDoseTotalCell_FractionFromSolution);
        h1.hDoseTotalCell_FractionFromMembrane->Add(h2.hDoseTotalCell_FractionFromMembrane);
        h1.hDoseTotalCell_FractionFromCytoplasm->Add(h2.hDoseTotalCell_FractionFromCytoplasm);

        //--------------------------
        // Hits by alpha alpha particle

        // Histograms for number of hits by alpha particles
        h1.hDoseTotalCell_HitsAlpha->Add(h2.hDoseTotalCell_HitsAlpha);
        h1.hDoseNucleus_HitsAlpha->Add(h2.hDoseNucleus_HitsAlpha);
        h1.hDoseMembrane_HitsAlpha->Add(h2.hDoseMembrane_HitsAlpha);
        h1.hDoseCytoplasm_HitsAlpha->Add(h2.hDoseCytoplasm_HitsAlpha);


        //-----------------------
        // Histrogram for kinetic energy of alpha particle before entering a volume
        h1.hKineticEnergyAlphaTotalCell_FromSolution->Add(h2.hKineticEnergyAlphaTotalCell_FromSolution);
        h1.hKineticEnergyAlphaTotalCell_FromMembrane->Add(h2.hKineticEnergyAlphaTotalCell_FromMembrane);
        h1.hKineticEnergyAlphaTotalCell_FromCytoplasm->Add(h2.hKineticEnergyAlphaTotalCell_FromCytoplasm);

        h1.hKineticEnergyAlphaNucleus_FromSolution->Add(h2.hKineticEnergyAlphaNucleus_FromSolution);
        h1.hKineticEnergyAlphaNucleus_FromMembrane->Add(h2.hKineticEnergyAlphaNucleus_FromMembrane);
        h1.hKineticEnergyAlphaNucleus_FromCytoplasm->Add(h2.hKineticEnergyAlphaNucleus_FromCytoplasm);

}
