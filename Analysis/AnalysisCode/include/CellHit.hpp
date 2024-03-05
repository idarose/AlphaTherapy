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

#ifndef CELLHIT_HPP
#define CELLHIT_HPP

//------------------–----------
class CellHit
{
    //----------------------
    //  Class to store information on the energy depositions in ONE specific cell.
    //  When a new cell is "registered" as hit by a particle, a new instance of this class is created.
    //  Any further energy depositions registered in this cell then gets added to the previously
    //  registered energy deposition. The energy deposition is stored according to which "part" of
    //  the cell is hit by the radiation

    public:
        CellHit(int cellID_in);
        // ~CellHit();

        int GetCellID(){return cellID;};

        double GetMassNucleus(){return massNucleus;};
        double GetMassMembrane(){return massMembrane;};
        double GetMassCytoplasm(){return massCytoplasm;};
        double GetMassCell(){return massCell;};

        std::vector<std::tuple<int,double,int,int>> GetEnergyDepsVec(){return energyDepsVec;};

        double GetEnergyDepositionMembrane(){return energyDepMembrane;};
        double GetEnergyDepositionCytoplasm(){return energyDepCytoplasm;};
        double GetEnergyDepositionNucleus(){return energyDepNucleus;};
        double GetEnergyDepositionTotalCell(){return energyDepTotalCell;};

        double GetEnergyDepositionMembrane_FromSolution(){return energyDepMembrane_FromSolution;};
        double GetEnergyDepositionMembrane_FromMembrane(){return energyDepMembrane_FromMembrane;};
        double GetEnergyDepositionMembrane_FromCytoplasm(){return energyDepMembrane_FromCytoplasm;};

        double GetEnergyDepositionCytoplasm_FromSolution(){return energyDepCytoplasm_FromSolution;};
        double GetEnergyDepositionCytoplasm_FromMembrane(){return energyDepCytoplasm_FromMembrane;};
        double GetEnergyDepositionCytoplasm_FromCytoplasm(){return energyDepCytoplasm_FromCytoplasm;};

        double GetEnergyDepositionNucleus_FromSolution(){return energyDepNucleus_FromSolution;};
        double GetEnergyDepositionNucleus_FromMembrane(){return energyDepNucleus_FromMembrane;};
        double GetEnergyDepositionNucleus_FromCytoplasm(){return energyDepNucleus_FromCytoplasm;};

        double GetEnergyDepositionTotalCell_FromSolution(){return energyDepTotalCell_FromSolution;};
        double GetEnergyDepositionTotalCell_FromMembrane(){return energyDepTotalCell_FromMembrane;};
        double GetEnergyDepositionTotalCell_FromCytoplasm(){return energyDepTotalCell_FromCytoplasm;};

        double GetEnergyDepositionNucleus_FromAlpha(){return energyDepNucleus_FromAlpha;};
        double GetEnergyDepositionMembrane_FromAlpha(){return energyDepMembrane_FromAlpha;};
        double GetEnergyDepositionCytoplasm_FromAlpha(){return energyDepCytoplasm_FromAlpha;};
        double GetEnergyDepositionTotalCell_FromAlpha(){return energyDepTotalCell_FromAlpha;};

        int GetNumberHitsAlphas_TotalCell(){return numberHitsAlphas_TotalCell;};
        int GetNumberHitsAlphas_Nucleus(){return nummberHitsAlphas_Nucleus;};
        int GetNumberHitsAlphas_Membrane(){return numberHitsAlphas_Membrane;};
        int GetNumberHitsAlphas_Cytoplasm(){return numberHitsAlphas_Cytoplasm;};

        int GetNumberHitsBeta(){return numberHitsBetas;};

        void AddEnergyDeposition(double energyDep_in, int volumeTypeInteraction_in, int volumeTypeOriginDecay_in, int particleType_in);
        void HitByAlphaParticle(int volumeTypeHit, bool firstTimeCountingAlpha);
        void HitByBetaParticle();

        void FinalizeCellHit();


    private:
        // Mass cell components
        double massNucleus;
        double massCytoplasm;
        double massCell;
        double massMembrane;

        int cellID;
        int numberHitsAlphas_TotalCell;
        int nummberHitsAlphas_Nucleus;
        int numberHitsAlphas_Membrane;
        int numberHitsAlphas_Cytoplasm;
        int numberHitsAlphas;
        int numberHitsBetas;

        // A vector of tuples for every interaction <volume type of interaction, energydep in that interaction, origin of decay volume type, track ID, parent track ID>
        std::vector<std::tuple<int,double, int, int>> energyDepsVec;

        double energyDepMembrane;
        double energyDepCytoplasm;
        double energyDepNucleus;
        double energyDepTotalCell;


        double energyDepMembrane_FromSolution;
        double energyDepMembrane_FromMembrane;
        double energyDepMembrane_FromCytoplasm;

        double energyDepCytoplasm_FromSolution;
        double energyDepCytoplasm_FromMembrane;
        double energyDepCytoplasm_FromCytoplasm;

        double energyDepNucleus_FromSolution;
        double energyDepNucleus_FromMembrane;
        double energyDepNucleus_FromCytoplasm;

        double energyDepTotalCell_FromSolution;
        double energyDepTotalCell_FromMembrane;
        double energyDepTotalCell_FromCytoplasm;

        double energyDepNucleus_FromAlpha;
        double energyDepMembrane_FromAlpha;
        double energyDepCytoplasm_FromAlpha;
        double energyDepTotalCell_FromAlpha;
};

#endif // CELLHIT_HPP