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
        std::vector<std::tuple<int,double,int,int>> GetEnergyDepsVec(){return energyDepsVec;};

        //---------------------------
        void AddEnergyDeposition(double energyDep_in, int volumeTypeInteraction_in, int volumeTypeOriginDecay_in, int particleType_in);
        void HitByAlphaParticle(int volumeTypeHit, bool firstTimeCountingAlpha, int volumeTypeOriginDecay_in, double kineticEnergyAlpha);
        // void HitByBetaParticle();
        void FinalizeCellHit(std::string cellGeometry);

        //---------------------------
        double GetMassNucleus(){return massNucleus;};
        double GetMassMembrane(){return massMembrane;};
        double GetMassCytoplasm(){return massCytoplasm;};
        double GetMassCell(){return massCell;};

        //---------------------------
        double GetEnergyDepositionMembrane(){return energyDepMembrane;};
        double GetEnergyDepositionCytoplasm(){return energyDepCytoplasm;};
        double GetEnergyDepositionNucleus(){return energyDepNucleus;};
        double GetEnergyDepositionTotalCell(){return energyDepTotalCell;};

        //---------------------------
        double GetEnergyDepositionMembrane_FromSolution(){return energyDepMembrane_FromSolution;};
        double GetEnergyDepositionMembrane_FromMembrane(){return energyDepMembrane_FromMembrane;};
        double GetEnergyDepositionMembrane_FromCytoplasm(){return energyDepMembrane_FromCytoplasm;};

        double GetFractionEnergyDepMembrane_FromSolution(){return fractionEnergyDepMembrane_FromSolution;}
        double GetFractionEnergyDepMembrane_FromMembrane(){return fractionEnergyDepMembrane_FromMembrane;}
        double GetFractionEnergyDepMembrane_FromCytoplasm(){return fractionEnergyDepMembrane_FromCytoplasm;}

        //---------------------------
        double GetEnergyDepositionCytoplasm_FromSolution(){return energyDepCytoplasm_FromSolution;};
        double GetEnergyDepositionCytoplasm_FromMembrane(){return energyDepCytoplasm_FromMembrane;};
        double GetEnergyDepositionCytoplasm_FromCytoplasm(){return energyDepCytoplasm_FromCytoplasm;};

        double GetFractionEnergyDepCytoplasm_FromSolution(){return fractionEnergyDepCytoplasm_FromSolution;}
        double GetFractionEnergyDepCytoplasm_FromMembrane(){return fractionEnergyDepCytoplasm_FromMembrane;}
        double GetFractionEnergyDepCytoplasm_FromCytoplasm(){return fractionEnergyDepCytoplasm_FromCytoplasm;}

        //---------------------------
        double GetEnergyDepositionNucleus_FromSolution(){return energyDepNucleus_FromSolution;};
        double GetEnergyDepositionNucleus_FromMembrane(){return energyDepNucleus_FromMembrane;};
        double GetEnergyDepositionNucleus_FromCytoplasm(){return energyDepNucleus_FromCytoplasm;};

        double GetFractionEnergyDepNucleus_FromSolution(){return fractionEnergyDepNucleus_FromSolution;}
        double GetFractionEnergyDepNucleus_FromMembrane(){return fractionEnergyDepNucleus_FromMembrane;}
        double GetFractionEnergyDepNucleus_FromCytoplasm(){return fractionEnergyDepNucleus_FromCytoplasm;}

        //---------------------------
        double GetEnergyDepositionTotalCell_FromSolution(){return energyDepTotalCell_FromSolution;};
        double GetEnergyDepositionTotalCell_FromMembrane(){return energyDepTotalCell_FromMembrane;};
        double GetEnergyDepositionTotalCell_FromCytoplasm(){return energyDepTotalCell_FromCytoplasm;};

        double GetFractionEnergyDepTotalCell_FromSolution(){return fractionEnergyDepTotalCell_FromSolution;}
        double GetFractionEnergyDepTotalCell_FromMembrane(){return fractionEnergyDepTotalCell_FromMembrane;}
        double GetFractionEnergyDepTotalCell_FromCytoplasm(){return fractionEnergyDepTotalCell_FromCytoplasm;}

        //---------------------------
        int GetNumberHitsAlphas_TotalCell(){return numberHitsAlphas_TotalCell;};
        int GetNumberHitsAlphas_Nucleus(){return numberHitsAlphas_Nucleus;};
        int GetNumberHitsAlphas_Membrane(){return numberHitsAlphas_Membrane;};
        int GetNumberHitsAlphas_Cytoplasm(){return numberHitsAlphas_Cytoplasm;};


        int GetNumberHitsAlphasTotalCell_FromSolution(){return numberHitsAlphasTotalCell_FromSolution;};
        int GetNumberHitsAlphasTotalCell_FromMembrane(){return numberHitsAlphasTotalCell_FromMembrane;};
        int GetNumberHitsAlphasTotalCell_FromCytoplasm(){return numberHitsAlphasTotalCell_FromCytoplasm;};

        int GetNumberHitsAlphasNucleus_FromSolution(){return numberHitsAlphasNucleus_FromSolution;};
        int GetNumberHitsAlphasNucleus_FromMembrane(){return numberHitsAlphasNucleus_FromMembrane;};
        int GetNumberHitsAlphasNucleus_FromCytoplasm(){return numberHitsAlphasNucleus_FromCytoplasm;};


        //---------------------------
        std::vector<double> GetKineticEnergyAlphaTotalCell_FromSolution_Vec(){return kineticEnergyAlphaTotalCell_FromSolution_Vec;};
        std::vector<double> GetKineticEnergyAlphaTotalCell_FromMembrane_Vec(){return kineticEnergyAlphaTotalCell_FromMembrane_Vec;};
        std::vector<double> GetKineticEnergyAlphaTotalCell_FromCytoplasm_Vec(){return kineticEnergyAlphaTotalCell_FromCytoplasm_Vec;};

        std::vector<double> GetKineticEnergyAlphaNucleus_FromSolution_Vec(){return kineticEnergyAlphaNucleus_FromSolution_Vec;};
        std::vector<double> GetKineticEnergyAlphaNucleus_FromMembrane_Vec(){return kineticEnergyAlphaNucleus_FromMembrane_Vec;};
        std::vector<double> GetKineticEnergyAlphaNucleus_FromCytoplasm_Vec(){return kineticEnergyAlphaNucleus_FromCytoplasm_Vec;};



    private:

        int cellID;

        // Mass cell components
        double massNucleus;
        double massCytoplasm;
        double massCell;
        double massMembrane;

        int numberHitsAlphas_TotalCell;
        int numberHitsAlphas_Nucleus;
        int numberHitsAlphas_Membrane;
        int numberHitsAlphas_Cytoplasm;

        int numberHitsAlphasTotalCell_FromSolution;
        int numberHitsAlphasTotalCell_FromMembrane;
        int numberHitsAlphasTotalCell_FromCytoplasm;

        int numberHitsAlphasNucleus_FromSolution;
        int numberHitsAlphasNucleus_FromMembrane;
        int numberHitsAlphasNucleus_FromCytoplasm;

        // A vector of tuples for every interaction <volume type of interaction, energydep in that interaction, origin of decay volume type, track ID, parent track ID>
        std::vector<std::tuple<int,double, int, int>> energyDepsVec;

        double energyDepMembrane;
        double energyDepCytoplasm;
        double energyDepNucleus;
        double energyDepTotalCell;

        //-------------------------
        double energyDepMembrane_FromSolution;
        double energyDepMembrane_FromMembrane;
        double energyDepMembrane_FromCytoplasm;

        double fractionEnergyDepMembrane_FromSolution;
        double fractionEnergyDepMembrane_FromMembrane;
        double fractionEnergyDepMembrane_FromCytoplasm;

        //-------------------------
        double energyDepCytoplasm_FromSolution;
        double energyDepCytoplasm_FromMembrane;
        double energyDepCytoplasm_FromCytoplasm;

        double fractionEnergyDepCytoplasm_FromSolution;
        double fractionEnergyDepCytoplasm_FromMembrane;
        double fractionEnergyDepCytoplasm_FromCytoplasm;

        //-------------------------
        double energyDepNucleus_FromSolution;
        double energyDepNucleus_FromMembrane;
        double energyDepNucleus_FromCytoplasm;

        double fractionEnergyDepNucleus_FromSolution;
        double fractionEnergyDepNucleus_FromMembrane;
        double fractionEnergyDepNucleus_FromCytoplasm;

        //-------------------------
        double energyDepTotalCell_FromSolution;
        double energyDepTotalCell_FromMembrane;
        double energyDepTotalCell_FromCytoplasm;

        double fractionEnergyDepTotalCell_FromSolution;
        double fractionEnergyDepTotalCell_FromMembrane;
        double fractionEnergyDepTotalCell_FromCytoplasm;

        //-------------------------
        std::vector<double> kineticEnergyAlphaTotalCell_FromSolution_Vec;
        std::vector<double> kineticEnergyAlphaTotalCell_FromMembrane_Vec;
        std::vector<double> kineticEnergyAlphaTotalCell_FromCytoplasm_Vec;

        std::vector<double> kineticEnergyAlphaNucleus_FromSolution_Vec;
        std::vector<double> kineticEnergyAlphaNucleus_FromMembrane_Vec;
        std::vector<double> kineticEnergyAlphaNucleus_FromCytoplasm_Vec;
};

#endif // CELLHIT_HPP