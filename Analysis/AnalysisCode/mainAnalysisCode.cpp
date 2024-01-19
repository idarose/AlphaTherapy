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


//------------------–----------
std::vector<double> ImportDataStringDouble(std::string filename)
{
    //----------------------
    //  Imports data from file where the first column contains data in string form
    //  and second column contains data in double form. Returns the columns of double
    //  values in the form of a vector

    std::vector<double> input_data;

    std::fstream myfile(filename, ios_base::in);


    if (myfile.is_open())
    {
        std::string line;
        std::string x;
        double y;

        while(std::getline(myfile,line))
        {
            std::stringstream mystream(line);
            mystream >> x >> y;
            input_data.push_back(y);
        }
    }
    else
    {
        std::cout << "Unable to open file " << filename << std::endl;
    }
    myfile.close();

    return input_data;

}

//------------------–----------
class DecayDynamics
{
    //----------------------
    //  Class to store calculated decay dynamics data from the Mathematica calculations.
    //  Assumes the output file from Mathematica is structured in a specific way

    public:
        DecayDynamics(int activitySample_in, double U0PerCell_in, double U0InternalizedPerCell_in, std::string cellLine_in);

        void LoadDataFromMathematicaCalculations(std::string filepathToMathematicaOutput);

        double GetNumberDecaysInSolutionFirstHour(){return numberDecays212PbInSolution1hTo2h;};
        double GetNumberDecaysInMembraneTotalTime(){return numberDecays212PbInMembrane1h2To26h;};
        double GetNumberDecaysInCytoplasmTotalTime(){return numberDecays212PbInCytoplasm1hTo26h;};

        int GetActivity(){return activitySample;};
        std::string GetCellLine(){return cellLine;};

    private:

        // These number of decays are calculated in Mathematica for a sample of 0.2mL in volume
        double numberDecays212PbInSolution1hTo2h;
        double numberDecays212PbInMembrane1h2To26h;
        double numberDecays212PbInCytoplasm1hTo26h;

        double U0PerCell; // Number radionuclides absorbed per cell
        double U0InternalizedPerCell; // Number of absorbed radionuclides that are found in cytoplasm
        int activitySample;  // Given in kBq/1mL
        std::string cellLine;
};


DecayDynamics::DecayDynamics(int activitySample_in, double U0PerCell_in, double U0InternalizedPerCell_in, std::string cellLine_in)
{
    activitySample = activitySample_in;
    U0PerCell = U0PerCell_in;
    U0InternalizedPerCell = U0InternalizedPerCell_in;
    cellLine = cellLine_in;

}



void DecayDynamics::LoadDataFromMathematicaCalculations(std::string filepathToMathematicaOutput)
{
    //----------------------
    // Loads data from calculations, assuming output files are structured as shown in the file "212PbDecayDynamics.nb"

    std::string filepathSolutionData = filepathToMathematicaOutput + "/" + cellLine + "/Solution/Activity_" + std::to_string(activitySample) + "kBq/NumberDecays.dat";

    // Importing data for decays occuring in solution
    std::vector<double> decayDataSolution = ImportDataStringDouble(filepathSolutionData);
    numberDecays212PbInSolution1hTo2h = decayDataSolution[0];


    // If there are no radionuclides internalized there is no output file for "Cells"
    // So only read "Cells" files if there is uptake
    if(U0PerCell>0.0)
    {
        // Importing data for decays occuring in cells
        std::string filepathCellData = filepathToMathematicaOutput + "/" + cellLine + "/Cells/Activity_" + std::to_string(activitySample) + "kBq/NumberDecays.dat";
        std::vector<double> decayDataCells = ImportDataStringDouble(filepathCellData);

        numberDecays212PbInMembrane1h2To26h = decayDataCells[7];
        numberDecays212PbInCytoplasm1hTo26h = decayDataCells[14];
    }
    else
    {
        numberDecays212PbInMembrane1h2To26h = 0.0;
        numberDecays212PbInCytoplasm1hTo26h = 0.0;
    }
}



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

        int GetCellID(){return cellID;};
        std::vector<std::tuple<int,double,int>> GetEnergyDepsVec(){return energyDepsVec;};
        // std::vector<std::tuple<int,int,int>> GetParticleHitsVec(){return particleHitsVec;};

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

        int GetNumberHitsAlpha(){return numberHitsAlphas;};
        int GetNumberHitsBeta(){return numberHitsBetas;};

        void AddEnergyDeposition(double energyDep_in, int volumeTypeInteraction_in, int volumeTypeOriginDecay_in);
        // void AddParticleHitsInfo(int trackID_in, int parentID_in, int pdgEncoding_in);
        void FinalizeCellHit();

        void HitByAlphaParticle();
        void HitByBetaParticle();


    private:
        int cellID;
        int numberHitsAlphas;
        int numberHitsBetas;

        // A vector of tuples for every interaction <volume type of interaction, energydep in that interaction, origin of decay volume type, track ID, parent track ID>
        std::vector<std::tuple<int,double, int>> energyDepsVec;
        // std::vector<std::tuple<int,int,int>> particleHitsVec;

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
};


//------------------–----------
CellHit::CellHit(int cellID_in)
{
    cellID = cellID_in;
    numberHitsAlphas = 0;

    energyDepMembrane = 0.0;
    energyDepCytoplasm = 0.0;
    energyDepNucleus = 0.0;

    energyDepMembrane_FromSolution = 0.0;
    energyDepMembrane_FromMembrane = 0.0;
    energyDepMembrane_FromCytoplasm = 0.0;

    energyDepCytoplasm_FromSolution = 0.0;
    energyDepCytoplasm_FromMembrane = 0.0;
    energyDepCytoplasm_FromCytoplasm = 0.0;

    energyDepNucleus_FromSolution = 0.0;
    energyDepNucleus_FromMembrane = 0.0;
    energyDepNucleus_FromCytoplasm = 0.0;

    energyDepTotalCell_FromSolution = 0.0;
    energyDepTotalCell_FromMembrane = 0.0;
    energyDepTotalCell_FromCytoplasm = 0.0;

    numberHitsAlphas = 0;
    numberHitsBetas = 0;
}


//------------------–----------
void CellHit::AddEnergyDeposition(double energyDep_in, int volumeTypeInteraction_in, int volumeTypeOriginDecay_in)
{
    std::tuple<int,double,int> hitTuple;
    hitTuple = std::make_tuple(volumeTypeInteraction_in, energyDep_in, volumeTypeOriginDecay_in);
    energyDepsVec.push_back(hitTuple);
}


//------------------–----------
void CellHit::HitByAlphaParticle()
{
    numberHitsAlphas ++;
}

//----------------------------
void CellHit::HitByBetaParticle()
{
    numberHitsBetas ++;
}


//------------------–----------
void CellHit::FinalizeCellHit()
{

    //----------------------
    // Sum the energy depositions per cell component

    int interactionVolume;
    int originVolume;
    double energyDepInteraction;

    for(int i=0; i<energyDepsVec.size(); i++)
    {
        interactionVolume = std::get<0>(energyDepsVec[i]);
        energyDepInteraction = std::get<1>(energyDepsVec[i]);
        originVolume = std::get<2>(energyDepsVec[i]);


        if(interactionVolume==1)
        {
            // std::cout << "In membrane: " << energyDepInteraction << std::endl;
            energyDepMembrane += energyDepInteraction;
            if(originVolume == 0){energyDepMembrane_FromSolution += energyDepInteraction;}
            if(originVolume == 1){energyDepMembrane_FromMembrane += energyDepInteraction;}
            if(originVolume == 2){energyDepMembrane_FromCytoplasm += energyDepInteraction;}
        }
        if(interactionVolume==2)
        {
            std::cout << "In cytoplasm : " << energyDepInteraction << std::endl;
            energyDepCytoplasm += energyDepInteraction;
            if(originVolume == 0){energyDepCytoplasm_FromSolution += energyDepInteraction;}
            if(originVolume == 1){energyDepCytoplasm_FromMembrane += energyDepInteraction;}
            if(originVolume == 2){energyDepCytoplasm_FromCytoplasm += energyDepInteraction;}
        }
        if(interactionVolume==3)
        {
            std::cout << "In nucleus : " << energyDepInteraction << std::endl;
            energyDepNucleus += energyDepInteraction;
            if(originVolume == 0){energyDepNucleus_FromSolution += energyDepInteraction;}
            if(originVolume == 1){energyDepNucleus_FromMembrane += energyDepInteraction;}
            if(originVolume == 2){energyDepNucleus_FromCytoplasm += energyDepInteraction;}
        }

    }

    energyDepTotalCell = energyDepMembrane + energyDepCytoplasm + energyDepNucleus;

    energyDepTotalCell_FromSolution = energyDepMembrane_FromSolution + energyDepCytoplasm_FromSolution + energyDepNucleus_FromSolution;
    energyDepTotalCell_FromMembrane = energyDepMembrane_FromMembrane + energyDepCytoplasm_FromMembrane + energyDepNucleus_FromMembrane;
    energyDepTotalCell_FromCytoplasm = energyDepMembrane_FromCytoplasm + energyDepCytoplasm_FromCytoplasm + energyDepNucleus_FromCytoplasm;

}


//------------------–----------
class EnergyDepositionHistograms
{
    public:
        EnergyDepositionHistograms(int NBins_in, double EMin_in, double EMax_in);

        void GenerateEmptyHistograms(DecayDynamics decayDynamicsInstance);
        void AddCellHitsToHistograms(CellHit cellHit);
        void ScaleHistograms(double factor);
        void WriteHistogramsToFile();

        TH1D* GetEnergyDepNucleusHist(){return hEnergyDepsNucleus;};
        TH1D* GetEnergyDepMembraneHist(){return hEnergyDepsMembrane;};
        TH1D* GetEnergyDepCytoplasmHist(){return hEnergyDepsCytoplasm;};
        TH1D* GetEnergyDepTotalCellHist(){return hEnergyDepsTotalCell;};

    private:
        int NBins;
        double EMin;
        double EMax;


        // Histogram for total energy deposited in one nuclei per number of cells
        TH1D *hEnergyDepsNucleus;
        TH1D *hEnergyDepsNucleus_FromSolution;
        TH1D *hEnergyDepsNucleus_FromMembrane;
        TH1D *hEnergyDepsNucleus_FromCytoplasm;

        // Histogram for total energy deposited in one membrane per number of cells
        TH1D *hEnergyDepsMembrane;
        TH1D *hEnergyDepsMembrane_FromSolution;
        TH1D *hEnergyDepsMembrane_FromMembrane;
        TH1D *hEnergyDepsMembrane_FromCytoplasm;


        // Histogram for total energy deposited in one cytoplasm per number of cells
        TH1D *hEnergyDepsCytoplasm;
        TH1D *hEnergyDepsCytoplasm_FromSolution;
        TH1D *hEnergyDepsCytoplasm_FromMembrane;
        TH1D *hEnergyDepsCytoplasm_FromCytoplasm;

        // Histogram for total energy deposited in one cell per number of cells
        TH1D *hEnergyDepsTotalCell;
        TH1D *hEnergyDepsTotalCell_FromSolution;
        TH1D *hEnergyDepsTotalCell_FromMembrane;
        TH1D *hEnergyDepsTotalCell_FromCytoplasm;

        TH2D* hEnergyDepNucleus_HitsAlpha;


        // // Histogram for total energy deposited in both membrane and cytoplasm of one cell, per number of cells
        // TH1D *hEnergyDepsMembraneAndCytoplasm;

        // // Histogram for total energy deposited in both membrane and nucleus of on cell per number of cells
        // TH1D *hEnergyDepsMembraneAndNucleus;

        // // Histogram for total energy deposited in both nucleus and cytoplasm of one cell per number of cells
        // TH1D *hEnergyDepsNucleusAndCytoplasm;
};



//------------------–----------
EnergyDepositionHistograms::EnergyDepositionHistograms(int NBins_in, double EMin_in, double EMax_in)
{
    NBins = NBins_in;
    EMin = EMin_in;
    EMax = EMax_in;
}


//------------------–----------
void EnergyDepositionHistograms::GenerateEmptyHistograms(DecayDynamics decayDynamicsInstance)
{
    std::string generalHistogramName = "hEnergyDeps_212Pb_" + decayDynamicsInstance.GetCellLine() + "_" + std::to_string(decayDynamicsInstance.GetActivity()) + "kBq_";

    std::string histogramNameNucleus = generalHistogramName + "Nucleus";
    std::string histogramNameMembrane = generalHistogramName + "Membrane";
    std::string histogramNameCytoplasm = generalHistogramName + "Cytoplasm";
    std::string histogramNameTotalCell = generalHistogramName + "TotalCell";

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


    // Double_t xEdges[XBINS + 1] = {0.0, 0.2, 0.3, 0.6, 0.8, 1.0};
    // TH1* h = new TH1D(
    //   /* name */ "h1",
    //   /* title */ "Hist with variable bin width",
    //   /* number of bins */ NBINS,
    //   /* edge array */ edges
    // )

    int LogBins = 40000;
    double LogWidth[LogBins+1];

    //calculate bins
    for(int i = 0; i <= LogBins; i++)
    {
        LogWidth[i] = std::pow(10.0,std::log10(EMin+0.01)+(std::log10(EMax) - std::log10(EMin+0.01))/double(LogBins)*double(i));
    }

    int hitBins = 50;
    double hitMin = 0.;
    double hitMax = 50.;
    double HitWidth[hitBins+1];

    for(int i=0; i<51; i++)
    {
        HitWidth[i] = i;
    }

    std::string histogramNameNucleus_HitsAlpha =  generalHistogramName + "Nucleus_HitsAlpha";

    hEnergyDepNucleus_HitsAlpha = new TH2D(histogramNameNucleus_HitsAlpha.c_str(), "Energy deposition in nucleus and number of hits by alphas", LogBins, LogWidth, hitBins, HitWidth);


    hEnergyDepsNucleus = new TH1D(histogramNameNucleus.c_str(), "Energy Deposition in Cell Nucleus / Number of Cells", NBins, EMin, EMax);
    hEnergyDepsMembrane = new TH1D(histogramNameMembrane.c_str(), "Energy Depsition in Cell Membrane / Number of Cells", NBins, EMin, EMax);
    hEnergyDepsCytoplasm = new TH1D(histogramNameCytoplasm.c_str(), "Energy Depsition in Cell Cytoplasm / Number of Cells", NBins, EMin, EMax);
    hEnergyDepsTotalCell = new TH1D(histogramNameTotalCell.c_str(), "Energy Depsition in Cell / Number of Cells", NBins, EMin, EMax);

    hEnergyDepsNucleus_FromSolution = new TH1D(histogramNameNucleus_FromSolution.c_str(), "Energy Deposition in Cell Nucleus (Decay originated in Solution) / Number of Cells", NBins, EMin, EMax);
    hEnergyDepsNucleus_FromMembrane = new TH1D(histogramNameNucleus_FromMembrane.c_str(), "Energy Deposition in Cell Nucleus (Decay originated in Membrane) / Number of Cells", NBins, EMin, EMax);
    hEnergyDepsNucleus_FromCytoplasm = new TH1D(histogramNameNucleus_FromCytoplasm.c_str(), "Energy Deposition in Cell Nucleus (Decay originated in Cytoplasm) / Number of Cells", NBins, EMin, EMax);

     hEnergyDepsMembrane_FromSolution = new TH1D(histogramNameMembrane_FromSolution.c_str(), "Energy Deposition in Cell Membrane (Decay originated in Solution) / Number of Cells", NBins, EMin, EMax);
    hEnergyDepsMembrane_FromMembrane = new TH1D(histogramNameMembrane_FromMembrane.c_str(), "Energy Deposition in Cell Membrane (Decay originated in Membrane) / Number of Cells", NBins, EMin, EMax);
    hEnergyDepsMembrane_FromCytoplasm = new TH1D(histogramNameMembrane_FromCytoplasm.c_str(), "Energy Deposition in Cell Membrane (Decay originated in Cytoplasm) / Number of Cells", NBins, EMin, EMax);

     hEnergyDepsCytoplasm_FromSolution = new TH1D(histogramNameCytoplasm_FromSolution.c_str(), "Energy Deposition in Cell Cytoplasm (Decay originated in Solution) / Number of Cells", NBins, EMin, EMax);
    hEnergyDepsCytoplasm_FromMembrane = new TH1D(histogramNameCytoplasm_FromMembrane.c_str(), "Energy Deposition in Cell Cytoplasm (Decay originated in Membrane) / Number of Cells", NBins, EMin, EMax);
    hEnergyDepsCytoplasm_FromCytoplasm = new TH1D(histogramNameCytoplasm_FromCytoplasm.c_str(), "Energy Deposition in Cell Cytoplasm (Decay originated in Cytoplasm) / Number of Cells", NBins, EMin, EMax);

     hEnergyDepsTotalCell_FromSolution = new TH1D(histogramNameTotalCell_FromSolution.c_str(), "Energy Deposition in Cell (Decay originated in Solution) / Number of Cells", NBins, EMin, EMax);
    hEnergyDepsTotalCell_FromMembrane = new TH1D(histogramNameTotalCell_FromMembrane.c_str(), "Energy Deposition in Cell (Decay originated in Membrane) / Number of Cells", NBins, EMin, EMax);
    hEnergyDepsTotalCell_FromCytoplasm = new TH1D(histogramNameTotalCell_FromCytoplasm.c_str(), "Energy Deposition in Cell (Decay originated in Cytoplasm) / Number of Cells", NBins, EMin, EMax);
}

//------------------–----------
void EnergyDepositionHistograms::AddCellHitsToHistograms(CellHit cellHit)
{
    // cellHit.FinalizeCellHit();
    // std::cout << cellHit.GetEnergyDepositionCytoplasm() << std::endl;

    if(cellHit.GetEnergyDepositionMembrane()>0.0){hEnergyDepsMembrane->Fill(cellHit.GetEnergyDepositionMembrane());}
    if(cellHit.GetEnergyDepositionCytoplasm()>0.0){hEnergyDepsCytoplasm->Fill(cellHit.GetEnergyDepositionCytoplasm());}
    if(cellHit.GetEnergyDepositionNucleus()>0.0){hEnergyDepsNucleus->Fill(cellHit.GetEnergyDepositionNucleus());}
    if(cellHit.GetEnergyDepositionTotalCell()>0.0){hEnergyDepsTotalCell->Fill(cellHit.GetEnergyDepositionTotalCell());}

    if(cellHit.GetEnergyDepositionNucleus_FromSolution()>0.0){hEnergyDepsNucleus_FromSolution->Fill(cellHit.GetEnergyDepositionNucleus_FromSolution());}
    if(cellHit.GetEnergyDepositionNucleus_FromMembrane()>0.0){hEnergyDepsNucleus_FromMembrane->Fill(cellHit.GetEnergyDepositionNucleus_FromMembrane());}
    if(cellHit.GetEnergyDepositionNucleus_FromCytoplasm()>0.0){hEnergyDepsNucleus_FromCytoplasm->Fill(cellHit.GetEnergyDepositionNucleus_FromCytoplasm());}

    if(cellHit.GetEnergyDepositionMembrane_FromSolution()>0.0){hEnergyDepsMembrane_FromSolution->Fill(cellHit.GetEnergyDepositionMembrane_FromSolution());}
    if(cellHit.GetEnergyDepositionMembrane_FromMembrane()>0.0){hEnergyDepsMembrane_FromMembrane->Fill(cellHit.GetEnergyDepositionMembrane_FromMembrane());}
    if(cellHit.GetEnergyDepositionMembrane_FromCytoplasm()>0.0){hEnergyDepsMembrane_FromCytoplasm->Fill(cellHit.GetEnergyDepositionMembrane_FromCytoplasm());}

    if(cellHit.GetEnergyDepositionCytoplasm_FromSolution()>0.0){hEnergyDepsCytoplasm_FromSolution->Fill(cellHit.GetEnergyDepositionCytoplasm_FromSolution());}
    if(cellHit.GetEnergyDepositionCytoplasm_FromMembrane()>0.0){hEnergyDepsCytoplasm_FromMembrane->Fill(cellHit.GetEnergyDepositionCytoplasm_FromMembrane());}
    if(cellHit.GetEnergyDepositionCytoplasm_FromCytoplasm()>0.0){hEnergyDepsCytoplasm_FromCytoplasm->Fill(cellHit.GetEnergyDepositionCytoplasm_FromCytoplasm());}

    if(cellHit.GetEnergyDepositionTotalCell_FromSolution()>0.0){hEnergyDepsTotalCell_FromSolution->Fill(cellHit.GetEnergyDepositionTotalCell_FromSolution());}
    if(cellHit.GetEnergyDepositionTotalCell_FromMembrane()>0.0){hEnergyDepsTotalCell_FromMembrane->Fill(cellHit.GetEnergyDepositionTotalCell_FromMembrane());}
    if(cellHit.GetEnergyDepositionTotalCell_FromCytoplasm()>0.0){hEnergyDepsTotalCell_FromCytoplasm->Fill(cellHit.GetEnergyDepositionTotalCell_FromCytoplasm());}

    // if(cellHit.GetEnergyDepositionNucleus()>0.0){std::cout << cellHit.GetNumberHitsAlpha() << std::endl;}
    // std::cout << cellHit.GetNumberHitsAlpha() << ", " << cellHit.GetEnergyDepositionTotalCell() << std::endl;
    if(cellHit.GetEnergyDepositionTotalCell()>0.0){hEnergyDepNucleus_HitsAlpha->Fill(cellHit.GetEnergyDepositionTotalCell(), cellHit.GetNumberHitsAlpha());}

}

//------------------–----------
void EnergyDepositionHistograms::ScaleHistograms(double factor)
{
    hEnergyDepsMembrane->Scale(factor);
    hEnergyDepsCytoplasm->Scale(factor);
    hEnergyDepsNucleus->Scale(factor);
    hEnergyDepsTotalCell->Scale(factor);

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


    hEnergyDepNucleus_HitsAlpha->Scale(factor);
}

//------------------–----------
void EnergyDepositionHistograms::WriteHistogramsToFile()
{
    //------------------–----------
    hEnergyDepsMembrane->Write();
    hEnergyDepsCytoplasm->Write();
    hEnergyDepsNucleus->Write();
    hEnergyDepsTotalCell->Write();

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

    hEnergyDepNucleus_HitsAlpha->Write();


    //------------------–----------
    hEnergyDepsMembrane->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    hEnergyDepsMembrane->GetYaxis()->SetTitle("Fraction of Cells hit");

    hEnergyDepsCytoplasm->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    hEnergyDepsCytoplasm->GetYaxis()->SetTitle("Fraction of Cells hit");

    hEnergyDepsNucleus->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    hEnergyDepsNucleus->GetYaxis()->SetTitle("Fraction of Cells hit");

    hEnergyDepsTotalCell->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    hEnergyDepsTotalCell->GetYaxis()->SetTitle("Fraction of Cells hit");


    hEnergyDepsNucleus_FromSolution->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    hEnergyDepsNucleus_FromSolution->GetYaxis()->SetTitle("Fraction of Cells hit");

    hEnergyDepsNucleus_FromMembrane->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    hEnergyDepsNucleus_FromMembrane->GetYaxis()->SetTitle("Fraction of Cells hit");

    hEnergyDepsNucleus_FromCytoplasm->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    hEnergyDepsNucleus_FromCytoplasm->GetYaxis()->SetTitle("Fraction of Cells hit");


    hEnergyDepsMembrane_FromSolution->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    hEnergyDepsMembrane_FromSolution->GetYaxis()->SetTitle("Fraction of Cells hit");

    hEnergyDepsMembrane_FromMembrane->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    hEnergyDepsMembrane_FromMembrane->GetYaxis()->SetTitle("Fraction of Cells hit");

    hEnergyDepsMembrane_FromCytoplasm->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    hEnergyDepsMembrane_FromCytoplasm->GetYaxis()->SetTitle("Fraction of Cells hit");


    hEnergyDepsCytoplasm_FromSolution->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    hEnergyDepsCytoplasm_FromSolution->GetYaxis()->SetTitle("Fraction of Cells hit");

    hEnergyDepsCytoplasm_FromMembrane->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    hEnergyDepsCytoplasm_FromMembrane->GetYaxis()->SetTitle("Fraction of Cells hit");

    hEnergyDepsCytoplasm_FromCytoplasm->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    hEnergyDepsCytoplasm_FromCytoplasm->GetYaxis()->SetTitle("Fraction of Cells hit");


    hEnergyDepsTotalCell_FromSolution->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    hEnergyDepsTotalCell_FromSolution->GetYaxis()->SetTitle("Fraction of Cells hit");

    hEnergyDepsTotalCell_FromMembrane->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    hEnergyDepsTotalCell_FromMembrane->GetYaxis()->SetTitle("Fraction of Cells hit");

    hEnergyDepsTotalCell_FromCytoplasm->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    hEnergyDepsTotalCell_FromCytoplasm->GetYaxis()->SetTitle("Fraction of Cells hit");

}


//------------------–----------
EnergyDepositionHistograms MakeHistograms(DecayDynamics decayDynamicsInstance, int numberIterations, double volumeRatio, int numberCells)
{
    std::cout << "-----------------------" << std::endl;
    std::cout << "Making Histogram: " << decayDynamicsInstance.GetCellLine() << ", Activity: " << decayDynamicsInstance.GetActivity() << std::endl;
    //------------------–----------
    // Loading decay dynamics
    // double numberDecays212PbInSolution1hTo2h = decayDynamicsInstance.GetNumberDecaysInSolutionFirstHour()*volumeRatio;
    // double numberDecays212PbInMembrane1hTo26h = decayDynamicsInstance.GetNumberDecaysInMembraneTotalTime()*volumeRatio;
    // double numberDecays212PbInCytoplasm1hTo26h = decayDynamicsInstance.GetNumberDecaysInCytoplasmTotalTime()*volumeRatio;
    double numberDecays212PbInSolution1hTo2h = 10.;
    double numberDecays212PbInMembrane1hTo26h = 0.;
    double numberDecays212PbInCytoplasm1hTo26h = 10.;
    double numberCells_Sim = numberCells*volumeRatio;

    // std::cout << "Num: " << numberDecays212PbInSolution1hTo2h << std::endl;
    //------------------–----------
    int NBins = 4000000;
    double EMin = 0.0;
    double EMax = 40.0;

    //------------------–----------
    // Generating empty histograms
    EnergyDepositionHistograms energyDepHistograms = EnergyDepositionHistograms(NBins, EMin, EMax);
    energyDepHistograms.GenerateEmptyHistograms(decayDynamicsInstance);



    std::vector<double> energyDepCytoplasmVec;

    //------------------–----------
    //  Function for filling histograms
    auto FillHistograms = [&](std::string filepathSolutionSim_i, std::string filepathMembraneSim_i, std::string filepathCytoplasmSim_i)
    // auto FillHistograms = [&](TChain* chainSolutionSim, TChain* chainMembraneSim, TChain* chainCytoplasmSim)
    {


        //------------------–----------
        // Opening TTree files and creating TTreeReaders

        // Reader for solution simulation
        std::shared_ptr<TFile> myFileSolutionSim(TFile::Open(filepathSolutionSim_i.c_str(), "READ"));
        auto treeSolutionSim = myFileSolutionSim->Get<TTree>("B4");
        TTreeReader myReaderSolutionSim(treeSolutionSim);


        // // Reader for membrane simulation
        // std::shared_ptr<TFile> myFileMembraneSim(TFile::Open(filepathMembraneSim_i.c_str(), "READ"));
        // auto treeMembraneSim = myFileMembraneSim->Get<TTree>("B4");
        // TTreeReader myReaderMembraneSim(treeMembraneSim);


        // //Reader for cytoplasm simulation
        // std::shared_ptr<TFile> myFileCytoplasmSim(TFile::Open(filepathCytoplasmSim_i.c_str(), "READ"));
        // auto treeCytoplasmSim = myFileCytoplasmSim->Get<TTree>("B4");
        // TTreeReader myReaderCytoplasmSim(treeCytoplasmSim);


        // Reader for membrane simulation
        std::shared_ptr<TFile> myFileMembraneSim(TFile::Open(filepathMembraneSim_i.c_str(), "READ"));
        TTree *treeMembraneSim = 0;

        if(myFileMembraneSim)
        {
            treeMembraneSim = myFileMembraneSim->Get<TTree>("B4");
        }

        TTreeReader myReaderMembraneSim(treeMembraneSim);


        //Reader for cytoplasm simulation
        std::shared_ptr<TFile> myFileCytoplasmSim(TFile::Open(filepathCytoplasmSim_i.c_str(), "READ"));
        TTree *treeCytoplasmSim = 0;

        if(myFileCytoplasmSim)
        {
            treeCytoplasmSim = myFileCytoplasmSim->Get<TTree>("B4");
        }

        TTreeReader myReaderCytoplasmSim(treeCytoplasmSim);


        //------------------–----------
        // Accessing brances of tree

        // Solution
        TTreeReaderArray<double> energyDepsSolutionSim(myReaderSolutionSim, "EnergyDeps");
        TTreeReaderArray<int> volumeTypesSolutionSim(myReaderSolutionSim, "VolumeTypes");
        TTreeReaderArray<int> cellIDsSolutionSim(myReaderSolutionSim, "CellIDs");
        TTreeReaderArray<double> kineticEnergySolutionSim(myReaderSolutionSim, "KineticEnergy");
        TTreeReaderArray<int> particleTypeSolutionSim(myReaderSolutionSim, "ParticleType");
        TTreeReaderArray<double> interactionTimeSolutionSim(myReaderSolutionSim, "InteractionTime");
        TTreeReaderArray<double> firstInteractionTimeSolutionSim(myReaderSolutionSim, "FirstInteractionTime");
        TTreeReaderArray<int> firstInteractionVolumeSolutionSim(myReaderSolutionSim, "FirstInteractionVolume");
        TTreeReaderArray<int> trackIDSolutionSim(myReaderSolutionSim, "TrackID");
        TTreeReaderArray<int> parentIDSolutionSim(myReaderSolutionSim, "ParentID");

        // Membrane
        TTreeReaderArray<double> energyDepsMembraneSim(myReaderMembraneSim, "EnergyDeps");
        TTreeReaderArray<int> volumeTypesMembraneSim(myReaderMembraneSim, "VolumeTypes");
        TTreeReaderArray<int> cellIDsMembraneSim(myReaderMembraneSim, "CellIDs");
        TTreeReaderArray<double> kineticEnergyMembraneSim(myReaderMembraneSim, "KineticEnergy");
        TTreeReaderArray<int> particleTypeMembraneSim(myReaderMembraneSim, "ParticleType");
        TTreeReaderArray<double> interactionTimeMembraneSim(myReaderMembraneSim, "InteractionTime");
        TTreeReaderArray<double> firstInteractionTimeMembraneSim(myReaderMembraneSim, "FirstInteractionTime");
        TTreeReaderArray<int> firstInteractionVolumeMembraneSim(myReaderMembraneSim, "FirstInteractionVolume");
        TTreeReaderArray<int> trackIDMembraneSim(myReaderMembraneSim, "TrackID");
        TTreeReaderArray<int> parentIDMembraneSim(myReaderMembraneSim, "ParentID");

        // Cytoplasm
        TTreeReaderArray<double> energyDepsCytoplasmSim(myReaderCytoplasmSim, "EnergyDeps");
        TTreeReaderArray<int> volumeTypesCytoplasmSim(myReaderCytoplasmSim, "VolumeTypes");
        TTreeReaderArray<int> cellIDsCytoplasmSim(myReaderCytoplasmSim, "CellIDs");
        TTreeReaderArray<double> kineticEnergyCytoplasmSim(myReaderCytoplasmSim, "KineticEnergy");
        TTreeReaderArray<int> particleTypeCytoplasmSim(myReaderCytoplasmSim, "ParticleType");
        TTreeReaderArray<double> interactionTimeCytoplasmSim(myReaderCytoplasmSim, "InteractionTime");
        TTreeReaderArray<double> firstInteractionTimeCytoplasmSim(myReaderCytoplasmSim, "FirstInteractionTime");
        TTreeReaderArray<int> firstInteractionVolumeCytoplasmSim(myReaderCytoplasmSim, "FirstInteractionVolume");
        TTreeReaderArray<int> trackIDCytoplasmSim(myReaderCytoplasmSim, "TrackID");
        TTreeReaderArray<int> parentIDCytoplasmSim(myReaderCytoplasmSim, "ParentID");


        // //---------------------------
        // // Using chains
        // TTreeReader myReaderSolutionSim(chainSolutionSim);
        // TTreeReader myReaderMembraneSim(chainMembraneSim);
        // TTreeReader myReaderCytoplasmSim(chainCytoplasmSim);



        //------------------–----------
        // Vector to store cell hits
        std::vector<CellHit> storedCellHits;

        int decayOriginSolution = 0;
        int decayOriginMembrane = 1;
        int decayOriginCytoplasm = 2;

        //------------------–----------
        //Looping through data for decays occurring in solution from hour 1 to hour 2

        // Counter variable
        int numberDecays212PbInSolution1hTo2h_counter = 0;
        bool whileLoopSolutionSimWasBroken = false;

        while(myReaderSolutionSim.Next())
        {

            //--------------------------------------
            // bool for if the decay occurred in solution from hour 1 to hour 2
            bool firstInteractionInSolution1hTo2h = false;

            // Checking when and where first interaction occured
            if(firstInteractionVolumeSolutionSim[0]==0)
            {

                if(firstInteractionTimeSolutionSim[0]/3600. >= 1. && firstInteractionTimeSolutionSim[0]/3600. <= 2.)
                {

                    // Sett bool true
                    firstInteractionInSolution1hTo2h = true;

                    // Update counter
                    numberDecays212PbInSolution1hTo2h_counter ++;
                }
            }

            //------------------–----------
            // Break loop if number of decays have been reached
            if(numberDecays212PbInSolution1hTo2h_counter >= numberDecays212PbInSolution1hTo2h)
            {
                std::cout << "Solution finished" << std::endl;
                whileLoopSolutionSimWasBroken = true;
                break;
            }

            //--------------------------------------
            // Vector to store <cellID, trackID> for particle hits to a cell for one event
            std::vector<std::tuple<int,int>> particleHitsInfoForEventVec;

            // Sorting through interactions
            if(firstInteractionInSolution1hTo2h)
            {
                //--------------------------------------
                // looping over all steps for one event/decay
                for(int i=0; i<energyDepsSolutionSim.GetSize(); i++)
                {

                    // Only add if actual energy deposition
                    if(energyDepsSolutionSim[i]>0.)
                    {

                        // Only add if interaction happened between hour 1 and hour 2
                        if(interactionTimeSolutionSim[i]/3600.0 >= 1. && interactionTimeSolutionSim[i]/3600.0 <= 2.)
                        {

                            // boolean: true if cell has not already been hit before
                            // false if already hit
                            bool cellAlreadyHit = false;


                            //--------------------------------------
                            // Loop over all previous cell hits stored
                            for(int ii=0; ii<storedCellHits.size(); ii++)
                            {

                                //--------------------------------------
                                // Add energy deposition to cell that has already been hit before
                                if(cellIDsSolutionSim[i]==storedCellHits[ii].GetCellID())
                                {
                                    // Add energy deposition
                                    // std::cout << "Interaction volume " << volumeTypesSolutionSim[i] << std::endl;
                                    storedCellHits[ii].AddEnergyDeposition(energyDepsSolutionSim[i],volumeTypesSolutionSim[i],decayOriginSolution);

                                    // Update boolean
                                    cellAlreadyHit = true;


                                    //---------------------------------------------
                                    // boolean: true if the particle track has already been stored for a specific cell ID
                                    bool trackAlreadyRegistered = false;

                                    // Looping over all previously stored track IDs and corresponding cell IDs for this event
                                    for(int iii=0; iii<particleHitsInfoForEventVec.size(); iii++)
                                    {

                                        // If track ID is already stored for this cell ID, do not update particle hit counts for this cell ID
                                        if(cellIDsSolutionSim[i] == std::get<0>(particleHitsInfoForEventVec[iii]))
                                        {
                                            if(trackIDSolutionSim[i]==std::get<1>(particleHitsInfoForEventVec[iii]))
                                            {
                                                trackAlreadyRegistered = true;
                                            }
                                        }
                                    }
                                    // If track ID has not already been registered or not been registered for this cell ID
                                    if(!trackAlreadyRegistered)
                                    {

                                        // If alpha particle type
                                        if(particleTypeSolutionSim[i]==1000020040)
                                        {
                                            storedCellHits[ii].HitByAlphaParticle();
                                        }

                                        // Save the hit and cellID
                                        std::tuple<int,int> particleHit = std::make_tuple(cellIDsSolutionSim[i],trackIDSolutionSim[i]);
                                        particleHitsInfoForEventVec.push_back(particleHit);

                                    }
                                }
                            }

                            //--------------------------------------
                            // Register a new cell hit, and add energy deposition
                            if(!cellAlreadyHit)
                            {
                                // std::cout << "Interaction volume " << volumeTypesSolutionSim[i] << std::endl;
                                CellHit aNewCellHit = CellHit(cellIDsSolutionSim[i]);
                                aNewCellHit.AddEnergyDeposition(energyDepsSolutionSim[i], volumeTypesSolutionSim[i],decayOriginSolution);
                                storedCellHits.push_back(aNewCellHit);


                                // If alpha particle
                                if(particleTypeSolutionSim[i]==1000020040)
                                {
                                    aNewCellHit.HitByAlphaParticle();
                                }

                                std::tuple<int,int> particleHit = std::make_tuple(cellIDsSolutionSim[i],trackIDSolutionSim[i]);
                                particleHitsInfoForEventVec.push_back(particleHit);

                            }
                        }
                    }
                }
            }
        }

        //--------------------------------------
        // Check if enough decays were processed
        if(!whileLoopSolutionSimWasBroken)
        {
            std::cout << "Not enough events in solution simulation file! Need " << numberDecays212PbInSolution1hTo2h << " number of decays. Only reached " << numberDecays212PbInSolution1hTo2h_counter << " number of decays at entry number " << myReaderSolutionSim.GetCurrentEntry() << std::endl;
        }


        //------------------–----------
        // Looping through data for decays occuring in Membrane from hour 1 to hour 26

        // Counter variable
        int numberDecays212PbInMembrane1hTo26h_counter = 0;
        bool whileLoopMembraneSimWasBroken = false;

        if(myFileMembraneSim)
        {
            while(myReaderMembraneSim.Next())
            {
                //--------------------------------------
                // Chacking if first interaction took place between 1 and 26 hours
                bool firstInteractionInMembrane1hTo26h = false;

                if(firstInteractionTimeMembraneSim[0]/3600. >= 1. && firstInteractionTimeMembraneSim[0]/3600. <= 26.)
                {
                    firstInteractionInMembrane1hTo26h = true;
                    numberDecays212PbInMembrane1hTo26h_counter ++;
                }

                //--------------------------------------
                // Break loop if number of decays have been reached
                if(numberDecays212PbInMembrane1hTo26h_counter >= numberDecays212PbInMembrane1hTo26h)
                {
                    whileLoopMembraneSimWasBroken = true;
                    std::cout << "Membrane finished" << std::endl;
                    break;
                }

                //--------------------------------------
                // Vector to store <cellID, trackID> for particle hits to a cell for one event
                std::vector<std::tuple<int,int>> particleHitsInfoForEventVec;

                if(firstInteractionInMembrane1hTo26h)
                {
                    //--------------------------------------
                    // looping over all steps for one event/decay
                    for(int i=0; i<energyDepsMembraneSim.GetSize(); i++)
                    {

                        // Only add if actual energy deposition
                        if(energyDepsMembraneSim[i]>0.)
                        {

                            // Only add if interaction happened between hour 1 and hour 26
                            if(interactionTimeMembraneSim[i]/3600.0 >= 1. && interactionTimeMembraneSim[i]/3600.0 <= 26.)
                            {

                                // boolean: true if cell has already been hit before
                                bool cellAlreadyHit = false;

                                //--------------------------------------
                                // Loop over all previous cell hits stored
                                for(int ii=0; ii<storedCellHits.size(); ii++)
                                {

                                    // Add energy deposition to cell that has already been hit before
                                    if(cellIDsMembraneSim[i]==storedCellHits[ii].GetCellID())
                                    {
                                        // std::cout << "Interaction volume " << volumeTypesSolutionSim[i] << std::endl;

                                        storedCellHits[ii].AddEnergyDeposition(energyDepsMembraneSim[i],volumeTypesMembraneSim[i],decayOriginMembrane);

                                        // Update boolean
                                        cellAlreadyHit = true;

                                        //--------------------------------------
                                        // boolean: true if the particle track has already been stored for a specific cell ID
                                        bool trackAlreadyRegistered = false;

                                        // Looping over all previously stored track IDs and corresponding cell IDs for this event
                                        for(int iii=0; iii<particleHitsInfoForEventVec.size(); iii++)
                                        {

                                            // If track ID is already stored for this cell ID, do not update particle hit counts for this cell ID
                                            if(cellIDsMembraneSim[i] == std::get<0>(particleHitsInfoForEventVec[iii]))
                                            {
                                                if(trackIDMembraneSim[i]==std::get<1>(particleHitsInfoForEventVec[iii]))
                                                {
                                                    trackAlreadyRegistered = true;
                                                }
                                            }
                                        }
                                        // If track ID has not already been registered or not been registered for this cell ID
                                        if(!trackAlreadyRegistered)
                                        {

                                            // If alpha particle type
                                            if(particleTypeMembraneSim[i]==1000020040)
                                            {
                                                storedCellHits[ii].HitByAlphaParticle();
                                            }

                                            // Save the hit and cellID
                                            std::tuple<int,int> particleHit = std::make_tuple(cellIDsMembraneSim[i],trackIDMembraneSim[i]);
                                            particleHitsInfoForEventVec.push_back(particleHit);

                                        }
                                    }
                                }

                                // Register a new cell hit, and add energy deposition
                                if(!cellAlreadyHit)
                                {
                                    // std::cout << "Interaction volume " << volumeTypesSolutionSim[i] << std::endl;

                                    CellHit aNewCellHit = CellHit(cellIDsMembraneSim[i]);
                                    aNewCellHit.AddEnergyDeposition(energyDepsMembraneSim[i], volumeTypesMembraneSim[i], decayOriginMembrane);
                                    storedCellHits.push_back(aNewCellHit);

                                    // If alpha particle
                                    if(particleTypeMembraneSim[i]==1000020040)
                                    {
                                        aNewCellHit.HitByAlphaParticle();
                                    }

                                    std::tuple<int,int> particleHit = std::make_tuple(cellIDsMembraneSim[i],trackIDMembraneSim[i]);
                                    particleHitsInfoForEventVec.push_back(particleHit);
                                }
                            }
                        }
                    }
                }
            }

            //--------------------------------------
            // Checking if enough decays were processed
            if(!whileLoopMembraneSimWasBroken)
            {
                std::cout << "Not enough events in membrane simulation file! Need " << numberDecays212PbInMembrane1hTo26h << " number of decays. Only reached " << numberDecays212PbInMembrane1hTo26h_counter << " number of decays at entry number " << myReaderMembraneSim.GetCurrentEntry() << std::endl;
            }
        }



        // ------------------–----------
        // Looping through data for decays occuring in Cytoplasm between hour 1 and hour 26

        // Counter variable
        int numberDecays212PbInCytoplasm1hTo26h_counter = 0;
        bool whileLoopCytoplasmSimWasBroken = false;

        if(myFileCytoplasmSim)
        {
            while(myReaderCytoplasmSim.Next())
            {

                //--------------------------------------
                // Chacking if first interaction took place between 1 and 26 hours
                bool firstInteractionInCytoplasm1hTo26h = false;

                if(firstInteractionVolumeCytoplasmSim[0] == 2)
                {
                    if(firstInteractionTimeCytoplasmSim[0]/3600. >= 1. && firstInteractionTimeCytoplasmSim[0]/3600. <= 26.)
                    {
                        firstInteractionInCytoplasm1hTo26h = true;
                        numberDecays212PbInCytoplasm1hTo26h_counter ++;
                    }
                }

                //--------------------------------------
                // Break loop if number of decays have been reached
                if(numberDecays212PbInCytoplasm1hTo26h_counter >= numberDecays212PbInCytoplasm1hTo26h)
                {
                    std::cout << "Cytoplasm finished" << std::endl;
                    whileLoopCytoplasmSimWasBroken = true;
                    break;
                }

                //--------------------------------------
                // Vector to store <cellID, trackID> for particle hits to a cell for one event
                std::vector<std::tuple<int,int>> particleHitsInfoForEventVec;

                if(firstInteractionInCytoplasm1hTo26h)
                {
                    // looping over all steps for one event/decay
                    for(int i=0; i<energyDepsCytoplasmSim.GetSize(); i++)
                    {

                        // Only add if actual energy deposition
                        if(energyDepsCytoplasmSim[i]>0.)
                        {

                            // Only add if interaction happened within between hour 1 and hour 26
                            if(interactionTimeCytoplasmSim[i]/3600.0 >= 1. && interactionTimeCytoplasmSim[i]/3600.0 <= 26.)
                            {

                                //--------------------------------------
                                // boolean: true if cell has already been hit before in this event
                                bool cellAlreadyHit = false;

                                // Loop over all previous cell hits stored
                                for(int ii=0; ii<storedCellHits.size(); ii++)
                                {

                                    // Add energy deposition to cell that has already been hit before
                                    if(cellIDsCytoplasmSim[i]==storedCellHits[ii].GetCellID())
                                    {
                                        // std::cout << "Interaction volume " << volumeTypesSolutionSim[i] << std::endl;

                                        storedCellHits[ii].AddEnergyDeposition(energyDepsCytoplasmSim[i],volumeTypesCytoplasmSim[i],decayOriginCytoplasm);

                                        // Update boolean
                                        cellAlreadyHit = true;

                                        // boolean: true if the particle track has already been stored for a specific cell ID
                                        bool trackAlreadyRegistered = false;

                                        // Looping over all previously stored track IDs and corresponding cell IDs for this event
                                        for(int iii=0; iii<particleHitsInfoForEventVec.size(); iii++)
                                        {

                                            // If track ID is already stored for this cell ID, do not update particle hit counts for this cell ID
                                            if(cellIDsCytoplasmSim[i] == std::get<0>(particleHitsInfoForEventVec[iii]))
                                            {
                                                if(trackIDCytoplasmSim[i]==std::get<1>(particleHitsInfoForEventVec[iii]))
                                                {
                                                    trackAlreadyRegistered = true;
                                                }
                                            }
                                        }
                                        // If track ID has not already been registered or not been registered for this cell ID
                                        if(!trackAlreadyRegistered)
                                        {

                                            // If alapha particle type
                                            if(particleTypeCytoplasmSim[i]==1000020040)
                                            {
                                                storedCellHits[ii].HitByAlphaParticle();
                                                // std::cout << "New track : "<< cellIDsCytoplasmSim[i] << " , " << trackIDCytoplasmSim[i] << " , " << particleTypeCytoplasmSim[i] << std::endl;
                                            }

                                            // Save the hit and cellID
                                            std::tuple<int,int> particleHit = std::make_tuple(cellIDsCytoplasmSim[i],trackIDCytoplasmSim[i]);
                                            particleHitsInfoForEventVec.push_back(particleHit);

                                        }
                                    }
                                }
                                // Register a new cell hit, and add energy deposition
                                if(!cellAlreadyHit)
                                {
                                    // std::cout << "Interaction volume " << volumeTypesSolutionSim[i] << std::endl;

                                    CellHit aNewCellHit = CellHit(cellIDsCytoplasmSim[i]);
                                    aNewCellHit.AddEnergyDeposition(energyDepsCytoplasmSim[i], volumeTypesCytoplasmSim[i],decayOriginCytoplasm);
                                    storedCellHits.push_back(aNewCellHit);

                                    // If alpha particle
                                    if(particleTypeCytoplasmSim[i]==1000020040)
                                    {
                                        aNewCellHit.HitByAlphaParticle();
                                    }

                                    std::tuple<int,int> particleHit = std::make_tuple(cellIDsCytoplasmSim[i],trackIDCytoplasmSim[i]);
                                    particleHitsInfoForEventVec.push_back(particleHit);
                                }
                            }
                        }
                    }

                }
            }

            //--------------------------------------
            // Checking if enough decays where processed
            if(!whileLoopCytoplasmSimWasBroken)
            {
                std::cout << "Not enough events in cytoplasm simulation file! Need " << numberDecays212PbInCytoplasm1hTo26h << " number of decays. Only reached " << numberDecays212PbInCytoplasm1hTo26h_counter << " number of decays at entry number " << myReaderCytoplasmSim.GetCurrentEntry() << std::endl;
            }
        }

        //------------------–----------
        // Looping over all stored cell energy depositions
        for(int i=0; i<storedCellHits.size(); i++)
        {
            // Adding cell energy depositions to histograms
            storedCellHits[i].FinalizeCellHit();
            energyDepCytoplasmVec.push_back(storedCellHits[i].GetEnergyDepositionCytoplasm());
            energyDepHistograms.AddCellHitsToHistograms(storedCellHits[i]);
        }
    };


    double maxE = *std::max_element(energyDepCytoplasmVec.begin(), energyDepCytoplasmVec.end());
    std::cout << "Max energy cytoplasm " << std::endl;
    // //----------------------------
    // // Making chain from files
    // TChain* chSolution = new TChain("B4");
    // TChain* chMembrane = new TChain("B4");
    // TChain* chCytoplasm = new TChain("B4");

    // std::string fileSolution;
    // std::string fileMembrane;
    // std::string fileCytoplasm;

    // for(int i=0; i<10; i++)
    // {
    //     fileSolution = "/Volumes/SamsungT7/Output_Solution_Thread_" + std::to_string(i) + ".root";
    //     fileMembrane = "/Volumes/SamsungT7/Output_Membrane_Thread_" + std::to_string(i) + ".root";
    //     fileCytoplasm = "/Volumes/SamsungT7/Output_Cytoplasm_Thread_" + std::to_string(i) + ".root";

    //     chSolution->Add(fileSolution.c_str());
    //     chMembrane->Add(fileMembrane.c_str());
    //     chCytoplasm->Add(fileCytoplasm.c_str());
    // }


    // std::string filepathSimulationOutput = "../../GEANT4Simulations/OutputFromSaga/";
    // std::string filepathSimulationOutput = "../../GEANT4Simulations/CellDamageSimulation-build/";
    std::string filepathSimulationOutput = "/Volumes/SamsungT7/";

    std::string filepathSolutionIteration_i;
    std::string filepathMembraneIteration_i;
    std::string filepathCytoplasmIteration_i;

    //------------------–----------
    // Filling histgrams
    for(int i=0; i<numberIterations; i++)
    {
        std::cout << "Iteration number : " << i << std::endl;
        filepathSolutionIteration_i = filepathSimulationOutput + "Output_Pb212_" + decayDynamicsInstance.GetCellLine() + "_Activity" + std::to_string(decayDynamicsInstance.GetActivity()) + "kBq_Solution_Thread_" + std::to_string(i) + ".root";
        filepathMembraneIteration_i = filepathSimulationOutput + "Output_Pb212_" + decayDynamicsInstance.GetCellLine() + "_Activity" + std::to_string(decayDynamicsInstance.GetActivity()) + "kBq_Membrane_Thread_" + std::to_string(i) + ".root";
        filepathCytoplasmIteration_i = filepathSimulationOutput + "Output_Pb212_" + decayDynamicsInstance.GetCellLine() + "_Activity" + std::to_string(decayDynamicsInstance.GetActivity()) + "kBq_Cytoplasm_Thread_" + std::to_string(i) + ".root";


        FillHistograms(filepathSolutionIteration_i, filepathMembraneIteration_i, filepathCytoplasmIteration_i);
        // // -----------------------------
        // // Filling using chains
        // FillHistograms(chSolution, chMembrane, chCytoplasm);
        // std::cout << " Activity " << decayDynamicsInstance.GetActivity() <<", finished iteration number: " << i+1 << std::endl;
    }


    //------------------–----------
    // Scaling histograms
    double scalingFactor = (numberCells_Sim)*((double) numberIterations);
    energyDepHistograms.ScaleHistograms(1./scalingFactor);



    //------------------–----------
    return energyDepHistograms;

}

void mainAnalysisCode()
{
    // Calculating volume ratio
    double VolumeSample = 0.2*1000; // mm^3
    double volumeCellTube = TMath::Pi()*std::pow(0.5,2.0)*1.0; // mm^3
    double volumeRatio = volumeCellTube/VolumeSample;
    int numberCells = 500000;

    int numberIterations = 1;

    // std::cout << volumeRatio << std::endl;


    //------------------–----------
    // Defining decay dynamics


    // C4-2 Cells

    DecayDynamics decays_A5kBq_C4_2 = DecayDynamics(5,1.14,1.16,"C4_2");
    DecayDynamics decays_A10kBq_C4_2 = DecayDynamics(10,1.97,1.98,"C4_2");
    DecayDynamics decays_A25kBq_C4_2 = DecayDynamics(25,4.78,5.91,"C4_2");
    DecayDynamics decays_A50kBq_C4_2 = DecayDynamics(50,8.94,11.07,"C4_2");
    DecayDynamics decays_A75kBq_C4_2 = DecayDynamics(75,10.79,13.12,"C4_2");
    DecayDynamics decays_A100kBq_C4_2 = DecayDynamics(100,13.16,22.72,"C4_2");
    DecayDynamics decays_A150kBq_C4_2 = DecayDynamics(150,16.40,23.56,"C4_2");


    // PC3 PIP Cells

    DecayDynamics decays_A10kBq_PC3_PIP = DecayDynamics(10, 47., 2., "PC3_PIP");
    DecayDynamics decays_A25kBq_PC3_PIP = DecayDynamics(25, 119., 229., "PC3_PIP");
    DecayDynamics decays_A50kBq_PC3_PIP = DecayDynamics(50, 229., 11., "PC3_PIP");
    DecayDynamics decays_A75kBq_PC3_PIP = DecayDynamics(75, 335., 17., "PC3_PIP");
    DecayDynamics decays_A100kBq_PC3_PIP = DecayDynamics(100, 448., 22., "PC3_PIP");
    DecayDynamics decays_A150kBq_PC3_PIP = DecayDynamics(150, 565., 28., "PC3_PIP");


    // PC3 Flu Cells
    DecayDynamics decays_A10kBq_PC3_Flu = DecayDynamics(10, 0., 0., "PC3_Flu");
    DecayDynamics decays_A25kBq_PC3_Flu = DecayDynamics(25, 0., 0., "PC3_Flu");
    DecayDynamics decays_A50kBq_PC3_Flu = DecayDynamics(50, 0., 0., "PC3_Flu");
    DecayDynamics decays_A75kBq_PC3_Flu = DecayDynamics(75, 0., 0., "PC3_Flu");
    DecayDynamics decays_A100kBq_PC3_Flu = DecayDynamics(100, 0., 0., "PC3_Flu");
    DecayDynamics decays_A150kBq_PC3_Flu = DecayDynamics(150, 0., 0., "PC3_Flu");
    // DecayDynamics decays_A1000kBq_PC3_Flu = DecayDynamics(1000,0.,0., "PC3_Flu");
    // DecayDynamics decays_A300kBq_PC3_Flu = DecayDynamics(300,0.,0., "PC3_Flu");
    // DecayDynamics decays_A500kBq_PC3_Flu = DecayDynamics(500,0.,0., "PC3_Flu");
    // DecayDynamics decays_A200kBq_PC3_Flu = DecayDynamics(200,0.,0., "PC3_Flu");
    // DecayDynamics decays_A250kBq_PC3_Flu = DecayDynamics(250,0.,0., "PC3_Flu");




    //------------------–----------
    // Loading decay dynamics calculations

    std::string mathematicaOutput = "../../Mathematica/Output";

    decays_A5kBq_C4_2.LoadDataFromMathematicaCalculations(mathematicaOutput.c_str());
    decays_A10kBq_C4_2.LoadDataFromMathematicaCalculations(mathematicaOutput.c_str());
    decays_A25kBq_C4_2.LoadDataFromMathematicaCalculations(mathematicaOutput.c_str());
    decays_A50kBq_C4_2.LoadDataFromMathematicaCalculations(mathematicaOutput.c_str());
    decays_A75kBq_C4_2.LoadDataFromMathematicaCalculations(mathematicaOutput.c_str());
    decays_A100kBq_C4_2.LoadDataFromMathematicaCalculations(mathematicaOutput.c_str());
    decays_A150kBq_C4_2.LoadDataFromMathematicaCalculations(mathematicaOutput.c_str());

    decays_A10kBq_PC3_PIP.LoadDataFromMathematicaCalculations(mathematicaOutput.c_str());
    decays_A25kBq_PC3_PIP.LoadDataFromMathematicaCalculations(mathematicaOutput.c_str());
    decays_A50kBq_PC3_PIP.LoadDataFromMathematicaCalculations(mathematicaOutput.c_str());
    decays_A75kBq_PC3_PIP.LoadDataFromMathematicaCalculations(mathematicaOutput.c_str());
    decays_A100kBq_PC3_PIP.LoadDataFromMathematicaCalculations(mathematicaOutput.c_str());
    decays_A150kBq_PC3_PIP.LoadDataFromMathematicaCalculations(mathematicaOutput.c_str());

    decays_A10kBq_PC3_Flu.LoadDataFromMathematicaCalculations(mathematicaOutput.c_str());
    decays_A25kBq_PC3_Flu.LoadDataFromMathematicaCalculations(mathematicaOutput.c_str());
    decays_A50kBq_PC3_Flu.LoadDataFromMathematicaCalculations(mathematicaOutput.c_str());
    decays_A75kBq_PC3_Flu.LoadDataFromMathematicaCalculations(mathematicaOutput.c_str());
    decays_A100kBq_PC3_Flu.LoadDataFromMathematicaCalculations(mathematicaOutput.c_str());
    decays_A150kBq_PC3_Flu.LoadDataFromMathematicaCalculations(mathematicaOutput.c_str());
    // decays_A1000kBq_PC3_Flu.LoadDataFromMathematicaCalculations(mathematicaOutput.c_str());
    // decays_A300kBq_PC3_Flu.LoadDataFromMathematicaCalculations(mathematicaOutput.c_str());
    // decays_A500kBq_PC3_Flu.LoadDataFromMathematicaCalculations(mathematicaOutput.c_str());
    // decays_A200kBq_PC3_Flu.LoadDataFromMathematicaCalculations(mathematicaOutput.c_str());
    // decays_A250kBq_PC3_Flu.LoadDataFromMathematicaCalculations(mathematicaOutput.c_str());



    //-------------------------------------
    // Creating energy deposition histograms

    // EnergyDepositionHistograms Hist_A150kBq_PC3_Flu = MakeHistograms(decays_A150kBq_PC3_Flu, numberIterations, volumeRatio, numberCells);
    // auto output = new TFile("Output_PC3_Flu_150kBq.root", "RECREATE");
    // Hist_A150kBq_PC3_Flu.WriteHistogramsToFile();
    // output->Write();
    // output->Close();

    // EnergyDepositionHistograms Hist_A100kBq_PC3_PIP = MakeHistograms(decays_A100kBq_PC3_PIP, numberIterations, volumeRatio, numberCells);
    // auto output = new TFile("Output_PC3_PIP_100kBq.root", "RECREATE");
    // Hist_A100kBq_PC3_PIP.WriteHistogramsToFile();
    // output->Write();
    // output->Close();

    EnergyDepositionHistograms Hist_A5kBq_C4_2 = MakeHistograms(decays_A5kBq_C4_2, numberIterations, volumeRatio, numberCells);
    auto output = new TFile("Output_C4_2_5kBq.root", "RECREATE");
    Hist_A5kBq_C4_2.WriteHistogramsToFile();
    output->Write();
    output->Close();
}

