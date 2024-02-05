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

        double GetVolumeRatio(){return volumeRatio;};
        double GetNumberCells(){return numberCells;};

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

        // Volume ratio between simulation and actual sample volume
        double volumeRatio;

        // Number of cells in simualted volume
        double numberCells;


};


DecayDynamics::DecayDynamics(int activitySample_in, double U0PerCell_in, double U0InternalizedPerCell_in, std::string cellLine_in)
{
    activitySample = activitySample_in;
    U0PerCell = U0PerCell_in;
    U0InternalizedPerCell = U0InternalizedPerCell_in;
    cellLine = cellLine_in;

    double VolumeSample = 0.2*1000; // mm^3
    double volumeCellTube = TMath::Pi()*std::pow(0.5,2.0)*1.0; // mm^3
    volumeRatio = volumeCellTube/VolumeSample;

    numberCells = 500000.*volumeRatio;
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


//------------------–----------
CellHit::CellHit(int cellID_in)
{
    cellID = cellID_in;
    nummberHitsAlphas_Nucleus = 0;
    numberHitsAlphas_Membrane = 0;
    numberHitsAlphas_Cytoplasm = 0;
    numberHitsAlphas_TotalCell = 0;

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

    energyDepNucleus_FromAlpha = 0.;
    energyDepMembrane_FromAlpha = 0.;
    energyDepCytoplasm_FromAlpha = 0.;
    energyDepTotalCell_FromAlpha = 0.;

    numberHitsAlphas = 0;
    numberHitsBetas = 0;

    double densityWater = 1000. ; // kg/m^3
    double radiusCytoplasm = 9.0e-6; // m
    double radiusNucleus = 2.5e-6; // m
    double radiusCell = radiusCytoplasm + 4.e-9; // m

    massNucleus = (4./3.)*TMath::Pi()*std::pow(radiusNucleus,3.)*densityWater; // kg
    massCytoplasm = (4./3.)*TMath::Pi()*std::pow(radiusCytoplasm,3.)*densityWater - massNucleus; // kg
    massCell = (4./3.)*TMath::Pi()*std::pow(radiusCell,3.)*densityWater; // kg
    massMembrane = massCell - (4./3.)*TMath::Pi()*std::pow(radiusCytoplasm,3.)*densityWater; // kg
}


//------------------–----------
void CellHit::AddEnergyDeposition(double energyDep_in, int volumeTypeInteraction_in, int volumeTypeOriginDecay_in, int particleType_in)
{
    std::tuple<int,double,int,int> hitTuple;
    hitTuple = std::make_tuple(volumeTypeInteraction_in, energyDep_in, volumeTypeOriginDecay_in, particleType_in);
    energyDepsVec.push_back(hitTuple);
}


//------------------–----------
void CellHit::HitByAlphaParticle(int volumeTypeHit, bool firstTimeCountingAlpha)
{
    if(volumeTypeHit==1){numberHitsAlphas_Membrane++;}
    if(volumeTypeHit==2){numberHitsAlphas_Cytoplasm++;}
    if(volumeTypeHit==3){nummberHitsAlphas_Nucleus++;}

    if(firstTimeCountingAlpha){numberHitsAlphas_TotalCell++;}
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
    int particleType;

    for(int i=0; i<energyDepsVec.size(); i++)
    {
        interactionVolume = std::get<0>(energyDepsVec[i]);
        energyDepInteraction = std::get<1>(energyDepsVec[i]);
        originVolume = std::get<2>(energyDepsVec[i]);
        particleType = std::get<3>(energyDepsVec[i]);


        if(interactionVolume==1)
        {
            // std::cout << "In membrane: " << energyDepInteraction << std::endl;
            energyDepMembrane += energyDepInteraction;
            if(originVolume == 0){energyDepMembrane_FromSolution += energyDepInteraction;}
            if(originVolume == 1){energyDepMembrane_FromMembrane += energyDepInteraction;}
            if(originVolume == 2){energyDepMembrane_FromCytoplasm += energyDepInteraction;}

            if(particleType==1000020040)
            {
                energyDepMembrane_FromAlpha += energyDepInteraction;
            }
        }
        if(interactionVolume==2)
        {
            // std::cout << "In cytoplasm : " << energyDepInteraction << std::endl;
            energyDepCytoplasm += energyDepInteraction;
            if(originVolume == 0){energyDepCytoplasm_FromSolution += energyDepInteraction;}
            if(originVolume == 1){energyDepCytoplasm_FromMembrane += energyDepInteraction;}
            if(originVolume == 2){energyDepCytoplasm_FromCytoplasm += energyDepInteraction;}

            if(particleType==1000020040)
            {
                energyDepCytoplasm_FromAlpha += energyDepInteraction;
            }
        }
        if(interactionVolume==3)
        {
            // std::cout << "In nucleus : " << energyDepInteraction << std::endl;
            energyDepNucleus += energyDepInteraction;
            if(originVolume == 0){energyDepNucleus_FromSolution += energyDepInteraction;}
            if(originVolume == 1){energyDepNucleus_FromMembrane += energyDepInteraction;}
            if(originVolume == 2){energyDepNucleus_FromCytoplasm += energyDepInteraction;}

            if(particleType==1000020040)
            {
                energyDepNucleus_FromAlpha += energyDepInteraction;
            }
        }

    }

    energyDepTotalCell = energyDepMembrane + energyDepCytoplasm + energyDepNucleus;

    energyDepTotalCell_FromSolution = energyDepMembrane_FromSolution + energyDepCytoplasm_FromSolution + energyDepNucleus_FromSolution;
    energyDepTotalCell_FromMembrane = energyDepMembrane_FromMembrane + energyDepCytoplasm_FromMembrane + energyDepNucleus_FromMembrane;
    energyDepTotalCell_FromCytoplasm = energyDepMembrane_FromCytoplasm + energyDepCytoplasm_FromCytoplasm + energyDepNucleus_FromCytoplasm;

    energyDepTotalCell_FromAlpha = energyDepMembrane_FromAlpha + energyDepNucleus_FromAlpha + energyDepCytoplasm_FromAlpha;

}


//------------------–----------
class EnergyDepositionHistograms
{
    public:
        EnergyDepositionHistograms(double EMax_in, int iteration);
        // ~EnergyDepositionHistograms();

        void GenerateEmptyHistograms(DecayDynamics decayDynamicsInstance);
        void AddCellHitsToHistograms(CellHit cellHit);
        void ScaleHistograms(double factor);
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

        // General hit histogram
        TH1D *hFractionHitsAlpha_TotalCell;
        TH1D *hFractionHitsAlpha_Nucleus;
        TH1D *hFractionHitsAlpha_Membrane;
        TH1D *hFractionHitsAlpha_Cytoplasm;

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

        // Histograms for number of hits by alpha particles
        TH2D* hDoseTotalCell_HitsAlpha;
        TH2D* hDoseNucleus_HitsAlpha;
        TH2D* hDoseMembrane_HitsAlpha;
        TH2D* hDoseCytoplasm_HitsAlpha;
};



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

    hDoseNucleus_FromSolution = new TH1D(histogramNameNucleus_FromSolution.c_str(), "Dose Delivered in Cell Nucleus (Decay originated in Solution) ", NBins_mGyBinning, MinDose_mGyBinning, MaxDose_uGyBinning);
    hDoseNucleus_FromMembrane = new TH1D(histogramNameNucleus_FromMembrane.c_str(), "Dose Delivered in Cell Nucleus (Decay originated in Membrane) ", NBins_mGyBinning, MinDose_mGyBinning, MaxDose_uGyBinning);
    hDoseNucleus_FromCytoplasm = new TH1D(histogramNameNucleus_FromCytoplasm.c_str(), "Dose Delivered in Cell Nucleus (Decay originated in Cytoplasm) ", NBins_mGyBinning, MinDose_mGyBinning, MaxDose_uGyBinning);

    hDoseMembrane_FromSolution = new TH1D(histogramNameMembrane_FromSolution.c_str(), "Dose Delivered in Cell Membrane (Decay originated in Solution) ", NBins_mGyBinning, MinDose_mGyBinning, MaxDose_uGyBinning);
    hDoseMembrane_FromMembrane = new TH1D(histogramNameMembrane_FromMembrane.c_str(), "Dose Delivered in Cell Membrane (Decay originated in Membrane) ", NBins_mGyBinning, MinDose_mGyBinning, MaxDose_uGyBinning);
    hDoseMembrane_FromCytoplasm = new TH1D(histogramNameMembrane_FromCytoplasm.c_str(), "Dose Delivered in Cell Membrane (Decay originated in Cytoplasm) ", NBins_mGyBinning, MinDose_mGyBinning, MaxDose_uGyBinning);

    hDoseCytoplasm_FromSolution = new TH1D(histogramNameCytoplasm_FromSolution.c_str(), "Dose Delivered in Cell Cytoplasm (Decay originated in Solution) ", NBins_mGyBinning, MinDose_mGyBinning, MaxDose_uGyBinning);
    hDoseCytoplasm_FromMembrane = new TH1D(histogramNameCytoplasm_FromMembrane.c_str(), "Dose Delivered in Cell Cytoplasm (Decay originated in Membrane) ", NBins_mGyBinning, MinDose_mGyBinning, MaxDose_uGyBinning);
    hDoseCytoplasm_FromCytoplasm = new TH1D(histogramNameCytoplasm_FromCytoplasm.c_str(), "Dose Delivered in Cell Cytoplasm (Decay originated in Cytoplasm) ", NBins_mGyBinning, MinDose_mGyBinning, MaxDose_uGyBinning);

    hDoseTotalCell_FromSolution = new TH1D(histogramNameTotalCell_FromSolution.c_str(), "Dose Delivered in Cell (Decay originated in Solution) ", NBins_mGyBinning, MinDose_mGyBinning, MaxDose_uGyBinning);
    hDoseTotalCell_FromMembrane = new TH1D(histogramNameTotalCell_FromMembrane.c_str(), "Dose Delivered in Cell (Decay originated in Membrane) ", NBins_mGyBinning, MinDose_mGyBinning, MaxDose_uGyBinning);
    hDoseTotalCell_FromCytoplasm = new TH1D(histogramNameTotalCell_FromCytoplasm.c_str(), "Dose Delivered in Cell (Decay originated in Cytoplasm) ", NBins_mGyBinning, MinDose_mGyBinning, MaxDose_uGyBinning);

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
    hDoseTotalCell_HitsAlpha = new TH2D(histogramNameTotalCell_HitsAlpha.c_str(), "Dose Delivered in Cell per N Number of hits by alphas", NBins_mGyBinning, MinDose_mGyBinning, MaxDose_uGyBinning, NBins_Hits, MinHits, MaxHits);
    hDoseNucleus_HitsAlpha = new TH2D(histogramNameNucleus_HitsAlpha.c_str(), "Dose Delivered in Cell Nucleus per N Number of hits by alphas", NBins_mGyBinning, MinDose_mGyBinning, MaxDose_uGyBinning, NBins_Hits, MinHits, MaxHits);
    hDoseMembrane_HitsAlpha = new TH2D(histogramNameMembrane_HitsAlpha.c_str(), "Dose Delivered in Cell Membrane per N Number of hits by alphas", NBins_mGyBinning, MinDose_mGyBinning, MaxDose_uGyBinning, NBins_Hits, MinHits, MaxHits);
    hDoseCytoplasm_HitsAlpha = new TH2D(histogramNameCytoplasm_HitsAlpha.c_str(), "Dose Delivered in Cell Cytoplasm per N Number of hits by alphas", NBins_mGyBinning, MinDose_mGyBinning, MaxDose_uGyBinning, NBins_Hits, MinHits, MaxHits);


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


class AddEnergyDepositionHistograms{
    public:
        void AddHistograms(EnergyDepositionHistograms& h1, EnergyDepositionHistograms& h2);
};


void AddEnergyDepositionHistograms::AddHistograms(EnergyDepositionHistograms& h1, EnergyDepositionHistograms& h2)
{
         // Histogram for total energy deposited in one nuclei per number of cells
        h1.hEnergyDepsNucleus_eVBinning->Add(h2.hEnergyDepsNucleus_eVBinning);
        h1.hEnergyDepsNucleus_keVBinning->Add(h2.hEnergyDepsNucleus_keVBinning);
        h1.hEnergyDepsNucleus_FromSolution->Add(h2.hEnergyDepsNucleus_FromSolution);
        h1.hEnergyDepsNucleus_FromMembrane->Add(h2.hEnergyDepsNucleus_FromMembrane);
        h1.hEnergyDepsNucleus_FromCytoplasm->Add(h2.hEnergyDepsNucleus_FromCytoplasm);

        // Histogram for total energy deposited in one membrane per number of cells
        h1.hEnergyDepsMembrane_eVBinning->Add(h2.hEnergyDepsMembrane_eVBinning);
        h1.hEnergyDepsMembrane_keVBinning->Add(h2.hEnergyDepsMembrane_keVBinning);
        h1.hEnergyDepsMembrane_FromSolution->Add(h2.hEnergyDepsMembrane_FromSolution);
        h1.hEnergyDepsMembrane_FromMembrane->Add(h2.hEnergyDepsMembrane_FromMembrane);
        h1.hEnergyDepsMembrane_FromCytoplasm->Add(h2.hEnergyDepsMembrane_FromCytoplasm);


        // Histogram for total energy deposited in one cytoplasm per number of cells
        h1.hEnergyDepsCytoplasm_eVBinning->Add(h2.hEnergyDepsCytoplasm_eVBinning);
        h1.hEnergyDepsCytoplasm_keVBinning->Add(h2.hEnergyDepsCytoplasm_keVBinning);
        h1.hEnergyDepsCytoplasm_FromSolution->Add(h2.hEnergyDepsCytoplasm_FromSolution);
        h1.hEnergyDepsCytoplasm_FromMembrane->Add(h2.hEnergyDepsCytoplasm_FromMembrane);
        h1.hEnergyDepsCytoplasm_FromCytoplasm->Add(h2.hEnergyDepsCytoplasm_FromCytoplasm);

        // Histogram for total energy deposited in one cell per number of cells
        h1.hEnergyDepsTotalCell_eVBinning->Add(h2.hEnergyDepsTotalCell_eVBinning);
        h1.hEnergyDepsTotalCell_keVBinning->Add(h2.hEnergyDepsTotalCell_keVBinning);
        h1.hEnergyDepsTotalCell_FromSolution->Add(h2.hEnergyDepsTotalCell_FromSolution);
        h1.hEnergyDepsTotalCell_FromMembrane->Add(h2.hEnergyDepsTotalCell_FromMembrane);
        h1.hEnergyDepsTotalCell_FromCytoplasm->Add(h2.hEnergyDepsTotalCell_FromCytoplasm);

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

        // Histogram for total energy deposited in one membrane per number of cells
        h1.hDoseMembrane_uGyBinning->Add(h2.hDoseMembrane_uGyBinning);
        h1.hDoseMembrane_mGyBinning->Add(h2.hDoseMembrane_mGyBinning);
        h1.hDoseMembrane_FromSolution->Add(h2.hDoseMembrane_FromSolution);
        h1.hDoseMembrane_FromMembrane->Add(h2.hDoseMembrane_FromMembrane);
        h1.hDoseMembrane_FromCytoplasm->Add(h2.hDoseMembrane_FromCytoplasm);


        // Histogram for total energy deposited in one cytoplasm per number of cells
        h1.hDoseCytoplasm_uGyBinning->Add(h2.hDoseCytoplasm_uGyBinning);
        h1.hDoseCytoplasm_mGyBinning->Add(h2.hDoseCytoplasm_mGyBinning);
        h1.hDoseCytoplasm_FromSolution->Add(h2.hDoseCytoplasm_FromSolution);
        h1.hDoseCytoplasm_FromMembrane->Add(h2.hDoseCytoplasm_FromMembrane);
        h1.hDoseCytoplasm_FromCytoplasm->Add(h2.hDoseCytoplasm_FromCytoplasm);

        // Histogram for total energy deposited in one cell per number of cells
        h1.hDoseTotalCell_uGyBinning->Add(h2.hDoseTotalCell_uGyBinning);
        h1.hDoseTotalCell_mGyBinning->Add(h2.hDoseTotalCell_mGyBinning);
        h1.hDoseTotalCell_FromSolution->Add(h2.hDoseTotalCell_FromSolution);
        h1.hDoseTotalCell_FromMembrane->Add(h2.hDoseTotalCell_FromMembrane);
        h1.hDoseTotalCell_FromCytoplasm->Add(h2.hDoseTotalCell_FromCytoplasm);

        // Histograms for number of hits by alpha particles
        h1.hDoseTotalCell_HitsAlpha->Add(h2.hDoseTotalCell_HitsAlpha);
        h1.hDoseNucleus_HitsAlpha->Add(h2.hDoseNucleus_HitsAlpha);
        h1.hDoseMembrane_HitsAlpha->Add(h2.hDoseMembrane_HitsAlpha);
        h1.hDoseCytoplasm_HitsAlpha->Add(h2.hDoseCytoplasm_HitsAlpha);
}



//------------------–----------
EnergyDepositionHistograms AnalyzeHistogramsFromSimulation(DecayDynamics decayDynamicsInstance, int numberIterations)
{
    std::cout << "-----------------------" << std::endl;
    std::cout << "Analyzing Histograms for : " << decayDynamicsInstance.GetCellLine() << ", Activity: " << decayDynamicsInstance.GetActivity() << std::endl;
    //------------------–----------
    // Loading decay dynamics
    double numberDecays212PbInSolution1hTo2h = decayDynamicsInstance.GetNumberDecaysInSolutionFirstHour()*decayDynamicsInstance.GetVolumeRatio();
    double numberDecays212PbInMembrane1hTo26h = decayDynamicsInstance.GetNumberDecaysInMembraneTotalTime()*decayDynamicsInstance.GetVolumeRatio();
    double numberDecays212PbInCytoplasm1hTo26h = decayDynamicsInstance.GetNumberDecaysInCytoplasmTotalTime()*decayDynamicsInstance.GetVolumeRatio();

    // double numberDecays212PbInSolution1hTo2h = 100;
    // double numberDecays212PbInMembrane1hTo26h = 100;
    // double numberDecays212PbInCytoplasm1hTo26h = 100;


    double MaxEnergyCytoplasm = 100.;




    //------------------–----------
    //  Function for filling histograms
    auto MakeHistogramOneIteration = [&](std::string filepathSolutionSim_i, std::string filepathMembraneSim_i, std::string filepathCytoplasmSim_i, int itNum)
    // auto FillHistograms = [&](TChain* chainSolutionSim, TChain* chainMembraneSim, TChain* chainCytoplasmSim)
    {

        //------------------–----------
        // Generating empty histograms
        EnergyDepositionHistograms energyDepHistograms = EnergyDepositionHistograms(MaxEnergyCytoplasm, itNum);
        energyDepHistograms.GenerateEmptyHistograms(decayDynamicsInstance);


        //------------------–----------
        // Opening TTree files and creating TTreeReaders

        // Reader for solution simulation
        std::shared_ptr<TFile> myFileSolutionSim(TFile::Open(filepathSolutionSim_i.c_str(), "READ"));
        auto treeSolutionSim = myFileSolutionSim->Get<TTree>("B4");
        TTreeReader myReaderSolutionSim(treeSolutionSim);


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
            // Vector to store <cellID, trackID, volume hit> for particle hits to a cell for one event
            std::vector<std::tuple<int,int,int>> particleHitsInfoForEventVec;

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
                                    storedCellHits[ii].AddEnergyDeposition(energyDepsSolutionSim[i],volumeTypesSolutionSim[i],decayOriginSolution, particleTypeSolutionSim[i]);

                                    // Update boolean
                                    cellAlreadyHit = true;


                                    //---------------------------------------------
                                    // boolean: true if the particle track has already been stored for a specific cell ID in this specific component of the cell
                                    bool trackAlreadyRegisteredInThisCellComponent = false;

                                    // boolean: true if it is the first time the track has been counted in a cell
                                    bool firstTimeCountingTrackInCell = true;

                                    // Looping over all previously stored track IDs and corresponding cell IDs for this event
                                    for(int iii=0; iii<particleHitsInfoForEventVec.size(); iii++)
                                    {

                                        // If track ID is already stored for this cell ID, do not update particle hit counts for this cell ID
                                        if(cellIDsSolutionSim[i] == std::get<0>(particleHitsInfoForEventVec[iii]))
                                        {
                                            if(trackIDSolutionSim[i]==std::get<1>(particleHitsInfoForEventVec[iii]))
                                            {
                                                // Track has been counted before in this cell, just in a different component of the cell
                                                firstTimeCountingTrackInCell = false;

                                                if(volumeTypesSolutionSim[i]==get<2>(particleHitsInfoForEventVec[iii]))
                                                {
                                                    // Track already counted in this cell and cell component
                                                    trackAlreadyRegisteredInThisCellComponent = true;
                                                }
                                            }
                                        }
                                    }
                                    // If track ID has not already been registered or not been registered for this cell ID
                                    if(!trackAlreadyRegisteredInThisCellComponent)
                                    {

                                        // If alpha particle type
                                        if(particleTypeSolutionSim[i]==1000020040)
                                        {
                                            storedCellHits[ii].HitByAlphaParticle(volumeTypesSolutionSim[i],firstTimeCountingTrackInCell);
                                        }

                                        // Save the hit and cellID
                                        std::tuple<int,int,int> particleHit = std::make_tuple(cellIDsSolutionSim[i],trackIDSolutionSim[i],volumeTypesSolutionSim[i]);
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
                                aNewCellHit.AddEnergyDeposition(energyDepsSolutionSim[i], volumeTypesSolutionSim[i],decayOriginSolution,particleTypeSolutionSim[i]);
                                storedCellHits.push_back(aNewCellHit);


                                // If alpha particle
                                if(particleTypeSolutionSim[i]==1000020040)
                                {
                                    aNewCellHit.HitByAlphaParticle(volumeTypesSolutionSim[i],true);
                                }

                                std::tuple<int,int,int> particleHit = std::make_tuple(cellIDsSolutionSim[i],trackIDSolutionSim[i],volumeTypesSolutionSim[i]);
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
                std::vector<std::tuple<int,int,int>> particleHitsInfoForEventVec;

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

                                        storedCellHits[ii].AddEnergyDeposition(energyDepsMembraneSim[i],volumeTypesMembraneSim[i],decayOriginMembrane, particleTypeMembraneSim[i]);

                                        // Update boolean
                                        cellAlreadyHit = true;

                                        //--------------------------------------
                                        // boolean: true if the particle track has already been stored for a specific cell ID
                                        bool trackAlreadyRegisteredInThisCellComponent = false;
                                        bool firstTimeCountingTrackInCell = true;

                                        // Looping over all previously stored track IDs and corresponding cell IDs for this event
                                        for(int iii=0; iii<particleHitsInfoForEventVec.size(); iii++)
                                        {

                                            // If track ID is already stored for this cell ID, do not update particle hit counts for this cell ID
                                            if(cellIDsMembraneSim[i] == std::get<0>(particleHitsInfoForEventVec[iii]))
                                            {
                                                if(trackIDMembraneSim[i]==std::get<1>(particleHitsInfoForEventVec[iii]))
                                                {
                                                    firstTimeCountingTrackInCell = false;
                                                    if(volumeTypesMembraneSim[i]==std::get<2>(particleHitsInfoForEventVec[iii]))
                                                    {
                                                        trackAlreadyRegisteredInThisCellComponent = true;
                                                    }
                                                }
                                            }
                                        }
                                        // If track ID has not already been registered or not been registered for this cell ID
                                        if(!trackAlreadyRegisteredInThisCellComponent)
                                        {

                                            // If alpha particle type
                                            if(particleTypeMembraneSim[i]==1000020040)
                                            {
                                                storedCellHits[ii].HitByAlphaParticle(volumeTypesMembraneSim[i],firstTimeCountingTrackInCell);
                                            }

                                            // Save the hit and cellID
                                            std::tuple<int,int,int> particleHit = std::make_tuple(cellIDsMembraneSim[i],trackIDMembraneSim[i],volumeTypesMembraneSim[i]);
                                            particleHitsInfoForEventVec.push_back(particleHit);

                                        }
                                    }
                                }

                                // Register a new cell hit, and add energy deposition
                                if(!cellAlreadyHit)
                                {
                                    // std::cout << "Interaction volume " << volumeTypesSolutionSim[i] << std::endl;

                                    CellHit aNewCellHit = CellHit(cellIDsMembraneSim[i]);
                                    aNewCellHit.AddEnergyDeposition(energyDepsMembraneSim[i], volumeTypesMembraneSim[i], decayOriginMembrane,particleTypeMembraneSim[i]);
                                    storedCellHits.push_back(aNewCellHit);

                                    // If alpha particle
                                    if(particleTypeMembraneSim[i]==1000020040)
                                    {
                                        aNewCellHit.HitByAlphaParticle(volumeTypesMembraneSim[i],true);
                                    }

                                    std::tuple<int,int,int> particleHit = std::make_tuple(cellIDsMembraneSim[i],trackIDMembraneSim[i],volumeTypesMembraneSim[i]);
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
                std::vector<std::tuple<int,int,int>> particleHitsInfoForEventVec;

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

                                        storedCellHits[ii].AddEnergyDeposition(energyDepsCytoplasmSim[i],volumeTypesCytoplasmSim[i],decayOriginCytoplasm,particleTypeCytoplasmSim[i]);

                                        // Update boolean
                                        cellAlreadyHit = true;

                                        // boolean: true if the particle track has already been stored for a specific cell ID
                                        bool trackAlreadyRegisteredInThisCellComponent = false;
                                        bool firstTimeCountingTrackInCell = true;

                                        // Looping over all previously stored track IDs and corresponding cell IDs for this event
                                        for(int iii=0; iii<particleHitsInfoForEventVec.size(); iii++)
                                        {

                                            // If track ID is already stored for this cell ID, do not update particle hit counts for this cell ID
                                            if(cellIDsCytoplasmSim[i] == std::get<0>(particleHitsInfoForEventVec[iii]))
                                            {
                                                if(trackIDCytoplasmSim[i]==std::get<1>(particleHitsInfoForEventVec[iii]))
                                                {
                                                    firstTimeCountingTrackInCell = false;
                                                    if(volumeTypesCytoplasmSim[i]==std::get<2>(particleHitsInfoForEventVec[iii]))
                                                    {
                                                        trackAlreadyRegisteredInThisCellComponent = true;
                                                    }
                                                }
                                            }
                                        }
                                        // If track ID has not already been registered or not been registered for this cell ID
                                        if(!trackAlreadyRegisteredInThisCellComponent)
                                        {

                                            // If alapha particle type
                                            if(particleTypeCytoplasmSim[i]==1000020040)
                                            {
                                                storedCellHits[ii].HitByAlphaParticle(volumeTypesCytoplasmSim[i],firstTimeCountingTrackInCell);
                                                // std::cout << "New track : "<< cellIDsCytoplasmSim[i] << " , " << trackIDCytoplasmSim[i] << " , " << particleTypeCytoplasmSim[i] << std::endl;
                                            }

                                            // Save the hit and cellID
                                            std::tuple<int,int,int> particleHit = std::make_tuple(cellIDsCytoplasmSim[i],trackIDCytoplasmSim[i],volumeTypesCytoplasmSim[i]);
                                            particleHitsInfoForEventVec.push_back(particleHit);

                                        }
                                    }
                                }
                                // Register a new cell hit, and add energy deposition
                                if(!cellAlreadyHit)
                                {
                                    // std::cout << "Interaction volume " << volumeTypesSolutionSim[i] << std::endl;

                                    CellHit aNewCellHit = CellHit(cellIDsCytoplasmSim[i]);
                                    aNewCellHit.AddEnergyDeposition(energyDepsCytoplasmSim[i], volumeTypesCytoplasmSim[i],decayOriginCytoplasm,particleTypeCytoplasmSim[i]);
                                    storedCellHits.push_back(aNewCellHit);

                                    // If alpha particle
                                    if(particleTypeCytoplasmSim[i]==1000020040)
                                    {
                                        aNewCellHit.HitByAlphaParticle(volumeTypesCytoplasmSim[i],true);
                                    }

                                    std::tuple<int,int,int> particleHit = std::make_tuple(cellIDsCytoplasmSim[i],trackIDCytoplasmSim[i],volumeTypesCytoplasmSim[i]);
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
            // energyDepCytoplasmVec.push_back(storedCellHits[i].GetEnergyDepositionCytoplasm());
            energyDepHistograms.AddCellHitsToHistograms(storedCellHits[i]);


        }
        return energyDepHistograms;
    };


    // double maxE = *std::max_element(energyDepCytoplasmVec.begin(), energyDepCytoplasmVec.end());


    // std::string filepathSimulationOutput = "../../GEANT4Simulations/OutputFromSaga/";
    // std::string filepathSimulationOutput = "../../GEANT4Simulations/CellDamageSimulation-build/";
    std::string filepathSimulationOutput = "/Volumes/SamsungT7/OutputFromSaga/";

    std::string filepathSolutionIteration_i;
    std::string filepathMembraneIteration_i;
    std::string filepathCytoplasmIteration_i;


    EnergyDepositionHistograms histMain = EnergyDepositionHistograms(100., 0);
    histMain.GenerateEmptyHistograms(decayDynamicsInstance);

    AddEnergyDepositionHistograms addHistograms;


    //------------------–----------
    // Filling histgrams
    for(int i=0; i<numberIterations; i++)
    {
        std::cout << "Iteration number : " << i+1 << std::endl;
        filepathSolutionIteration_i = filepathSimulationOutput + "Output_Pb212_" + decayDynamicsInstance.GetCellLine() + "_Activity" + std::to_string(decayDynamicsInstance.GetActivity()) + "kBq_Solution_Thread_" + std::to_string(i) + ".root";
        filepathMembraneIteration_i = filepathSimulationOutput + "Output_Pb212_" + decayDynamicsInstance.GetCellLine() + "_Activity" + std::to_string(decayDynamicsInstance.GetActivity()) + "kBq_Membrane_Thread_" + std::to_string(i) + ".root";
        filepathCytoplasmIteration_i = filepathSimulationOutput + "Output_Pb212_" + decayDynamicsInstance.GetCellLine() + "_Activity" + std::to_string(decayDynamicsInstance.GetActivity()) + "kBq_Cytoplasm_Thread_" + std::to_string(i) + ".root";

        EnergyDepositionHistograms histIteration = MakeHistogramOneIteration(filepathSolutionIteration_i, filepathMembraneIteration_i, filepathCytoplasmIteration_i, i+1);
        addHistograms.AddHistograms(histMain,histIteration);
        // histIteration.ResetHistograms();
    }


    //------------------–----------
    // Scaling histograms
    double scalingFactor = decayDynamicsInstance.GetNumberCells()*((double) numberIterations);
    histMain.ScaleHistograms(1./scalingFactor);



    //------------------–----------
    return histMain;

}

void mainAnalysisCode()
{

    int numberIterations = 3;



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
    // DecayDynamics decays_A210kBq_PC3_Flu = DecayDynamics(250,0.,0., "PC3_Flu");




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
    // decays_A210kBq_PC3_Flu.LoadDataFromMathematicaCalculations(mathematicaOutput.c_str());



    //-------------------------------------
    // Creating energy deposition histograms

    // EnergyDepositionHistograms Hist_A150kBq_PC3_Flu = AnalyzeHistogramsFromSimulation(decays_A150kBq_PC3_Flu, numberIterations);
    // auto output = new TFile("/Volumes/SamsungT7/OutputFromAnalysis/Output_PC3_Flu_150kBq.root", "RECREATE");
    // Hist_A150kBq_PC3_Flu.WriteHistogramsToFile();
    // output->Write();
    // output->Close();

    // EnergyDepositionHistograms Hist_A150kBq_PC3_PIP = AnalyzeHistogramsFromSimulation(decays_A150kBq_PC3_PIP, numberIterations);
    // auto output = new TFile("/Volumes/SamsungT7/OutputFromAnalysis/Output_PC3_PIP_150kBq.root", "RECREATE");
    // Hist_A150kBq_PC3_PIP.WriteHistogramsToFile();
    // output->Write();
    // output->Close();

    EnergyDepositionHistograms Hist_A150kBq_C4_2 = AnalyzeHistogramsFromSimulation(decays_A150kBq_C4_2, numberIterations);
    auto output = new TFile("/Volumes/SamsungT7/OutputFromAnalysis/Output_C4_2_150kBq.root", "RECREATE");
    Hist_A150kBq_C4_2.WriteHistogramsToFile();
    output->Write();
    output->Close();
}


