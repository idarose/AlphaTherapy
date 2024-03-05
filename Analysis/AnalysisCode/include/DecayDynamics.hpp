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

#ifndef DECAYDYNAMICS_HPP
#define DECAYDYNAMICS_HPP

//------------------–----------
class DecayDynamics
{
    //----------------------
    //  Class to store calculated decay dynamics data from the Mathematica calculations.
    //  Assumes the output file from Mathematica is structured in a specific way

    public:
        DecayDynamics(int activitySample_in, std::string cellLine_in);

        std::vector<double> ReadFileFromMathematica(std::string filename);
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

        // double U0PerCell; // Number radionuclides absorbed per cell
        // double U0InternalizedPerCell; // Number of absorbed radionuclides that are found in cytoplasm
        int activitySample;  // Given in kBq/1mL
        std::string cellLine;

        // Volume ratio between simulation and actual sample volume
        double volumeRatio;

        // Number of cells in simualted volume
        double numberCells;


};

#endif // DECAYDYNAMICS_HPP