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
        DecayDynamics(int activitySample_in, std::string cellLine_in, std::string cellGeometry_in);

        std::vector<double> ReadFileFromMathematica(std::string filename);
        void LoadDataFromMathematicaCalculations(std::string filepathToMathematicaOutput);

        double GetNumberDecaysInSolutionFirstHour(){return numberDecays212PbInSolution1hTo2h;};
        double GetNumberDecaysInMembraneTotalTime(){return numberDecays212PbInMembrane1h2To26h;};
        double GetNumberDecaysInCytoplasmTotalTime(){return numberDecays212PbInCytoplasm1hTo26h;};

        double GetVolumeRatio(){return volumeRatio;};
        double GetNumberCells(){return numberCells;};

        double GetMassNucleus(){return massNucleus;};
        double GetMassMembrane(){return massMembrane;};
        double GetMassCytoplasm(){return massCytoplasm;};
        double GetMassCell(){return massCell;};


        int GetActivity(){return activitySample;};
        std::string GetCellLine(){return cellLine;};
        std::string GetCellGeometry(){return cellGeometry;};

    private:

        int activitySample;

        // These number of decays are calculated in Mathematica for a sample of 0.2mL in volume
        double numberDecays212PbInSolution1hTo2h;
        double numberDecays212PbInMembrane1h2To26h;
        double numberDecays212PbInCytoplasm1hTo26h;

        //------------------------------
        std::string cellGeometry;
        std::string cellLine;

        // Volume ratio between simulation and actual sample volume
        double volumeRatio;

        // Number of cells in simualted volume
        double numberCells;

        double massNucleus;
        double massCell;
        double massMembrane;
        double massCytoplasm;


};

#endif // DECAYDYNAMICS_HPP