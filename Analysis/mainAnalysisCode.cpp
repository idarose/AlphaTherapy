#include <TFile.h>
#include <TH1.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TMath.h>
#include <fstream>
#include <sstream>

#include <vector>
#include <iostream>

#include <tuple>


//------------------–----------
std::vector<double> ImportDataStringDouble(std::string filename)
{
    //----------------------
    //  Imports data from file where the first column contains data in string form
    //  and second column contains data in double form. Returns the columns of double
    //  values in the form of a vector
    std::vector<double> input_data;

    // std::fstream myfile;
    // myfile.open(filename);

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
        DecayDynamics(int activitySample_in, double U0InternalizedPerCell_in, double U0SurfaceBoundPerCell_in, std::string cellLine_in);

        void LoadDataFromMathematicaCalculations(std::string filepathToMathematicaOutput);

        double GetNumberDecaysInSolutionFirstHour(){return numberDecays212PbInSolutionFirstHour;};
        double GetNumberDecaysInMembraneTotalTime(){return numberDecays212PbInMembraneTotalTime;};
        double GetNumberDecaysInCytoplasmTotalTime(){return numberDecays212PbInCytoplasmTotalTime;};

        int GetActivity(){return activitySample;};
        std::string GetCellLine(){return cellLine;};

    private:

        // These number of decays are calculated in Mathematica for a sample of 0.2mL in volume
        double numberDecays212PbInSolutionFirstHour;
        double numberDecays212PbInMembraneTotalTime;
        double numberDecays212PbInCytoplasmTotalTime;

        double U0InternalizedPerCell;
        double U0SurfaceBoundPerCell;
        int activitySample;  // Given in kBq/1mL
        std::string cellLine;
};


DecayDynamics::DecayDynamics(int activitySample_in, double U0InternalizedPerCell_in, double U0SurfaceBoundPerCell_in, std::string cellLine_in)
{
    activitySample = activitySample_in;
    U0InternalizedPerCell = U0InternalizedPerCell_in;
    U0SurfaceBoundPerCell = U0SurfaceBoundPerCell_in;
    cellLine = cellLine_in;

}



void DecayDynamics::LoadDataFromMathematicaCalculations(std::string filepathToMathematicaOutput)
{
    //----------------------
    // Loads data from calculations, assuming output files are structured as shown in the file "212PbDecayDynamics.nb"

    std::string filepathSolutionData = filepathToMathematicaOutput + "/" + cellLine + "/Solution/Activity_" + std::to_string(activitySample) + "kBq/Decays.dat";

    // Importing data for decays occuring in solution
    std::vector<double> decayDataSolution = ImportDataStringDouble(filepathSolutionData);
    numberDecays212PbInSolutionFirstHour = decayDataSolution[0];

    // Importing data for decays occuring in cells
    std::string filepathCellData = filepathToMathematicaOutput + "/" + cellLine + "/Cells/Activity_" + std::to_string(activitySample) + "kBq/Decays.dat";
    std::vector<double> decayDataCells = ImportDataStringDouble(filepathCellData);

    numberDecays212PbInMembraneTotalTime = decayDataCells[7];
    numberDecays212PbInCytoplasmTotalTime = decayDataCells[14];

}



//------------------–----------
class CellHit
{
    //----------------------
    //  Class to store information on the energy depositions in ONE specific cell.
    //  When a new cell is "registered" as hit by a particle, a new instance of this class is created.
    //  Any further energy depositions registered in this cell then gets added to the previously
    //  registered energy deposition. The energy deposition is stored according to which "part" of
    //  the cell is hit by the.

    public:
        CellHit(int cellID_in);

        void AddEnergyDeposition(double energyDep_in, double volumeType_in);

        double GetSumEnergyDepositions(){return energyDepMembrane + energyDepCytoplasm + energyDepNucleus;};

        double GetEnergyDepositionMembrane(){return energyDepMembrane;};
        double GetEnergyDepositionCytoplasm(){return energyDepCytoplasm;};
        double GetEnergyDepositionNucleus(){return energyDepNucleus;};

        int GetCellID(){return cellID;};

    private:
        int cellID;
        int volumeType;

        double energyDepMembrane;
        double energyDepCytoplasm;
        double energyDepNucleus;

        double sumEnergyDepositions;
};


//------------------–----------
CellHit::CellHit(int cellID_in)
{
    cellID = cellID_in;
    energyDepMembrane = 0.0;
    energyDepCytoplasm = 0.0;
    energyDepNucleus = 0.0;
}


//------------------–----------
void CellHit::AddEnergyDeposition(double energyDep_in, double volumeType_in)
{
    if(volumeType_in==1)
    {
        energyDepMembrane += energyDep_in;
    }
    else if(volumeType_in==2)
    {
        energyDepCytoplasm += energyDep_in;
    }
    else if(volumeType_in==3)
    {
        energyDepNucleus += energyDep_in;
    }
}



//------------------–----------
class EnergyDepositionHistograms
{
    public:
        EnergyDepositionHistograms(int NBins_in, double EMin_in, double EMax_in);

        void GenerateEmptyHistograms(DecayDynamics decayDynamicsInstance);
        void AddCellHitsToHistograms(CellHit cellHits);
        void ScaleHistograms(double factor);
        void WriteHistogramsToFile();

        TH1D* GetEnergyDepNucleusHist(){return hEnergyDepsNucleus;};
        TH1D* GetEnergyDepMembraneHist(){return hEnergyDepsMembrane;};
        TH1D* GetEnergyDepCytoplasmHist(){return hEnergyDepsCytoplasm;};
        TH1D* GetEnergyDepCellTotalHist(){return hEnergyDepsCellTotal;};

    private:
        int NBins;
        double EMin;
        double EMax;


        // Histogram for total energy deposited in one nuclei per number of cells
        TH1D *hEnergyDepsNucleus;

        // Histogram for total energy deposited in one membrane per number of cells
        TH1D *hEnergyDepsMembrane;

        // Histogram for total energy deposited in one cytoplasm per number of cells
        TH1D *hEnergyDepsCytoplasm;

        // Histogram for total energy deposited in one cell per number of cells
        TH1D *hEnergyDepsCellTotal;

        // Histogram for total energy deposited in both membrane and cytoplasm of one cell, per number of cells
        TH1D *hEnergyDepsMembraneAndCytoplasm;

        // Histogram for total energy deposited in both membrane and nucleus of on cell per number of cells
        TH1D *hEnergyDepsMembraneAndNucleus;

        // Histogram for total energy deposited in both nucleus and cytoplasm of one cell per number of cells
        TH1D *hEnergyDepsNucleusAndCytoplasm;
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
    std::string histogramNameCellTotal = generalHistogramName + "CellTotal";
    std::string histogramNameMembraneAndCytoplasm = generalHistogramName + "MembraneAndCytoplasm";
    std::string histogramNameMembraneAndNucleus = generalHistogramName + "MembraneAndNucleus";
    std::string histogramNameNucleusAndCytoplasm = generalHistogramName + "NucleusAndCytoplasm";

    hEnergyDepsNucleus = new TH1D(histogramNameNucleus.c_str(), "Energy Deposition in Cell Nucleus / Number of Cells", NBins, EMin, EMax);
    hEnergyDepsMembrane = new TH1D(histogramNameMembrane.c_str(), "Energy Depsition in Cell Membrane / Number of Cells", NBins, EMin, EMax);
    hEnergyDepsCytoplasm = new TH1D(histogramNameCytoplasm.c_str(), "Energy Depsition in Cell Cytoplasm / Number of Cells", NBins, EMin, EMax);
    hEnergyDepsCellTotal = new TH1D(histogramNameCellTotal.c_str(), "Energy Depsition in Cell / Number of Cells", NBins, EMin, EMax);
    hEnergyDepsMembraneAndCytoplasm = new TH1D(histogramNameMembraneAndCytoplasm.c_str(), "Energy Depsition in Cell Membrane and Cell Cytoplasm / Number of Cells", NBins, EMin, EMax);
    hEnergyDepsMembraneAndNucleus = new TH1D(histogramNameMembraneAndNucleus.c_str(), "Energy Depsition in Cell Membrane and Cell Nucleus / Number of Cells", NBins, EMin, EMax);
    hEnergyDepsNucleusAndCytoplasm = new TH1D(histogramNameNucleusAndCytoplasm.c_str(), "Energy Depsition in Cell Nucleus and Cell Cytoplasm / Number of Cells", NBins, EMin, EMax);
}

//------------------–----------
void EnergyDepositionHistograms::AddCellHitsToHistograms(CellHit cellHits)
{
    if(cellHits.GetEnergyDepositionMembrane()>0.0)
        {hEnergyDepsMembrane->Fill(cellHits.GetEnergyDepositionMembrane());}
    if(cellHits.GetEnergyDepositionCytoplasm()>0.0)
        {hEnergyDepsCytoplasm->Fill(cellHits.GetEnergyDepositionCytoplasm());}
    if(cellHits.GetEnergyDepositionNucleus()>0.0)
        {hEnergyDepsNucleus->Fill(cellHits.GetEnergyDepositionNucleus());}
    if(cellHits.GetSumEnergyDepositions()>0.0)
        {hEnergyDepsCellTotal->Fill(cellHits.GetSumEnergyDepositions());}
    if((cellHits.GetEnergyDepositionMembrane()+cellHits.GetEnergyDepositionCytoplasm())>0.0)
        {hEnergyDepsMembraneAndCytoplasm->Fill(cellHits.GetEnergyDepositionMembrane()+cellHits.GetEnergyDepositionCytoplasm());}
    if((cellHits.GetEnergyDepositionMembrane()+cellHits.GetEnergyDepositionNucleus())>0.0)
    {hEnergyDepsMembraneAndNucleus->Fill(cellHits.GetEnergyDepositionMembrane()+cellHits.GetEnergyDepositionNucleus());}
    if((cellHits.GetEnergyDepositionNucleus()+cellHits.GetEnergyDepositionCytoplasm())>0.0){hEnergyDepsNucleusAndCytoplasm->Fill(cellHits.GetEnergyDepositionNucleus()+cellHits.GetEnergyDepositionCytoplasm());}
}

//------------------–----------
void EnergyDepositionHistograms::ScaleHistograms(double factor)
{


    // double binCon = hEnergyDepsMembrane->GetBinContent(723);
    // double error = std::pow(binCon, 1.0/2.0);
    // std::cout << "Bin content: " << binCon << " Error: " << error << std::endl;
    // std::cout << "Bin error (%): " << 100.0*error/binCon  << std::endl;

    hEnergyDepsMembrane->Scale(factor);
    hEnergyDepsCytoplasm->Scale(factor);
    hEnergyDepsNucleus->Scale(factor);
    hEnergyDepsCellTotal->Scale(factor);
    hEnergyDepsMembraneAndCytoplasm->Scale(factor);
    hEnergyDepsMembraneAndNucleus->Scale(factor);
    hEnergyDepsNucleusAndCytoplasm->Scale(factor);

    // std:: cout << "Scaling factor = " << factor << " Scaling factor * error = " << error*factor << std::endl;
    // double binCon2 = hEnergyDepsMembrane->GetBinContent(723);
    // double error2 = hEnergyDepsMembrane->GetBinError(723);
    // double error2 = std::pow(binCon2, 1.0/2.0);
    // std::cout << "Bin content: " << binCon2 << " Error: " << error2 << std::endl;
    // std::cout << "Bin error2 (%): " << 100.0*error2/binCon2 << std::endl;

    // std::cout << "error_scaled / error_Unscaled = " << error2/error << std::endl;

}

//------------------–----------
void EnergyDepositionHistograms::WriteHistogramsToFile()
{
    // outputFile->cd();

    hEnergyDepsMembrane->Write();
    hEnergyDepsCytoplasm->Write();
    hEnergyDepsNucleus->Write();
    hEnergyDepsCellTotal->Write();
    hEnergyDepsMembraneAndCytoplasm->Write();
    hEnergyDepsMembraneAndNucleus->Write();
    hEnergyDepsNucleusAndCytoplasm->Write();

    //------------------–----------
    hEnergyDepsMembrane->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    hEnergyDepsMembrane->GetYaxis()->SetTitle("Hits / Cell");
    hEnergyDepsCytoplasm->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    hEnergyDepsCytoplasm->GetYaxis()->SetTitle("Hits / Cell");
    hEnergyDepsNucleus->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    hEnergyDepsNucleus->GetYaxis()->SetTitle("Hits / Cell");
    hEnergyDepsCellTotal->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    hEnergyDepsCellTotal->GetYaxis()->SetTitle("Hits / Cell");
    hEnergyDepsMembraneAndCytoplasm->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    hEnergyDepsMembraneAndCytoplasm->GetYaxis()->SetTitle("Hits / Cell");
    hEnergyDepsMembraneAndNucleus->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    hEnergyDepsMembraneAndNucleus->GetYaxis()->SetTitle("Hits / Cell");
    hEnergyDepsNucleusAndCytoplasm->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    hEnergyDepsNucleusAndCytoplasm->GetYaxis()->SetTitle("Hits / Cell");


}


//------------------–----------
EnergyDepositionHistograms MakeHistograms(DecayDynamics decayDynamicsInstance, int numberIterations, double volumeRatio, int numberCells)
{
    //------------------–----------
    // Loading decay dynamics
    double numberDecays212PbInSolutionFirstHour = decayDynamicsInstance.GetNumberDecaysInSolutionFirstHour()*volumeRatio;
    double numberDecays212PbInMembraneTotalTime = decayDynamicsInstance.GetNumberDecaysInMembraneTotalTime()*volumeRatio;
    double numberDecays212PbInCytoplasmTotalTime = decayDynamicsInstance.GetNumberDecaysInCytoplasmTotalTime()*volumeRatio;


    //------------------–----------
    int NBins = 2000000;
    double EMin = 0.0;
    double EMax = 20.0;

    //------------------–----------
    // Making energy deposition histogram instance
    EnergyDepositionHistograms energyDepHistograms = EnergyDepositionHistograms(NBins, EMin, EMax);
    energyDepHistograms.GenerateEmptyHistograms(decayDynamicsInstance);



    //------------------–----------
    // Opening TTree files and creating TTreeReaders

    // Reader for solution simulation
    // std::shared_ptr<TFile> myFileSolutionSim(TFile::Open("../GEANT4Simulations/OutputFromSaga/Output_212Pb_C4_2_Solution.root", "READ"));
    std::shared_ptr<TFile> myFileSolutionSim(TFile::Open("../GEANT4Simulations/B4aSolution-build/Output_212Pb_C4_2_Solution.root", "READ"));
    auto treeSolutionSim = myFileSolutionSim->Get<TTree>("B4");
    TTreeReader myReaderSolutionSim(treeSolutionSim);

    // Reader for membrane simulation
    // std::shared_ptr<TFile> myFileMembraneSim(TFile::Open("../GEANT4Simulations/OutputFromSaga/Output_212Pb_C4_2_Membrane.root", "READ"));
    std::shared_ptr<TFile> myFileMembraneSim(TFile::Open("../GEANT4Simulations/B4aMembrane-build/Output_212Pb_C4_2_Membrane.root", "READ"));
    auto treeMembraneSim = myFileMembraneSim->Get<TTree>("B4");
    TTreeReader myReaderMembraneSim(treeMembraneSim);

        // Reader for cytoplasm simulation
    // std::shared_ptr<TFile> myFileCytoplasmSim(TFile::Open("../GEANT4Simulations/OutputFromSaga/Output_212Pb_C4_2_Cytoplasm.root", "READ"));
    std::shared_ptr<TFile> myFileCytoplasmSim(TFile::Open("../GEANT4Simulations/B4aCytoplasm-build/Output_212Pb_C4_2_Cytoplasm.root", "READ"));
    auto treeCytoplasmSim = myFileCytoplasmSim->Get<TTree>("B4");
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

    // Membrane
    TTreeReaderArray<double> energyDepsMembraneSim(myReaderMembraneSim, "EnergyDeps");
    TTreeReaderArray<int> volumeTypesMembraneSim(myReaderMembraneSim, "VolumeTypes");
    TTreeReaderArray<int> cellIDsMembraneSim(myReaderMembraneSim, "CellIDs");
    TTreeReaderArray<double> kineticEnergyMembraneSim(myReaderMembraneSim, "KineticEnergy");
    TTreeReaderArray<int> particleTypeMembraneSim(myReaderMembraneSim, "ParticleType");
    TTreeReaderArray<double> interactionTimeMembraneSim(myReaderMembraneSim, "InteractionTime");
    TTreeReaderArray<double> firstInteractionTimeMembraneSim(myReaderSolutionSim, "FirstInteractionTime");
    TTreeReaderArray<int> firstInteractionVolumeMembraneSim(myReaderSolutionSim, "FirstInteractionVolume");


    // Cytoplasm
    TTreeReaderArray<double> energyDepsCytoplasmSim(myReaderCytoplasmSim, "EnergyDeps");
    TTreeReaderArray<int> volumeTypesCytoplasmSim(myReaderCytoplasmSim, "VolumeTypes");
    TTreeReaderArray<int> cellIDsCytoplasmSim(myReaderCytoplasmSim, "CellIDs");
    TTreeReaderArray<double> kineticEnergyCytoplasmSim(myReaderCytoplasmSim, "KineticEnergy");
    TTreeReaderArray<int> particleTypeCytoplasmSim(myReaderCytoplasmSim, "ParticleType");
    TTreeReaderArray<double> interactionTimeCytoplasmSim(myReaderCytoplasmSim, "InteractionTime");
    TTreeReaderArray<double> firstInteractionTimeCytoplasmSim(myReaderCytoplasmSim, "FirstInteractionTime");
    TTreeReaderArray<int> firstInteractionVolumeCytoplasmSim(myReaderCytoplasmSim, "FirstInteractionVolume");


    double entries = myReaderSolutionSim.GetEntries();

    double neededEventsSol = numberDecays212PbInSolutionFirstHour*2.8/0.98;
    double neededEventsCyt = numberDecays212PbInCytoplasmTotalTime*1/0.97;
    // std::cout << "Solution : Activity " << decayDynamicsInstance.GetActivity() << " Needed Events = " << neededEventsSol << std::endl;
    // std::cout << "Cytoplasm : Activity " << decayDynamicsInstance.GetActivity() << " Needed Events = " << neededEventsCyt << std::endl;

    // std::cout << "Activity : " << decayDynamicsInstance.GetActivity() << " Decays in solution " << numberDecays212PbInSolutionFirstHour << std::endl;

    //------------------–----------
    //  Function for filling histograms
    auto FillHistograms = [&]()
    {

        // Vector to store cell hits
        std::vector<CellHit> storedCellHits;

        //------------------–----------
        // Looping through data for decays occurring in solution in first hour

        // Counter variable
        int numberDecays212PbInSolutionFirstHour_counter = 0;

        while(myReaderSolutionSim.Next())
        {

            // bool for if the decay occurred in solution in first hour
            bool firstInteractionInSolutionInFirstHour = false;

            // Checking when and where first interaction occured
            if(firstInteractionVolumeSolutionSim[0]==0)
            {
                if(firstInteractionTimeSolutionSim[0]/3600. < 1.0)
                {

                    // Sett bool true
                    firstInteractionInSolutionInFirstHour = true;

                    // Update counter
                    numberDecays212PbInSolutionFirstHour_counter ++;
                }
            }

            // Sorting through interactions
            if(firstInteractionInSolutionInFirstHour)
            {
                // looping over all steps for one event/decay
                for(int i=0; i<energyDepsSolutionSim.GetSize(); i++)
                {

                    // Only add if actual energy deposition
                    if(energyDepsSolutionSim[i]>0.)
                    {

                        // Only add if interaction happened within first hour
                        if(interactionTimeSolutionSim[i]/3600.0 < 1.0)
                        {

                            // boolean: true if cell has not already been hit before in this event
                            // false if already hit this event
                            bool cellAlreadyHit = false;

                            // Loop over all previous cell hits stored
                            for(int ii=0; ii<storedCellHits.size(); ii++)
                            {

                                // Add energy deposition to cell that has already been hit before
                                if(cellIDsSolutionSim[i]==storedCellHits[ii].GetCellID())
                                {
                                    storedCellHits[ii].AddEnergyDeposition(energyDepsSolutionSim[i],volumeTypesSolutionSim[i]);

                                    // Update boolean
                                    cellAlreadyHit = true;
                                }
                            }

                            // Register a new cell hit, and add energy deposition
                            if(!cellAlreadyHit)
                            {
                                CellHit aNewCellHit = CellHit(cellIDsSolutionSim[i]);
                                aNewCellHit.AddEnergyDeposition(energyDepsSolutionSim[i], volumeTypesSolutionSim[i]);
                                storedCellHits.push_back(aNewCellHit);
                            }
                        }
                    }
                }
            }

            //------------------–----------
            // Break loop if number of decays have been reached
            if(numberDecays212PbInSolutionFirstHour_counter >= numberDecays212PbInSolutionFirstHour)
            {
                break;
            }
        }



        //------------------–----------
        // Looping through data for decays occuring in Membrane in 24 hours

        // Counter variable
        int numberDecays212PbInMembrane24h_counter = 0;


        while(myReaderMembraneSim.Next())
        {
            bool firstInteractionInMembraneFirst24h = false;

            if(firstInteractionTimeMembraneSim[0]/3600. < 24.)
            {
                firstInteractionInMembraneFirst24h = true;
                numberDecays212PbInMembrane24h_counter ++;
            }

            if(firstInteractionInMembraneFirst24h)
            {
                // looping over all steps for one event/decay
                for(int i=0; i<energyDepsMembraneSim.GetSize(); i++)
                {

                    // Only add if actual energy deposition
                    if(energyDepsMembraneSim[i]>0.)
                    {

                        // Only add if interaction happened within first 24 hours
                        if(interactionTimeMembraneSim[i]/3600.0 < 24.0)
                        {

                            // boolean: true if cell has not already been hit before in this event
                            // false if already hit this event
                            bool cellNotAlreadyHit = true;

                            // Loop over all previous cell hits stored
                            for(int ii=0; ii<storedCellHits.size(); ii++)
                            {

                                // Add energy deposition to cell that has already been hit before
                                if(cellIDsMembraneSim[i]==storedCellHits[ii].GetCellID())
                                {
                                    storedCellHits[ii].AddEnergyDeposition(energyDepsMembraneSim[i],volumeTypesMembraneSim[i]);

                                    // Update boolean
                                    cellNotAlreadyHit = false;
                                }
                            }

                            // Register a new cell hit, and add energy deposition
                            if(cellNotAlreadyHit)
                            {
                                CellHit aNewCellHit = CellHit(cellIDsMembraneSim[i]);
                                aNewCellHit.AddEnergyDeposition(energyDepsMembraneSim[i], volumeTypesMembraneSim[i]);
                                storedCellHits.push_back(aNewCellHit);
                            }
                        }
                    }
                }
            }

            // Break loop if number of decays have been reached
            if(numberDecays212PbInMembrane24h_counter >= numberDecays212PbInMembraneTotalTime)
            {
                break;
            }
        }



        //------------------–----------
        // Looping through data for decays occuring in Cytoplasm in 24 hours

        // Counter variable
        int numberDecays212PbInCytoplasm24h_counter = 0;


        while(myReaderCytoplasmSim.Next())
        {

            bool firstInteractionInCytoplasmFirst24h = false;

            if(firstInteractionVolumeCytoplasmSim[0] == 2)
            {
                if(firstInteractionTimeCytoplasmSim[0]/3600. < 24.)
                {
                    firstInteractionInCytoplasmFirst24h = true;
                    numberDecays212PbInCytoplasm24h_counter ++;
                }
            }

            if(firstInteractionInCytoplasmFirst24h)
            {
                // looping over all steps for one event/decay
                for(int i=0; i<energyDepsCytoplasmSim.GetSize(); i++)
                {

                    // Only add if actual energy deposition
                    if(energyDepsCytoplasmSim[i]>0.)
                    {

                        // Only add if interaction happened within first 24 hours
                        if(interactionTimeCytoplasmSim[i]/3600.0 < 24.0)
                        {

                            // boolean: true if cell has not already been hit before in this event
                            // false if already hit this event
                            bool cellNotAlreadyHit = true;

                            // Loop over all previous cell hits stored
                            for(int ii=0; ii<storedCellHits.size(); ii++)
                            {

                                // Add energy deposition to cell that has already been hit before
                                if(cellIDsCytoplasmSim[i]==storedCellHits[ii].GetCellID())
                                {
                                    storedCellHits[ii].AddEnergyDeposition(energyDepsCytoplasmSim[i],volumeTypesCytoplasmSim[i]);

                                    // Update boolean
                                    cellNotAlreadyHit = false;
                                }
                            }

                            // Register a new cell hit, and add energy deposition
                            if(cellNotAlreadyHit)
                            {
                                CellHit aNewCellHit = CellHit(cellIDsCytoplasmSim[i]);
                                aNewCellHit.AddEnergyDeposition(energyDepsCytoplasmSim[i], volumeTypesCytoplasmSim[i]);
                                storedCellHits.push_back(aNewCellHit);
                            }
                        }
                    }
                }

            }

            // Break loop if number of decays have been reached
            if(numberDecays212PbInCytoplasm24h_counter >= numberDecays212PbInCytoplasmTotalTime)
            {
                break;
            }
        }


        //------------------–----------
        // Looping over all stored cell energy depositions
        for(int i=0; i<storedCellHits.size(); i++)
        {
            // Adding energy cell energy depositions to histograms
            energyDepHistograms.AddCellHitsToHistograms(storedCellHits[i]);
        }
    };




    //------------------–----------
    // Filling histgrams
    for(int i=0; i<numberIterations; i++)
    {
        FillHistograms();
    }


    //------------------–----------
    // Scaling histograms
    double scalingFactor = (numberCells*volumeRatio)*((double) numberIterations);
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
    int numberCells = 1000000;

    std::cout << "Number cells: " << numberCells*volumeRatio << std::endl;

    //------------------–----------
    // Defining decay dynamics
    DecayDynamics decays_A5kBq_C4_2 = DecayDynamics(5,1.14,1.16,"C4-2");
    DecayDynamics decays_A10kBq_C4_2 = DecayDynamics(10,1.97,1.98,"C4-2");
    DecayDynamics decays_A25kBq_C4_2 = DecayDynamics(25,4.78,5.91,"C4-2");
    DecayDynamics decays_A50kBq_C4_2 = DecayDynamics(50,8.94,11.07,"C4-2");
    DecayDynamics decays_A75kBq_C4_2 = DecayDynamics(75,10.79,13.12,"C4-2");
    DecayDynamics decays_A100kBq_C4_2 = DecayDynamics(100,13.16,22.72,"C4-2");
    DecayDynamics decays_A150kBq_C4_2 = DecayDynamics(150,16.40,23.56,"C4-2");

    //------------------–----------
    // Loading decay dynamics calculations
    decays_A5kBq_C4_2.LoadDataFromMathematicaCalculations("../Mathematica/Output");
    decays_A10kBq_C4_2.LoadDataFromMathematicaCalculations("../Mathematica/Output");
    decays_A25kBq_C4_2.LoadDataFromMathematicaCalculations("../Mathematica/Output");
    decays_A50kBq_C4_2.LoadDataFromMathematicaCalculations("../Mathematica/Output");
    decays_A75kBq_C4_2.LoadDataFromMathematicaCalculations("../Mathematica/Output");
    decays_A100kBq_C4_2.LoadDataFromMathematicaCalculations("../Mathematica/Output");
    decays_A150kBq_C4_2.LoadDataFromMathematicaCalculations("../Mathematica/Output");


    std::cout << "150kBq - For 10 iterations need number of events :" << std::endl;

    /*
    solution : Submitted batch job 9239377
    membrane : Submitted batch job 9239390
    cytoplasm : Submitted batch job 9239395
    */


    /* Submitted batch job 9373971, 6 runs

    solution 10k : 41 s
    solution 20k : 48 s
    membrane 10k : 46 s
    membrane 20k : 21 s
    cytoplasm 10k : 21 s
    cytoplasm 20k : 21 s
    */

    // std::cout << "Solution: " << 10.0*decays_A150kBq_C4_2.GetNumberDecaysInSolutionFirstHour()/0.98 <<  " Time using 50 jobs : " << 10.0*(decays_A150kBq_C4_2.GetNumberDecaysInSolutionFirstHour()/0.98)*(60.0/10000.0)*(1.0/3600.)*(1.0/50.0) << " hours" <<  std::endl;
    std::cout << "Solution: 20 000" <<  " Time using 1 job : " << 20000.0*(60.0/10000.0) << " seconds" <<  std::endl;
    // = 120 seconds
    // From slurm : 38.0 s

    // std::cout << "Membrane: " << 10.0*decays_A150kBq_C4_2.GetNumberDecaysInMembraneTotalTime() << " Time using 50 jobs :" << 10.0*decays_A150kBq_C4_2.GetNumberDecaysInMembraneTotalTime()*(52.0/10000.0)*(1.0/3600.)*(1.0/50.0) << " hours"  << std::endl;
    std::cout << "Membrane: 20 000" << " Time using 1 job : " << 20000.0*(52.0/10000.0) << " seconds" <<  std::endl;
    // = 104 seconds
    // From slurm 17.0 s

    // std::cout << "Cytoplasm: " << (10.0*decays_A150kBq_C4_2.GetNumberDecaysInCytoplasmTotalTime()/0.97) << " Time using 50 jobs :" << (10.0*decays_A150kBq_C4_2.GetNumberDecaysInCytoplasmTotalTime()/0.97)*(45.0/10000.0)*(1.0/3600.)*(1.0/50.0) << " hours"  << std::endl;

    std::cout << "Cytoplasm: 20 000" <<  " Time using 1 job : " << 20000.0*(45.0/10000.0) << " seconds" <<  std::endl;
    // = 90 seconds
    // From slurm 19.0 s

    // std::cout << "In total: " << 10.0*decays_A150kBq_C4_2.GetNumberDecaysInSolutionFirstHour() + 0.0*decays_A150kBq_C4_2.GetNumberDecaysInMembraneTotalTime() + 10.0*decays_A150kBq_C4_2.GetNumberDecaysInCytoplasmTotalTime() << std::endl;

    // std::cout << "Volume Ratio " << volumeRatio << std::endl;
    // std::cout << "Decays in 0.2 mL" << decays_A10_C4_2.GetNumberDecaysInSolutionFirstHour() << std::endl;
    // std::cout << "Decays in cell tube " << decays_A10_C4_2.GetNumberDecaysInSolutionFirstHour()*volumeRatio << std::endl;


    /*
    Solution: 1e+06 Time using 50 jobs : 0.0333333 hours
    Membrane: 1e+06 Time using 50 jobs : 0.0288889 hours
    Membrane: 1e+06 Time using 50 jobs : 0.025 hours
    */

    /*

    int numberIterations = 1;


    // ------------------–----------
    // Creating "average energy deposition hisograms"
    EnergyDepositionHistograms Hist_A5_C4_2 = MakeHistograms(decays_A5_C4_2, numberIterations, volumeRatio, numberCells);
    std::cout << "Activity 5kBq finished" << std::endl;
    EnergyDepositionHistograms Hist_A10_C4_2 = MakeHistograms(decays_A10_C4_2, numberIterations, volumeRatio, numberCells);
    std::cout << "Activity 10kBq finished" << std::endl;
    EnergyDepositionHistograms Hist_A25_C4_2 = MakeHistograms(decays_A25_C4_2, numberIterations, volumeRatio, numberCells);
    std::cout << "Activity 25kBq finished" << std::endl;
    EnergyDepositionHistograms Hist_A50_C4_2 = MakeHistograms(decays_A50_C4_2, numberIterations, volumeRatio, numberCells);
    std::cout << "Activity 50kBq finished" << std::endl;
    EnergyDepositionHistograms Hist_A75_C4_2 = MakeHistograms(decays_A75_C4_2, numberIterations, volumeRatio, numberCells);
    std::cout << "Activity 75kBq finished" << std::endl;
    EnergyDepositionHistograms Hist_A100_C4_2 = MakeHistograms(decays_A100_C4_2, numberIterations, volumeRatio, numberCells);
    std::cout << "Activity 100kBq finished" << std::endl;
    EnergyDepositionHistograms Hist_A150_C4_2 = MakeHistograms(decays_A150_C4_2, numberIterations, volumeRatio, numberCells);
    std::cout << "Activity 150kBq finished" << std::endl;



    // double integral = Hist_A150_C4_2.GetEnergyDepNucleusHist()->Integral();
    // std::cout << "integral:" << integral << std::endl;

    //------------------–----------
    // Writing histograms to file
    Hist_A150_C4_2.WriteHistogramsToFile();
    // Hist_A10_C4_2.WriteHistogramsToFile();
    // Hist_A25_C4_2.WriteHistogramsToFile();
    // Hist_A50_C4_2.WriteHistogramsToFile();
    // Hist_A75_C4_2.WriteHistogramsToFile();
    // Hist_A100_C4_2.WriteHistogramsToFile();
    // Hist_A150_C4_2.WriteHistogramsToFile();

    //------------------–----------
    outputMainAnalysis->Write();
    outputMainAnalysis->Close();




    // //------------------–----------
    // // Making outputfile
    // EnergyDepositionHistograms Hist_A150kBq_C4_2 = MakeHistograms(decays_A150kBq_C4_2, numberIterations, volumeRatio, numberCells);
    // auto outputMainAnalysis_150kBq = new TFile("outputMainAnalysisCode_150kBq.root", "RECREATE");
    // Hist_A150kBq_C4_2.WriteHistogramsToFile();
    // outputMainAnalysis_150kBq->Write();
    // outputMainAnalysis_150kBq->Close();

    */

}

