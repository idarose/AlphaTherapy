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
std::vector<double> importDataStringDouble(std::string filename)
{
    //----------------------
    //  Imports data from file where the first column contains data in string form
    //  and second column contains data in double form. Returns the columns of double
    //  values in the form of a vector
    std::vector<double> input_data;

    std::fstream myfile;
    myfile.open(filename);

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
class decayDynamics
{
    //----------------------
    //  Class to store calculated decay dynamics data from the Mathematica calculations.
    //  Assumes the output file from Mathematica is structured in a specific way

    public:
        decayDynamics(int activitySample_in, double U0InternalizedPerCell_in, double U0SurfaceBoundPerCell_in, std::string cellLine_in);

        void loadDataFromMathematicaCalculations(std::string filepathToMathematicaOutput);

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


decayDynamics::decayDynamics(int activitySample_in, double U0InternalizedPerCell_in, double U0SurfaceBoundPerCell_in, std::string cellLine_in)
{
    activitySample = activitySample_in;
    U0InternalizedPerCell = U0InternalizedPerCell_in;
    U0SurfaceBoundPerCell = U0SurfaceBoundPerCell_in;
    cellLine = cellLine_in;

}



void decayDynamics::loadDataFromMathematicaCalculations(std::string filepathToMathematicaOutput)
{
    //----------------------
    // Loads data from calculations, assuming output files are structured as shown in the file "212PbDecayDynamics.nb"

    std::string filepathSolutionData = filepathToMathematicaOutput + "/" + cellLine + "/Solution/Activity_" + std::to_string(activitySample) + "kBq/Decays.dat";

    // Importing data for decays occuring in solution
    std::vector<double> decayDataSolution = importDataStringDouble(filepathSolutionData);
    numberDecays212PbInSolutionFirstHour = decayDataSolution[0];

    // Importing data for decays occuring in cells
    std::string filepathCellData = filepathToMathematicaOutput + "/" + cellLine + "/Cells/Activity_" + std::to_string(activitySample) + "kBq/Decays.dat";
    std::vector<double> decayDataCells = importDataStringDouble(filepathCellData);

    numberDecays212PbInMembraneTotalTime = decayDataCells[7];
    numberDecays212PbInCytoplasmTotalTime = decayDataCells[14];

}



//------------------–----------
class cellHit
{
    //----------------------
    //  Class to store information on the energy depositions in ONE specific cell.
    //  When a new cell is "registered" as hit by a particle, a new instance of this class is created.
    //  Any further energy depositions registered in this cell then gets added to the previously
    //  registered energy deposition. The energy deposition is stored according to which "part" of
    //  the cell is hit by the.

    public:
        cellHit(int cellID_in);

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
cellHit::cellHit(int cellID_in)
{
    cellID = cellID_in;
    energyDepMembrane = 0.0;
    energyDepCytoplasm = 0.0;
    energyDepNucleus = 0.0;
}


//------------------–----------
void cellHit::AddEnergyDeposition(double energyDep_in, double volumeType_in)
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
class energyDepositionHistograms
{
    public:
        energyDepositionHistograms(int NBins_in, double EMin_in, double EMax_in);

        void GenerateEmptyHistograms(decayDynamics decayDynamicsInstance);
        void AddCellHitsToHistograms(cellHit cellHitPerEvent);
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
energyDepositionHistograms::energyDepositionHistograms(int NBins_in, double EMin_in, double EMax_in)
{
    NBins = NBins_in;
    EMin = EMin_in;
    EMax = EMax_in;
}


//------------------–----------
void energyDepositionHistograms::GenerateEmptyHistograms(decayDynamics decayDynamicsInstance)
{
    std::string generalHistogramName = "hEnergyDeps_212Pb_" + decayDynamicsInstance.GetCellLine() + "_" + std::to_string(decayDynamicsInstance.GetActivity()) + "kBq_";

    std::string histogramNameNucleus = generalHistogramName + "_Nucleus";

    hEnergyDepsNucleus = new TH1D("hEnergyDepsNucleus", "Energy Deposition in Cell Nucleus / Number of Cells", NBins, EMin, EMax);
    hEnergyDepsMembrane = new TH1D("hEnergyDepsMembrane", "Energy Depsition in Cell Membrane / Number of Cells", NBins, EMin, EMax);
    hEnergyDepsCytoplasm = new TH1D("hEnergyDepsCytoplasm", "Energy Depsition in Cell Cytoplasm / Number of Cells", NBins, EMin, EMax);
    hEnergyDepsCellTotal = new TH1D("hEnergyDepsCellTotal", "Energy Depsition in Cell / Number of Cells", NBins, EMin, EMax);
    hEnergyDepsMembraneAndCytoplasm = new TH1D("hEnergyDepsMembraneCytoplasm", "Energy Depsition in Cell Membrane and Cell Cytoplasm / Number of Cells", NBins, EMin, EMax);
    hEnergyDepsMembraneAndNucleus = new TH1D("hEnergyDepsMembraneAndNucleus", "Energy Depsition in Cell Membrane and Cell Nucleus / Number of Cells", NBins, EMin, EMax);
    hEnergyDepsNucleusAndCytoplasm = new TH1D("hEnergyDepsNucleusAndCytoplasm", "Energy Depsition in Cell Nucleus and Cell Cytoplasm / Number of Cells", NBins, EMin, EMax);
}

//------------------–----------
void energyDepositionHistograms::AddCellHitsToHistograms(cellHit cellHitFromEvent)
{
    if(cellHitFromEvent.GetEnergyDepositionMembrane()>0.0)
        {hEnergyDepsMembrane->Fill(cellHitFromEvent.GetEnergyDepositionMembrane());}
    if(cellHitFromEvent.GetEnergyDepositionCytoplasm()>0.0)
        {hEnergyDepsCytoplasm->Fill(cellHitFromEvent.GetEnergyDepositionCytoplasm());}
    if(cellHitFromEvent.GetEnergyDepositionNucleus()>0.0)
        {hEnergyDepsNucleus->Fill(cellHitFromEvent.GetEnergyDepositionNucleus());}
    if(cellHitFromEvent.GetSumEnergyDepositions()>0.0)
        {hEnergyDepsCellTotal->Fill(cellHitFromEvent.GetSumEnergyDepositions());}
    if((cellHitFromEvent.GetEnergyDepositionMembrane()+cellHitFromEvent.GetEnergyDepositionCytoplasm())>0.0)
        {hEnergyDepsMembraneAndCytoplasm->Fill(cellHitFromEvent.GetEnergyDepositionMembrane()+cellHitFromEvent.GetEnergyDepositionCytoplasm());}
    if((cellHitFromEvent.GetEnergyDepositionMembrane()+cellHitFromEvent.GetEnergyDepositionNucleus())>0.0)
    {hEnergyDepsMembraneAndNucleus->Fill(cellHitFromEvent.GetEnergyDepositionMembrane()+cellHitFromEvent.GetEnergyDepositionNucleus());}
    if((cellHitFromEvent.GetEnergyDepositionNucleus()+cellHitFromEvent.GetEnergyDepositionCytoplasm())>0.0){hEnergyDepsNucleusAndCytoplasm->Fill(cellHitFromEvent.GetEnergyDepositionNucleus()+cellHitFromEvent.GetEnergyDepositionCytoplasm());}
}

//------------------–----------
void energyDepositionHistograms::ScaleHistograms(double factor)
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
void energyDepositionHistograms::WriteHistogramsToFile()
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
energyDepositionHistograms makeHistograms(decayDynamics decayDynamicsInstance, int numberIterations, double volumeRatio, int numberCells)
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
    energyDepositionHistograms energyDepHistograms = energyDepositionHistograms(NBins, EMin, EMax);
    energyDepHistograms.GenerateEmptyHistograms(decayDynamicsInstance);



    //------------------–----------
    // Opening TTree files and creating TTreeReaders

    // Reader for solution simulation
    std::shared_ptr<TFile> myFileSolutionSim(TFile::Open("../GEANT4Simulations/B4aSolution-build/B4.root", "READ"));
    auto treeSolutionSim = myFileSolutionSim->Get<TTree>("B4");
    TTreeReader myReaderSolutionSim(treeSolutionSim);

    // Reader for membrane simulation
    std::shared_ptr<TFile> myFileMembraneSim(TFile::Open("../GEANT4Simulations/B4aMembrane-build/B4.root", "READ"));
    auto treeMembraneSim = myFileMembraneSim->Get<TTree>("B4");
    TTreeReader myReaderMembraneSim(treeMembraneSim);

        // Reader for cytoplasm simulation
    std::shared_ptr<TFile> myFileCytoplasmSim(TFile::Open("../GEANT4Simulations/B4aCytoplasm-build/B4.root", "READ"));
    auto treeCytoplasmSim = myFileCytoplasmSim->Get<TTree>("B4");
    TTreeReader myReaderCytoplasmSim(treeCytoplasmSim);

    // std::cout << numberDecays212PbInSolutionFirstHour << std::endl;
    // std::cout << myReaderSolutionSim.GetEntries() << std::endl;
    // std::cout << myReaderSolutionSim.GetEntries()/numberDecays212PbInSolutionFirstHour << std::endl;

    //------------------–----------
    // Accessing brances of tree

    // TTreeReaderArray<int> initialVolumeTypeID(myReaderSolutionSim, "initialVolumeTypeID");

    // Solution
    TTreeReaderArray<double> energyDepsSolutionSim(myReaderSolutionSim, "EnergyDeps");
    TTreeReaderArray<int> volumeTypesSolutionSim(myReaderSolutionSim, "VolumeTypes");
    TTreeReaderArray<int> cellIDsSolutionSim(myReaderSolutionSim, "CellIDs");
    TTreeReaderArray<double> kineticEnergySolutionSim(myReaderSolutionSim, "KineticEnergy");
    TTreeReaderArray<int> particleTypeSolutionSim(myReaderSolutionSim, "ParticleType");
    TTreeReaderArray<double> interactionTimeSolutionSim(myReaderSolutionSim, "InteractionTime");

    // Membrane
    TTreeReaderArray<double> energyDepsMembraneSim(myReaderMembraneSim, "EnergyDeps");
    TTreeReaderArray<int> volumeTypesMembraneSim(myReaderMembraneSim, "VolumeTypes");
    TTreeReaderArray<int> cellIDsMembraneSim(myReaderMembraneSim, "CellIDs");
    TTreeReaderArray<double> kineticEnergyMembraneSim(myReaderMembraneSim, "KineticEnergy");
    TTreeReaderArray<int> particleTypeMembraneSim(myReaderMembraneSim, "ParticleType");
    TTreeReaderArray<double> interactionTimeMembraneSim(myReaderMembraneSim, "InteractionTime");


    // Cytoplasm
    TTreeReaderArray<double> energyDepsCytoplasmSim(myReaderCytoplasmSim, "EnergyDeps");
    TTreeReaderArray<int> volumeTypesCytoplasmSim(myReaderCytoplasmSim, "VolumeTypes");
    TTreeReaderArray<int> cellIDsCytoplasmSim(myReaderCytoplasmSim, "CellIDs");
    TTreeReaderArray<double> kineticEnergyCytoplasmSim(myReaderCytoplasmSim, "KineticEnergy");
    TTreeReaderArray<int> particleTypeCytoplasmSim(myReaderCytoplasmSim, "ParticleType");
    TTreeReaderArray<double> interactionTimeCytoplasmSim(myReaderCytoplasmSim, "InteractionTime");

    double j = myReaderCytoplasmSim.GetEntries();
    std::cout <<"Events run in G4: " << 500000 << " Entries cytoplasm tree : " << myReaderCytoplasmSim.GetEntries() << " Ratio : " << j/500000. << std::endl;

    double g = 500000.;
    double e = myReaderSolutionSim.GetEntries();
    std::cout << "Events run in G4: " << 500000 <<  " Entries in colution tree: " << myReaderSolutionSim.GetEntries() << " Ratio: " << e/g << std::endl;
    int timesThroughLoop = 0;

    std::cout << "Decays in solution " << numberDecays212PbInSolutionFirstHour << std::endl;

    //------------------–----------
    //  Function for filling histograms
    auto fillHistograms = [&]()
    {

        int numberDecays212PbInSolutionFirstHour_counter = 0;

        //------------------–----------
        // Looping through data for decays occuring in solution in first hour
        while(myReaderSolutionSim.Next())
        {

            // std::cout << numberDecays212PbInSolutionFirstHour_counter <<std::endl;
            // Vector to store all cell hits for one event/decay
            std::vector<cellHit> storedInfoForEvent;


            // If first interaction took place in first hour update counter
            if(energyDepsSolutionSim.GetSize()>0)
            {
                if(interactionTimeSolutionSim[0]/3600.0 < 1.0)
                {
                    numberDecays212PbInSolutionFirstHour_counter++;
                }
            }

            // looping over all steps for one event/decay
            for(int i=0; i<energyDepsSolutionSim.GetSize(); i++)
            {

                // Only add if actual energy deposition
                if(energyDepsSolutionSim[i]>0.)
                {

                    // Only add if interaction happened within first hour
                    if(interactionTimeSolutionSim[i]/3600.0 < 1.0)
                    {

                        // Register first cell hit
                        if(storedInfoForEvent.size()==0)
                        {
                            cellHit aNewCellHit = cellHit(cellIDsSolutionSim[i]);
                            aNewCellHit.AddEnergyDeposition(energyDepsSolutionSim[i],volumeTypesSolutionSim[i]);
                            storedInfoForEvent.push_back(aNewCellHit);
                        }

                        // Register remaining cell hits
                        else
                        {
                            // boolean: true if cell has not already been hit before in this event
                            // false if already hit this event
                            bool cellNotAlreadyHit = true;

                            // Loop over all previous cell hits stored
                            for(int ii=0; ii<storedInfoForEvent.size(); ii++)
                            {
                                // Add energy deposition to cell that has already been hit before
                                if(cellIDsSolutionSim[i]==storedInfoForEvent[ii].GetCellID())
                                {
                                    storedInfoForEvent[ii].AddEnergyDeposition(energyDepsSolutionSim[i],volumeTypesSolutionSim[i]);

                                    // Update boolean
                                    cellNotAlreadyHit = false;
                                }
                            }

                            // Register a new cell hit, and add energy deposition
                            if(cellNotAlreadyHit)
                            {
                                cellHit aNewCellHit = cellHit(cellIDsSolutionSim[i]);
                                aNewCellHit.AddEnergyDeposition(energyDepsSolutionSim[i], volumeTypesSolutionSim[i]);
                                storedInfoForEvent.push_back(aNewCellHit);
                            }
                        }
                    }
                }
            }

            // Looping over all stored cell energy depositions
            for(int ii=0; ii<storedInfoForEvent.size(); ii++)
            {
                // Adding energy cell energy depositions to histograms
                energyDepHistograms.AddCellHitsToHistograms(storedInfoForEvent[ii]);
            }

            // Break loop if number of decays have been reached
            if(numberDecays212PbInSolutionFirstHour_counter >= numberDecays212PbInSolutionFirstHour)
            {
                timesThroughLoop++;
                std::cout << "Round: " << timesThroughLoop << " Current Entry: " << myReaderSolutionSim.GetCurrentEntry() <<  std::endl;
                break;
            }
        }

        //------------------–----------
        // Looping through data for decays occuring in Membrane in 24 hours

        int numberDecays212PbInMembraneTotalTime_counter = 0;
        while(myReaderMembraneSim.Next())
        {
            // Vector to store all cell hits for one event/decay
            std::vector<cellHit> storedInfoForEvent;

            // If first interaction took place in first 24 hours update counter
            if(energyDepsMembraneSim.GetSize()>0)
            {
                if(interactionTimeMembraneSim[0]/3600.0 < 1.0)
                {
                    numberDecays212PbInMembraneTotalTime_counter ++;
                }
            }

            // looping over all steps for one event/decay
            for(int i=0; i<energyDepsMembraneSim.GetSize(); i++)
            {

                // Only add if actual energy deposition
                if(energyDepsMembraneSim[i]>0.)
                {

                    // Only add if interaction happened within first 24 hours
                    if(interactionTimeMembraneSim[i]/3600.0 < 24.0)
                    {

                        // Register first cell hit
                        if(storedInfoForEvent.size()==0)
                        {
                            cellHit aNewCellHit = cellHit(cellIDsMembraneSim[i]);
                            aNewCellHit.AddEnergyDeposition(energyDepsMembraneSim[i],volumeTypesMembraneSim[i]);
                            storedInfoForEvent.push_back(aNewCellHit);

                        }

                        // Register remaining cell hits
                        else
                        {
                            // boolean: true if cell has not already been hit before in this event
                            // false if already hit this event
                            bool cellNotAlreadyHit = true;

                            // Loop over all previous cell hits stored
                            for(int ii=0; ii<storedInfoForEvent.size(); ii++)
                            {

                                // Add energy deposition to cell that has already been hit before
                                if(cellIDsMembraneSim[i]==storedInfoForEvent[ii].GetCellID())
                                {
                                    storedInfoForEvent[ii].AddEnergyDeposition(energyDepsMembraneSim[i],volumeTypesMembraneSim[i]);

                                    // Update boolean
                                    cellNotAlreadyHit = false;
                                }
                            }

                            // Register a new cell hit, and add energy deposition
                            if(cellNotAlreadyHit)
                            {
                                cellHit aNewCellHit = cellHit(cellIDsMembraneSim[i]);
                                aNewCellHit.AddEnergyDeposition(energyDepsMembraneSim[i], volumeTypesMembraneSim[i]);
                                storedInfoForEvent.push_back(aNewCellHit);
                            }
                        }
                    }
                }
            }

            // Looping over all stored cell energy depositions
            for(int ii=0; ii<storedInfoForEvent.size(); ii++)
            {
                // Adding energy cell energy depositions to histograms
                energyDepHistograms.AddCellHitsToHistograms(storedInfoForEvent[ii]);
            }

            // Break loop if number of decays have been reached
            if(numberDecays212PbInMembraneTotalTime_counter >= numberDecays212PbInMembraneTotalTime)
            {
                break;
            }
        }


        //------------------–----------
        // Looping through data for decays occuring in Cytoplasm in 24 hours

        int numberDecays212PbInCytoplasmTotalTime_counter = 0;
        while(myReaderCytoplasmSim.Next())
        {
            // Vector to store all cell hits for one event/decay
            std::vector<cellHit> storedInfoForEvent;

            // If first interaction took place in first 24 hours update counter
            if(energyDepsCytoplasmSim.GetSize()>0)
            {
                if(interactionTimeCytoplasmSim[0]/3600.0 < 1.0)
                {
                    numberDecays212PbInCytoplasmTotalTime_counter ++;
                }
            }

            // looping over all steps for one event/decay
            for(int i=0; i<energyDepsCytoplasmSim.GetSize(); i++)
            {

                // Only add if actual energy deposition
                if(energyDepsCytoplasmSim[i]>0.)
                {

                    // Only add if interaction happened within first 24 hours
                    if(interactionTimeCytoplasmSim[i]/3600.0 < 24.0)
                    {

                        // Register first cell hit
                        if(storedInfoForEvent.size()==0)
                        {
                            cellHit aNewCellHit = cellHit(cellIDsCytoplasmSim[i]);
                            aNewCellHit.AddEnergyDeposition(energyDepsCytoplasmSim[i],volumeTypesCytoplasmSim[i]);
                            storedInfoForEvent.push_back(aNewCellHit);

                        }

                        // Register remaining cell hits
                        else
                        {
                            // boolean: true if cell has not already been hit before in this event
                            // false if already hit this event
                            bool cellNotAlreadyHit = true;

                            // Loop over all previous cell hits stored
                            for(int ii=0; ii<storedInfoForEvent.size(); ii++)
                            {

                                // Add energy deposition to cell that has already been hit before
                                if(cellIDsCytoplasmSim[i]==storedInfoForEvent[ii].GetCellID())
                                {
                                    storedInfoForEvent[ii].AddEnergyDeposition(energyDepsCytoplasmSim[i],volumeTypesCytoplasmSim[i]);

                                    // Update boolean
                                    cellNotAlreadyHit = false;
                                }
                            }

                            // Register a new cell hit, and add energy deposition
                            if(cellNotAlreadyHit)
                            {
                                cellHit aNewCellHit = cellHit(cellIDsCytoplasmSim[i]);
                                aNewCellHit.AddEnergyDeposition(energyDepsCytoplasmSim[i], volumeTypesCytoplasmSim[i]);
                                storedInfoForEvent.push_back(aNewCellHit);
                            }
                        }
                    }
                }
            }

            // Looping over all stored cell energy depositions
            for(int ii=0; ii<storedInfoForEvent.size(); ii++)
            {
                // Adding energy cell energy depositions to histograms
                energyDepHistograms.AddCellHitsToHistograms(storedInfoForEvent[ii]);
            }

            // Break loop if number of decays have been reached
            if(numberDecays212PbInCytoplasmTotalTime_counter >= numberDecays212PbInCytoplasmTotalTime)
            {
                break;
            }
        }
    };

    for(int i=0; i<numberIterations; i++)
    {
        fillHistograms();
    }

    double a = myReaderSolutionSim.GetCurrentEntry();
    double b = timesThroughLoop;
    std::cout << "Entries needed per loop: " << a/b << std::endl;
    std::cout << "For " << timesThroughLoop << " iterations need " << 100*a/b << " entries = " << (100*a/b)/(e/g) << " G4 events" << std::endl;

    double scalingFactor = (numberCells*volumeRatio + 0.5)*numberIterations;
    energyDepHistograms.ScaleHistograms(1/scalingFactor);

    // outputTest->Write();
    // outputTest->Close();

    return energyDepHistograms;

}

void mainAnalysisCode()
{
    // Calculating volume ratio
    double VolumeSample = 0.2*1000; // mm^3
    double volumeCellTube = TMath::Pi()*std::pow(0.5,2.0)*1.0; // mm^3
    double volumeRatio = volumeCellTube/VolumeSample;
    int numberCells = 1000000;

    //------------------–----------
    // Defining decay dynamics
    decayDynamics decays_A10_C4_2 = decayDynamics(10,1.97,1.98,"C4-2");
    decayDynamics decays_A25_C4_2 = decayDynamics(25,4.78,5.91,"C4-2");

    //------------------–----------
    // Loading decay dynamics calculations
    decays_A10_C4_2.loadDataFromMathematicaCalculations("../Mathematica/Output");
    decays_A25_C4_2.loadDataFromMathematicaCalculations("../Mathematica/Output");

    // std::cout << "Volume Ratio " << volumeRatio << std::endl;
    // std::cout << "Decays in 0.2 mL" << decays_A10_C4_2.GetNumberDecaysInSolutionFirstHour() << std::endl;
    // std::cout << "Decays in cell tube " << decays_A10_C4_2.GetNumberDecaysInSolutionFirstHour()*volumeRatio << std::endl;


    int numberIterations = 10;

    //------------------–----------
    // Creating "average energy deposition hisograms"
    energyDepositionHistograms Hist_A10_C4_2 = makeHistograms(decays_A10_C4_2, numberIterations, volumeRatio, numberCells);


    //------------------–----------
    // Making outputfile
    auto outputMainAnalysis = new TFile("outputMainAnalysisCode.root", "RECREATE");


    //------------------–----------
    // Writing histograms to file
    Hist_A10_C4_2.WriteHistogramsToFile();


    //------------------–----------
    outputMainAnalysis->Write();
    outputMainAnalysis->Close();

}

