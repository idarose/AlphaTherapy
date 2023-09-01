#include <TFile.h>
#include <TH1.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

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
            std::stringstream mysstream(line);
            mysstream >> x >> y;
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
        decayDynamics(int activitySample_in, int nuclidesInternalizedCell_in, int nuclidesInternalizedCytoplasm_in, std::string cellLine_in);

        void loadDataFromMathematicaCalculations(std::string filepathToMathematicaOutput);

        int GetNumberDecaysInSolutionFirstHour(){return numberDecays212PbInSolutionFirstHour;};
        int GetNumberDecaysInMembraneTotalTime(){return numberDecays212PbInMembraneTotalTime;};
        int GetNumberDecaysInCytoplasmTotalTime(){return numberDecays212PbInCytoplasmTotalTime;};
        // int GetActivtySample(){return activity;};

    private:
        double volumeCellSample = 0.2*1000.0; // mm^3
        double volumeCellTube = TMath::Pi()*std::pow(0.5,2.0)*1.0; //mm^3
        double volumeRatio = volumeCellTube/volumeCellSample;

        int numberDecays212PbInSolutionFirstHour;
        int numberDecays212PbInMembraneTotalTime;
        int numberDecays212PbInCytoplasmTotalTime;

        int nuclidesInternalizedCell;
        int nuclidesInternalizedCytoplasm;
        int activitySample;  // Bq/1mL
        std::string cellLine;
};


decayDynamics::decayDynamics(int activitySample_in, int nuclidesInternalizedCell_in, int nuclidesInternalizedCytoplasm_in, std::string cellLine_in)
{
    activitySample = activitySample_in;
    nuclidesInternalizedCell = nuclidesInternalizedCell_in;
    nuclidesInternalizedCytoplasm = nuclidesInternalizedCytoplasm_in;
    cellLine = cellLine_in;
}



void decayDynamics::loadDataFromMathematicaCalculations(std::string filepathToMathematicaOutput)
{
    //----------------------
    // Loads data from calculations, assuming output files are structured as in the file "212PbDecayDynamics.nb"

    std::string filepathSolutionData = filepathToMathematicaOutput + "/" + cellLine + std::to_string(nuclidesInternalizedCell) + "-" + std::to_string(nuclidesInternalizedCytoplasm) + "/Solution/Decays_Activity10_Solution.dat";

    std::vector<double> decayDataSolution = importDataStringDouble(filepathSolutionData);
    numberDecays212PbInSolutionFirstHour = decayDataSolution[0]*volumeRatio;

    std::string filepathCellData = filepathToMathematicaOutput + "/" + cellLine + std::to_string(nuclidesInternalizedCell) + "-" + std::to_string(nuclidesInternalizedCytoplasm) + "/Cells/Decays_Activity10_Cells.dat";
    std::vector<double> decayDataCells = importDataStringDouble(filepathCellData);

    numberDecays212PbInMembraneTotalTime = decayDataCells[7]*volumeRatio;
    numberDecays212PbInCytoplasmTotalTime = decayDataCells[14]*volumeRatio;

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

        void AddCellHitsHistograms(cellHit cellHitPerEvent);
        void ScaleHistograms(int factor);
        void generateHistograms();

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

void energyDepositionHistograms::generateHistograms()
{
    hEnergyDepsNucleus = new TH1D("hEnergyDepsNucleus", "Energy Deposition in Cell Nucleus / Number of Cells", NBins, EMin, EMax);
    hEnergyDepsMembrane = new TH1D("hEnergyDepsMembrane", "Energy Depsition in Cell Membrane / Number of Cells", NBins, EMin, EMax);
    hEnergyDepsCytoplasm = new TH1D("hEnergyDepsCytoplasm", "Energy Depsition in Cell Cytoplasm / Number of Cells", NBins, EMin, EMax);
    hEnergyDepsCellTotal = new TH1D("hEnergyDepsCellTotal", "Energy Depsition in Cell / Number of Cells", NBins, EMin, EMax);
    hEnergyDepsMembraneAndCytoplasm = new TH1D("hEnergyDepsMembraneCytoplasm", "Energy Depsition in Cell Membrane and Cell Cytoplasm / Number of Cells", NBins, EMin, EMax);
    hEnergyDepsMembraneAndNucleus = new TH1D("hEnergyDepsMembraneAndNucleus", "Energy Depsition in Cell Membrane and Cell Nucleus / Number of Cells", NBins, EMin, EMax);
    hEnergyDepsNucleusAndCytoplasm = new TH1D("hEnergyDepsNucleusAndCytoplasm", "Energy Depsition in Cell Nucleus and Cell Cytoplasm / Number of Cells", NBins, EMin, EMax);
}

//------------------–----------
void energyDepositionHistograms::AddCellHitsHistograms(cellHit cellHitFromEvent)
{
    hEnergyDepsMembrane->Fill(cellHitFromEvent.GetEnergyDepositionMembrane());
    hEnergyDepsCytoplasm->Fill(cellHitFromEvent.GetEnergyDepositionCytoplasm());
    hEnergyDepsNucleus->Fill(cellHitFromEvent.GetEnergyDepositionNucleus());
    hEnergyDepsCellTotal->Fill(cellHitFromEvent.GetSumEnergyDepositions());
    hEnergyDepsMembraneAndCytoplasm->Fill(cellHitFromEvent.GetEnergyDepositionMembrane()+cellHitFromEvent.GetEnergyDepositionCytoplasm());
    hEnergyDepsMembraneAndNucleus->Fill(cellHitFromEvent.GetEnergyDepositionMembrane()+cellHitFromEvent.GetEnergyDepositionNucleus());
    hEnergyDepsNucleusAndCytoplasm->Fill(cellHitFromEvent.GetEnergyDepositionNucleus()+cellHitFromEvent.GetEnergyDepositionCytoplasm());
}

//------------------–----------
void energyDepositionHistograms::ScaleHistograms(int factor)
{
    hEnergyDepsMembrane->Scale(1/factor);
    hEnergyDepsCytoplasm->Scale(1/factor);
    hEnergyDepsNucleus->Scale(1/factor);
    hEnergyDepsCellTotal->Scale(1/factor);
    hEnergyDepsMembraneAndCytoplasm->Scale(1/factor);
    hEnergyDepsMembraneAndNucleus->Scale(1/factor);
    hEnergyDepsNucleusAndCytoplasm->Scale(1/factor);
}


//------------------–----------
energyDepositionHistograms makeHistograms(decayDynamics decayDynamicsInstance, int numberIterations, int numberCells)
{
    //------------------–----------
    // Loading decay dynamics
    int numberDecays212PbInSolutionFirstHour = decayDynamicsInstance.GetNumberDecaysInSolutionFirstHour();
    int numberDecays212PbInMembraneTotalTime = decayDynamicsInstance.GetNumberDecaysInMembraneTotalTime();
    int numberDecays212PbInCytoplasmTotalTime = decayDynamicsInstance.GetNumberDecaysInCytoplasmTotalTime();

    //------------------–----------
    // Making energy deposition histogram instance
    energyDepositionHistograms energyDepHistograms = energyDepositionHistograms(2000000, 0.0, 20.0);

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

    //------------------–----------
    // Accessing brances of tree
    TTreeReaderArray<double> energyDepsSolutionSim(myReaderSolutionSim, "EnergyDeps");
    TTreeReaderArray<int> volumeTypesSolutionSim(myReaderSolutionSim, "VolumeTypes");
    TTreeReaderArray<int> cellIDsSolutionSim(myReaderSolutionSim, "CellIDs");
    TTreeReaderArray<double> kineticEnergySolutionSim(myReaderSolutionSim, "KineticEnergy");
    TTreeReaderArray<int> particleTypeSolutionSim(myReaderSolutionSim, "ParticleType");
    TTreeReaderArray<double> interactionTimeSolutionSim(myReaderSolutionSim, "InteractionTime");

    //------------------–----------
    int NBins = 2000000;
    double EMin = 0.0;
    double EMax = 20.0;

    //------------------–----------
    // Lambdas for filling histograms

    int numberDecays212PbInSolutionFirstHour_counter = 0;

    auto fillHistograms = [&]()
    {
        while(myReaderSolutionSim.Next())
        {
            std::vector<cellHit> storedInfoForEvent;
            for(int i=0; i<energyDepsSolutionSim.GetSize(); i++)
            {
                if(energyDepsSolutionSim[i]!=0.)
                {
                    if(interactionTimeSolutionSim[i]/3600.0 < 1.0)
                    {
                        if(storedInfoForEvent.size()==0)
                        {
                            numberDecays212PbInSolutionFirstHour_counter ++;

                            cellHit aNewCellHit = cellHit(cellIDsSolutionSim[i]);
                            aNewCellHit.AddEnergyDeposition(energyDepsSolutionSim[i],volumeTypesSolutionSim[i]);
                            // storedInfoForEvent.push_back(aNewCellHit);
                        }
                        else
                        {
                            for(int ii=0; ii<storedInfoForEvent.size(); ii++)
                            {
                                if(cellIDsSolutionSim[i]==storedInfoForEvent[ii].GetCellID())
                                {
                                    storedInfoForEvent[ii].AddEnergyDeposition(energyDepsSolutionSim[i],volumeTypesSolutionSim[i]);
                                }
                                else
                                {
                                    cellHit aNewCellHit = cellHit(cellIDsSolutionSim[i]);
                                    aNewCellHit.AddEnergyDeposition(energyDepsSolutionSim[i], volumeTypesSolutionSim[i]);
                                    // storedInfoForEvent.push_back(aNewCellHit);
                                }
                            }
                        }
                    }
                }
            }
            for(int ii=0; ii<storedInfoForEvent.size(); ii++)
            {
                energyDepHistograms.AddCellHitsHistograms(storedInfoForEvent[ii]);
            }
            if(numberDecays212PbInSolutionFirstHour_counter>=numberDecays212PbInSolutionFirstHour)
            {
                break;
            }
        }
    };

    for(int i=0; i<numberIterations; i++)
    {
        fillHistograms();
    }

    energyDepHistograms.ScaleHistograms(numberCells*numberIterations);

    return energyDepHistograms;

}

void mainAnalysisCode()
{
    //------------------–----------
    // Calculating number of cells
    double cellsSample = 1000000;
    double volumeCellSample = 0.2*1000.0; // mm^3
    double volumeCellTube = TMath::Pi()*std::pow(0.5,2.0)*1.0; //mm^3
    int numberCells = cellsSample*volumeCellTube/volumeCellSample;


    decayDynamics A10_C4_2 = decayDynamics(10,4,2,"C");
    energyDepositionHistograms Hist_A10_C4_2 = makeHistograms(A10_C4_2, 10, numberCells);

    // // std::cout << "Ratio volumes sample(0.2mL)/tube = " << volumeCellTube/volumeCellSample << std::endl;

    // //------------------–----------
    // // Importing information from Mathematica calculations

    // decayDynamics Activity10_C4_2 = decayDynamics(10,4,2,"C");
    // int decays212PbSolutionFirstHourIn2mLSample = Activity10_C4_2.GetNumberDecaysSolutionFirstHour();
    // // int decays212PbSolutionFirstHourCellTube = decays212PbSolutionFirstHourIn2mLSample*volumeCellTube/volumeCellSample;

    // //------------------–----------
    // // Opening TTree file and creating TTreeReader
    // std::shared_ptr<TFile> myFile(TFile::Open("../GEANT4Simulations/B4aSolution-build/B4.root", "READ"));
    // auto tree = myFile->Get<TTree>("B4");
    // TTreeReader myReader(tree);


    // //------------------–----------
    // // Accessing brances of tree
    // TTreeReaderArray<double> energyDeps(myReader, "EnergyDeps");
    // TTreeReaderArray<int> volumeTypes(myReader, "VolumeTypes");
    // TTreeReaderArray<int> cellIDs(myReader, "CellIDs");
    // TTreeReaderArray<double> kineticEnergy(myReader, "KineticEnergy");
    // TTreeReaderArray<int> particleType(myReader, "ParticleType");
    // TTreeReaderArray<double> interactionTime(myReader, "InteractionTime");


    // //------------------–----------
    // int NBins = 2000000;
    // double EMin = 0.0;
    // double EMax = 20.0;
    // // Histogram for total energy deposited in one nuclei per decay
    // TH1D *hEnergyDepsNucleus = new TH1D("hEnergyDepsNucleus", "Energy Deposition in Cell Nucleus / Decay", NBins, EMin, EMax);

    // // Histogram for total energy deposited in one membrane per decay
    // TH1D *hEnergyDepsMembrane = new TH1D("hEnergyDepsMembrane", "Energy Depsition in Cell Membrane / Decay", NBins, EMin, EMax);

    // // Histogram for total energy deposited in one cytoplasm per decay
    // TH1D *hEnergyDepsCytoplasm = new TH1D("hEnergyDepsCytoplasm", "Energy Depsition in Cell Cytoplasm / Decay", NBins, EMin, EMax);

    // // Histogram for total energy deposited in one cytoplasm per decay
    // TH1D *hEnergyDepsCellTotal = new TH1D("hEnergyDepsCellTotal", "Energy Depsition in Cell / Decay", NBins, EMin, EMax);



    // //------------------–----------
    // // Making outputfile
    // auto OutputTuplesAnalysis = new TFile("outputMainAnalysisCode.root", "RECREATE");


    // int numberDecaysSolution_counter = 0;

    // auto calculateOneHistogram = [&]()
    // {
    //     while(myReader.Next())
    //     {
    //         std::cout << myReader.GetCurrentEntry() << std::endl;
    //         numberDecaysSolution_counter ++;

    //         if(numberDecaysSolution_counter > 5)
    //         {
    //             std::cout << "-----------------" << std::endl;
    //             numberDecaysSolution_counter=0;
    //             break;
    //         }
    //     }
    // };

    // for(int i=0;i<2;i++)
    // {
    //     calculateOneHistogram();
    // }

    //------------------–----------
    // DECAYS IN SOLUTION


    // Counter to break loop when number of decays have been reached
    // int numberDecaysSolution_counter = 0;

    // for(int i=0; i<2; i++)
    // {
    //     while(myReader.Next())
    //     {
    //         std::cout << myReader.GetCurrentEntry() << std::endl;
    //         numberDecaysSolution_counter ++;

    //         if(numberDecaysSolution_counter > 5){
    //             std::cout << "-----------------" << std::endl;
    //             numberDecaysSolution_counter=0;
    //             break;
    //         }
    //     }
    // }
    // while(myReader.Next())
    // {
    //     std::vector<cellHit> storedInfoForEvent;
    //     for(int i=0; i<energyDeps.GetSize(); i++)
    //     {
    //         if(energyDeps[i]!=0.)
    //         {
    //             if(interactionTime[i]/3600.0 < 1.0)
    //             {
    //                 if(storedInfoForEvent.size()==0)
    //                 {
    //                     numberDecaysSolution_counter ++;

    //                     cellHit aNewCellHit = cellHit(cellIDs[i]);
    //                     aNewCellHit.AddEnergyDeposition(energyDeps[i],volumeTypes[i]);
    //                     storedInfoForEvent.push_back(aNewCellHit);
    //                 }
    //                 else
    //                 {
    //                     for(int ii=0; ii<storedInfoForEvent.size(); ii++)
    //                     {
    //                         if(cellIDs[i]==storedInfoForEvent[ii].GetCellID())
    //                         {
    //                             storedInfoForEvent[ii].AddEnergyDeposition(energyDeps[i],volumeTypes[i]);
    //                         }
    //                         else
    //                         {
    //                             cellHit aNewCellHit = cellHit(cellIDs[i]);
    //                             aNewCellHit.AddEnergyDeposition(energyDeps[i], volumeTypes[i]);
    //                             storedInfoForEvent.push_back(aNewCellHit);
    //                         }
    //                     }
    //                 }
    //             }
    //         }
    //     }
    //     for(int ii=0; ii<storedInfoForEvent.size(); ii++)
    //     {
    //         hEnergyDepsMembrane->Fill(storedInfoForEvent[ii].GetEnergyDepositionMembrane());
    //         hEnergyDepsCytoplasm->Fill(storedInfoForEvent[ii].GetEnergyDepositionCytoplasm());
    //         hEnergyDepsNucleus->Fill(storedInfoForEvent[ii].GetEnergyDepositionNucleus());
    //         hEnergyDepsCellTotal->Fill(storedInfoForEvent[ii].GetSumEnergyDepositions());
    //     }
    //     if(numberDecaysSolution_counter>=decays212PbSolutionFirstHourCellTube)
    //     {
    //         break;
    //     }
    // }



    // //------------------–----------
    // // Scaling histograms
    // hEnergyDepsMembrane->Scale(1.0/numberCells);
    // hEnergyDepsCytoplasm->Scale(1.0/numberCells);
    // hEnergyDepsNucleus->Scale(1.0/numberCells);
    // hEnergyDepsCellTotal->Scale(1.0/numberCells);


    // //------------------–----------
    // hEnergyDepsMembrane->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    // hEnergyDepsMembrane->GetYaxis()->SetTitle("Hits / Cell");
    // hEnergyDepsCytoplasm->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    // hEnergyDepsCytoplasm->GetYaxis()->SetTitle("Hits / Cell");
    // hEnergyDepsNucleus->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    // hEnergyDepsNucleus->GetYaxis()->SetTitle("Hits / Cell");
    // hEnergyDepsCellTotal->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    // hEnergyDepsCellTotal->GetYaxis()->SetTitle("Hits / Cell");


    // //------------------–----------
    // hEnergyDepsMembrane->Write();
    // hEnergyDepsCytoplasm->Write();
    // hEnergyDepsNucleus->Write();
    // hEnergyDepsCellTotal->Write();


    // //------------------–----------
    // OutputTuplesAnalysis->Write();
    // OutputTuplesAnalysis->Close();


}

