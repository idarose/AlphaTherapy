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
    public:
        decayDynamics(int activitySample_in, int nuclidesInternalizedCell_in, int nuclidesInternalizedCytoplasm_in, std::string cellLine_in);

        void loadDataFromMathematicaCalculations(std::string filepathToMathematicaOutput);

        int GetNumberDecaysSolutionFirstHour(){return numberDecays212PbSolutionFirstHour;};
        int GetNumberDecaysMembraneTotalTime(){return numberDecays212PbMembraneTotalTime;};
        int GetNumberDecaysCytoplasmTotalTime(){return numberDecays212PbCytoplasmTotalTime;};
        // int GetActivtySample(){return activity;};

    private:
        int numberDecays212PbSolutionFirstHour;
        int numberDecays212PbMembraneTotalTime;
        int numberDecays212PbCytoplasmTotalTime;

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

/*
The way the Mathematica output should be structure:


Output/CellLine/NumberRadionuclidesMembrane-NumberRadionuclidesCytoplasm/

and then either
/Cells/ or /Solution/

Then files of the form

/Decays_Activity#_Cells(Solution).dat
*/
void decayDynamics::loadDataFromMathematicaCalculations(std::string filepathToMathematicaOutput)
{
    std::string filepathSolutionData = filepathToMathematicaOutput + "/" + cellLine + std::to_string(nuclidesInternalizedCell) + "-" + std::to_string(nuclidesInternalizedCytoplasm) + "/Solution/Decays_Activity10_Solution.dat";

    std::vector<double> decayDataSolution = importDataStringDouble(filepathSolutionData);
    numberDecays212PbSolutionFirstHour = decayDataSolution[0];

    std::string filepathCellData = filepathToMathematicaOutput + "/" + cellLine + std::to_string(nuclidesInternalizedCell) + "-" + std::to_string(nuclidesInternalizedCytoplasm) + "/Cells/Decays_Activity10_Cells.dat";
    std::vector<double> decayDataCells = importDataStringDouble(filepathCellData);

    numberDecays212PbMembraneTotalTime = decayDataCells[7];
    numberDecays212PbCytoplasmTotalTime = decayDataCells[14];

}



//------------------–----------
class cellHit
{
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


void mainAnalysisCode()
{
    //------------------–----------
    // Calculating number of cells
    double cellsSample = 1000000;
    double volumeCellSample = 0.2*1000.0; // mm^3
    double volumeCellTube = TMath::Pi()*std::pow(0.5,2.0)*1.0; //mm^3

    int numberCells = cellsSample*volumeCellTube/volumeCellSample;

    // std::cout << "Ratio volumes sample(0.2mL)/tube = " << volumeCellTube/volumeCellSample << std::endl;

    //------------------–----------
    // Importing information from Mathematica calculations

    decayDynamics Activity10_C4_2 = decayDynamics(10,4,2,"C");
    int decays212PbSolutionFirstHourIn2mLSample = Activity10_C4_2.GetNumberDecaysSolutionFirstHour();
    int decays212PbSolutionFirstHourCellTube = decays212PbSolutionFirstHourIn2mLSample*volumeCellTube/volumeCellSample;

    //------------------–----------
    // Opening TTree file and creating TTreeReader
    std::unique_ptr<TFile> myFile(TFile::Open("../GEANT4Simulations/B4aSolution-build/B4.root", "READ"));
    auto tree = myFile->Get<TTree>("B4");
    TTreeReader myReader(tree);


    //------------------–----------
    // Accessing brances of tree
    TTreeReaderArray<double> energyDeps(myReader, "EnergyDeps");
    TTreeReaderArray<int> volumeTypes(myReader, "VolumeTypes");
    TTreeReaderArray<int> cellIDs(myReader, "CellIDs");
    TTreeReaderArray<double> kineticEnergy(myReader, "KineticEnergy");
    TTreeReaderArray<int> particleType(myReader, "ParticleType");
    TTreeReaderArray<double> interactionTime(myReader, "InteractionTime");


    //------------------–----------
    int NBins = 2000000;
    double EMin = 0.0;
    double EMax = 20.0;
    // Histogram for total energy deposited in one nuclei per decay
    TH1D *hEnergyDepsNucleus = new TH1D("hEnergyDepsNucleus", "Energy Deposition in Cell Nucleus / Decay", NBins, EMin, EMax);

    // Histogram for total energy deposited in one membrane per decay
    TH1D *hEnergyDepsMembrane = new TH1D("hEnergyDepsMembrane", "Energy Depsition in Cell Membrane / Decay", NBins, EMin, EMax);

    // Histogram for total energy deposited in one cytoplasm per decay
    TH1D *hEnergyDepsCytoplasm = new TH1D("hEnergyDepsCytoplasm", "Energy Depsition in Cell Cytoplasm / Decay", NBins, EMin, EMax);

    // Histogram for total energy deposited in one cytoplasm per decay
    TH1D *hEnergyDepsCellTotal = new TH1D("hEnergyDepsCellTotal", "Energy Depsition in Cell / Decay", NBins, EMin, EMax);



    //------------------–----------
    // Making outputfile
    auto OutputTuplesAnalysis = new TFile("outputMainAnalysisCode.root", "RECREATE");



    //------------------–----------
    // DECAYS IN SOLUTION


    // Counter to break loop when number of decays have been reached
    int numberDecaysSolution_counter = 0;

    while(myReader.Next())
    {
        std::vector<cellHit> storedInfoForEvent;
        for(int i=0; i<energyDeps.GetSize(); i++)
        {
            if(energyDeps[i]!=0.)
            {
                if(interactionTime[i]/3600.0 < 1.0)
                {
                    if(storedInfoForEvent.size()==0)
                    {
                        numberDecaysSolution_counter ++;

                        cellHit aNewCellHit = cellHit(cellIDs[i]);
                        aNewCellHit.AddEnergyDeposition(energyDeps[i],volumeTypes[i]);
                        storedInfoForEvent.push_back(aNewCellHit);
                    }
                    else
                    {
                        for(int ii=0; ii<storedInfoForEvent.size(); ii++)
                        {
                            if(cellIDs[i]==storedInfoForEvent[ii].GetCellID())
                            {
                                storedInfoForEvent[ii].AddEnergyDeposition(energyDeps[i],volumeTypes[i]);
                            }
                            else
                            {
                                cellHit aNewCellHit = cellHit(cellIDs[i]);
                                aNewCellHit.AddEnergyDeposition(energyDeps[i], volumeTypes[i]);
                                storedInfoForEvent.push_back(aNewCellHit);
                            }
                        }
                    }
                }
            }
        }
        for(int ii=0; ii<storedInfoForEvent.size(); ii++)
        {
            hEnergyDepsMembrane->Fill(storedInfoForEvent[ii].GetEnergyDepositionMembrane());
            hEnergyDepsCytoplasm->Fill(storedInfoForEvent[ii].GetEnergyDepositionCytoplasm());
            hEnergyDepsNucleus->Fill(storedInfoForEvent[ii].GetEnergyDepositionNucleus());
            hEnergyDepsCellTotal->Fill(storedInfoForEvent[ii].GetSumEnergyDepositions());
        }
        if(numberDecaysSolution_counter>=decays212PbSolutionFirstHourCellTube)
        {
            break;
        }

    }



    //------------------–----------
    // Scaling histograms
    hEnergyDepsMembrane->Scale(1.0/numberCells);
    hEnergyDepsCytoplasm->Scale(1.0/numberCells);
    hEnergyDepsNucleus->Scale(1.0/numberCells);
    hEnergyDepsCellTotal->Scale(1.0/numberCells);


    //------------------–----------
    hEnergyDepsMembrane->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    hEnergyDepsMembrane->GetYaxis()->SetTitle("Hits / Cell");
    hEnergyDepsCytoplasm->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    hEnergyDepsCytoplasm->GetYaxis()->SetTitle("Hits / Cell");
    hEnergyDepsNucleus->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    hEnergyDepsNucleus->GetYaxis()->SetTitle("Hits / Cell");
    hEnergyDepsCellTotal->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    hEnergyDepsCellTotal->GetYaxis()->SetTitle("Hits / Cell");


    //------------------–----------
    hEnergyDepsMembrane->Write();
    hEnergyDepsCytoplasm->Write();
    hEnergyDepsNucleus->Write();
    hEnergyDepsCellTotal->Write();


    //------------------–----------
    OutputTuplesAnalysis->Write();
    OutputTuplesAnalysis->Close();


}

