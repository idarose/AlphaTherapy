#include <TFile.h>
#include <TH1.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

#include <vector>
#include <iostream>
#include <tuple>


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


void tuplesAnalysis()
{
    std::string cellLine = "C4-2";
    std::string activity = "10";


    //------------------–----------
    // Importing information from Mathematica calculations


    // Number of decays of 212Pb in Solution first hour
    std::vector<double> decayInfoSolution = importDataStringDouble("../Mathematica/Output/C4-2/Solution/Decays_Activity10_Solution.dat");
    int decays212PbSolutionFirstHour = decayInfoSolution[0];

    // Number of decays in Membrane in 24h
    std::vector<double> decayInfoCells = importDataStringDouble("../Mathematica/Output/C4-2/Cells/Decays_Activity10_Cells.dat");
    int decays212PbCellMembrane = decayInfoCells[7];

    // Number of decays in cytoplasm in 24h
    int decays212PbCellCytoplasm = decayInfoCells[14];

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
    // Histogram for total energy deposited in one nuclei per decay
    TH1D *hEnergyDepsNucleus = new TH1D("hEnergyDepsNucleus", "Energy Deposition in Cell Nuclei / Decay", 10000, 0.0, 20.0);

    // Histogram for total energy deposited in one membrane per decay
    TH1D *hEnergyDepsMembrane = new TH1D("hEnergyDepsMembrane", "Energy Depsition in Cell Membrane / Decay", 10000, 0.0, 20.0);

    // Histogram for total energy deposited in one cytoplasm per decay
    TH1D *hEnergyDepsCytoplasm = new TH1D("hEnergyDepsCytoplasm", "Energy Depsition in Cell Cytoplasm / Decay", 10000, 0.0, 20.0);



    //------------------–----------
    // Making outputfile
    auto OutputTuplesAnalysis = new TFile("OutputTuplesAnalysis.root", "RECREATE");


    // Counter to break loop when number of decays have been reached
    int numberDecaysCounter = 0;

    //------------------–----------
    // Looping over every branch/event
    while(myReader.Next())
    {
        // Vector to store energy deposition information per event/decay
        std::vector<std::tuple<int,double,double,double>> storedInfoEvent;


        // Interaction time in hours
        double interactionTimeHours;

        //------------------–----------
        // Looping over every leaf/step
        for(int i=0; i<energyDeps.GetSize(); i++)
        {
            interactionTimeHours = interactionTime[i]/3600.;

            // Only store interactions happening in first hour
            if(interactionTimeHours <= 1.0)
            {
                // If first step store the energy
                if(storedInfoEvent.size()==0)
                {
                    // Count decay if first interactiontime is in first hour
                    numberDecaysCounter ++;

                    // Making tuples
                    if(volumeTypes[i]==1)
                    {
                        //Membrane
                        storedInfoEvent.push_back(make_tuple(cellIDs[i],energyDeps[i],0.0,0.0));
                    }
                    else if(volumeTypes[i]==2)
                    {
                        // Cytoplasm
                        storedInfoEvent.push_back(make_tuple(cellIDs[i],0.0,energyDeps[i], 0.0));
                    }
                    if(volumeTypes[i]==3)
                    {
                        // Nucleus
                        storedInfoEvent.push_back(make_tuple(cellIDs[i],0.0,0.0,energyDeps[i]));
                    }
                }
                // If not first step check for same cellIDS
                else
                {
                    // Looping over every tuple stored
                    for(int j=0;j<storedInfoEvent.size();j++)
                    {
                        // If the cellID has already been stored before then add the enery
                        if(get<0>(storedInfoEvent[j])==cellIDs[i])
                        {
                            // If in membrane
                            if(volumeTypes[i]==1)
                            {
                                get<1>(storedInfoEvent[j]) += energyDeps[i];
                            }
                            // If in cytoplasm
                            else if(volumeTypes[i]==2)
                            {
                                get<2>(storedInfoEvent[j]) += energyDeps[i];
                            }
                            // If in nucleus
                            else if(volumeTypes[i]==3)
                            {
                                get<3>(storedInfoEvent[j]) += energyDeps[i];
                            }
                        }
                    }
                }
            }

        }

        //------------------–----------
        // Break loop when number of decays reached
        if(numberDecaysCounter >= decays212PbSolutionFirstHour)
        {
            break;
        }

        //------------------–----------
        // Storing info from one event to histograms
        for(int k=0; k<storedInfoEvent.size(); k++)
        {
            hEnergyDepsMembrane->Fill(get<1>(storedInfoEvent[k]));
            hEnergyDepsCytoplasm->Fill(get<2>(storedInfoEvent[k]));
            hEnergyDepsNucleus->Fill(get<3>(storedInfoEvent[k]));
        }
    }

    // double ratio = 0.3596396;
    // std::cout << "Number of runs needed: " << decays212PbSolutionFirstHour/ratio << std::endl;

    //------------------–----------
    // Scaling histograms
    hEnergyDepsMembrane->Scale(1.0/numberDecaysCounter);
    hEnergyDepsCytoplasm->Scale(1.0/numberDecaysCounter);
    hEnergyDepsNucleus->Scale(1.0/numberDecaysCounter);


    //------------------–----------
    hEnergyDepsMembrane->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    hEnergyDepsMembrane->GetYaxis()->SetTitle("Counts / Num. Decays");
    hEnergyDepsCytoplasm->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    hEnergyDepsCytoplasm->GetYaxis()->SetTitle("Counts / Num. Decays");
    hEnergyDepsNucleus->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    hEnergyDepsNucleus->GetYaxis()->SetTitle("Counts / Num. Decays");


    //------------------–----------
    hEnergyDepsMembrane->Write();
    hEnergyDepsCytoplasm->Write();
    hEnergyDepsNucleus->Write();


    //------------------–----------
    OutputTuplesAnalysis->Write();
    OutputTuplesAnalysis->Close();


}
