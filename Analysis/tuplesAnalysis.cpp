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

    // std::vector<double> col1;
    // std::vector<double> col2;

    if (myfile.is_open())
    {
        std::string line;
        std::string x;
        double y;

        while(std::getline(myfile,line))
        {
            if(line.at(0)=='#')
            {
                continue;
            }
            else
            {
                std::stringstream mysstream(line);
                mysstream >> x >> y;
                input_data.push_back(y);
            }
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
    //------------------–----------
    // Importing data on decay dynamics from Mathematica calculations
    std::vector<double> decayData;

    decayData = importDataStringDouble("../Mathematica/Output/C4-2/Solution/Decays_Activity10_Solution.dat");

    // double totAlphaDecays1h = decayData[decayData.size()];
    double totDecays212PbSolutionFirstHour = decayData[0];

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
    // Making Histograms

    // // // Histogram for average energy deposited in one cell per decay
    // TH1D *hEnergyDepsTotal = new TH1D("hEnergyDepsTotal", "Energy Deposition", 10000, 0.0, 20.0);

    // Histogram for total energy deposited in one nuclei per decay
    TH1D *hEnergyDepsNucleus = new TH1D("hEnergyDepsNucleus", "Energy Deposition in Cell Nuclei / Decay", 10000, 0.0, 20.0);

    // Histogram for total energy deposited in one membrane per decay
    TH1D *hEnergyDepsMembrane = new TH1D("hEnergyDepsMembrane", "Energy Depsition in Cell Membrane / Decay", 10000, 0.0, 20.0);

    // Histogram for total energy deposited in one cytoplasm per decay
    TH1D *hEnergyDepsCytoplasm = new TH1D("hEnergyDepsCytoplasm", "Energy Depsition in Cell Cytoplasm / Decay", 10000, 0.0, 20.0);

    //------------------–----------
    // Making outputfile
    auto OutputTuplesAnalysis = new TFile("OutputTuplesAnalysis.root", "RECREATE");


    //------------------–----------
    // Data analysis

    // Make summation variable to stop the processing of the TTree data once the amount 212Pb decays reaches amount calculated by mathematica
    int decaysProcessed = 0;


    // Looping over every branch/event
    while(myReader.Next())
    {
        decaysProcessed ++;
        //Tuple <cellID, energyDepMembrane, energyDepCytoplasm, energyDepNucleus>; <0,1,2,3>

        std::vector<std::tuple<int,double,double,double>> storedInfo;

        // Looping over every leaf/step
        for(int i=0; i<energyDeps.GetSize(); i++)
        {
            // If first step store the energy
            if(storedInfo.size()==0)
            {
                if(volumeTypes[i]==1)
                {
                    //Membrane
                    storedInfo.push_back(make_tuple(cellIDs[i],energyDeps[i],0.0,0.0));
                }
                else if(volumeTypes[i]==2)
                {
                    // Cytoplasm
                    storedInfo.push_back(make_tuple(cellIDs[i],0.0,energyDeps[i], 0.0));
                }
                if(volumeTypes[i]==3)
                {
                    // Nucleus
                    storedInfo.push_back(make_tuple(cellIDs[i],0.0,0.0,energyDeps[i]));
                }
            }
            // If not first step check for same cellIDS
            else
            {
                // Looping over every tuple stored
                for(int j=0;j<storedInfo.size();j++)
                {
                    // If the cellID has already been stored before then add the enery
                    if(get<0>(storedInfo[j])==cellIDs[i])
                    {
                        // If in membrane
                        if(volumeTypes[i]==1)
                        {
                            get<1>(storedInfo[j]) += energyDeps[i];
                        }
                        // If in cytoplasm
                        else if(volumeTypes[i]==2)
                        {
                            get<2>(storedInfo[j]) += energyDeps[i];
                        }
                        // If in nucleus
                        else if(volumeTypes[i]==3)
                        {
                            get<3>(storedInfo[j]) += energyDeps[i];
                        }
                    }
                }
            }
        }
        // Filling histograms
        for(int k=0;k<storedInfo.size();k++)
        {
           hEnergyDepsMembrane->Fill(get<1>(storedInfo[k]));
           hEnergyDepsCytoplasm->Fill(get<2>(storedInfo[k]));
           hEnergyDepsNucleus->Fill(get<3>(storedInfo[k]));
        }

        if(decaysProcessed >= totDecays212PbSolutionFirstHour)
        {
            break;
        }
    }

    std::cout << decaysProcessed << " , " << totDecays212PbSolutionFirstHour << std::endl;
    // hEnergyDepsTotal->Write();
    hEnergyDepsMembrane->Write();
    hEnergyDepsCytoplasm->Write();
    hEnergyDepsNucleus->Write();
    OutputTuplesAnalysis->Write();
    OutputTuplesAnalysis->Close();

}
