#include <TFile.h>
#include <TH1.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

#include <vector>
#include <iostream>
#include <tuple>

void energyDepsHist()
{
    //------------------–----------
    // Opening TTree file and creating TTreeReader
    std::unique_ptr<TFile> myFile(TFile::Open("../GEANT4Simulations/B4a", "READ"));
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

    // // Histogram for average energy deposited in one cell per decay
    TH1D *hEnergyDepsTotal = new TH1D("hEnergyDepsTotal", "Energy Deposition", 10000, 0.0, 20.0);

    // Histogram for total energy deposited in one nuclei per decay
    TH1D *hEnergyDepsNucleus = new TH1D("hEnergyDepsNucleus", "Energy Deposition in Cell Nuclei / Decay", 10000, 0.0, 20.0);

    // Histogram for total energy deposited in one membrane per decay
    TH1D *hEnergyDepsMembrane = new TH1D("hEnergyDepsMembrane", "Energy Depsition in Cell Membrane / Decay", 10000, 0.0, 20.0);

    // Histogram for total energy deposited in one cytoplasm per decay
    TH1D *hEnergyDepsCytoplasm = new TH1D("hEnergyDepsCytoplasm", "Energy Depsition in Cell Cytoplasm / Decay", 10000, 0.0, 20.0);

    //------------------–----------
    // Making outputfile
    auto outputFile2 = new TFile("Output2.root", "RECREATE");


    //------------------–----------
    // Data analysis

    std::vector<std::tuple<int,int,double>> info;

    // Looping over every branch/event
    while(myReader.Next())
    {
        // Making vectors to store the cell ID of the nuclei where the energy is deposited, and vectors to store the total energy deposited
        std::vector<int> cellIDsNucleusStored;
        std::vector<int> cellIDsMembraneStored;
        std::vector<int> cellIDsCytoplasmStored;
        std::vector<int> cellIDsTotalStored;

        std::vector<double> energyDepsNucleusStored;
        std::vector<double> energyDepsMembraneStored;
        std::vector<double> energyDepsCytoplasmStored;
        std::vector<double> energyDepsTotalStored;

        // Looping over every leaf/step
        for(int i=0; i<energyDeps.GetSize(); i++)
        {
            //------------------–----------
            // Storing energy deposited in any part of the cell
            if(cellIDsTotalStored.size()==0)
            {
                cellIDsTotalStored.push_back(cellIDs[i]);
                energyDepsTotalStored.push_back(energyDeps[i]);
            }
            else
            {
                for(int j=0; j<cellIDsTotalStored.size(); j++)
                {
                    if(cellIDsTotalStored[j]==cellIDs[i])
                    {
                        energyDepsTotalStored[j] += energyDeps[i];
                    }
                    else
                    {
                        cellIDsTotalStored.push_back(cellIDs[i]);
                        energyDepsTotalStored.push_back(energyDeps[i]);
                    }
                }
            }

            //------------------–----------
            // Checking if the energy was deposited in a cell nucleus
            if(volumeTypes[i] == 3)
            {
                // Storing first energy deposition
                if(cellIDsNucleusStored.size()==0)
                {
                    cellIDsNucleusStored.push_back(cellIDs[i]);
                    energyDepsNucleusStored.push_back(energyDeps[i]);
                }
                // Storing further energy depositions
                else
                {
                    // Looping over the previously stored cell IDs
                    for(int j=0; j<cellIDsNucleusStored.size(); j++)
                    {
                        // Updating energy deposited to previously stored nucleus
                        if(cellIDsNucleusStored[j]==cellIDs[i])
                        {
                            energyDepsNucleusStored[j] += energyDeps[i];
                        }
                        // Storing energy deposited to not previously stored nucleus
                        else
                        {
                            energyDepsNucleusStored.push_back(energyDeps[i]);
                            cellIDsNucleusStored.push_back(cellIDs[i]);
                        }
                    }
                }
            }

            //------------------–----------
            // Checking if the energy was deposited in a cell cytoplasm
            else if(volumeTypes[i]==2)
            {
                // Storing first energy deposition
                if(cellIDsCytoplasmStored.size()==0)
                {
                    cellIDsCytoplasmStored.push_back(cellIDs[i]);
                    energyDepsCytoplasmStored.push_back(energyDeps[i]);
                }
                // Storing further energy depositions
                else
                {
                    for(int j=0; j<cellIDsCytoplasmStored.size(); j++)
                    {
                        // Updating energy deposited to previously stored nucleus
                        if(cellIDsCytoplasmStored[j]==cellIDs[i])
                        {
                            energyDepsCytoplasmStored[j] += energyDeps[i];
                        }
                        else
                        {
                            energyDepsCytoplasmStored.push_back(energyDeps[i]);
                            cellIDsCytoplasmStored.push_back(cellIDs[i]);
                        }
                    }
                }
            }

            //------------------–----------
            // Checking if the energy was deposited in a cell membrane
            else if(volumeTypes[i]==1)
            {
                // Storing first energy deposition
                if(cellIDsMembraneStored.size()==0)
                {
                    cellIDsMembraneStored.push_back(cellIDs[i]);
                    energyDepsMembraneStored.push_back(energyDeps[i]);
                }
                // Storing further energy depositions
                else
                {
                    for(int j=0; j<cellIDsMembraneStored.size(); j++)
                    {
                        // Updating energy deposited to previously stored nucleus
                        if(cellIDsMembraneStored[j]==cellIDs[i])
                        {
                            energyDepsMembraneStored[j] += energyDeps[i];
                        }
                        else
                        {
                            energyDepsMembraneStored.push_back(energyDeps[i]);
                            cellIDsMembraneStored.push_back(cellIDs[i]);
                        }
                    }
                }
            }

        }
        //------------------–----------
        // Filling histogram with energy deposited in one cell
        for(int k=0; k<cellIDsTotalStored.size(); k++)
        {
            hEnergyDepsTotal->Fill(energyDepsTotalStored[k]);
        }

        //------------------–----------
        // Filling histogram with energy deposited in nuclei
        for(int k=0; k<cellIDsNucleusStored.size(); k++)
        {
            hEnergyDepsNucleus->Fill(energyDepsNucleusStored[k]);
        }

        //------------------–----------
        // Filling histogram with energy deposited in cytoplasm
        for(int k=0; k<cellIDsCytoplasmStored.size(); k++)
        {
            hEnergyDepsCytoplasm->Fill(energyDepsCytoplasmStored[k]);
        }

        //------------------–----------
        // Filling histogram with energy deposited in memebrane
        for(int k=0; k<cellIDsMembraneStored.size(); k++)
        {
            hEnergyDepsMembrane->Fill(energyDepsMembraneStored[k]);
        }


    }

    // int f = tree->GetEntries();
    // std::cout << f << std::endl;
    // hEnergyDepsTotal->Scale(1/tree.GetEntries());
    //------------------–----------
    // Writing to output file
    hEnergyDepsTotal->Write();
    hEnergyDepsMembrane->Write();
    hEnergyDepsCytoplasm->Write();
    hEnergyDepsNucleus->Write();
    outputFile2->Write();
    outputFile2->Close();
}
