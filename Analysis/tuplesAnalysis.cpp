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

std::vector<std::tuple<int,double,double,double>> organizeDataFromTTree(const char* filename)
{
    /*Returns Tuple  with data organized as : <cellID, energyDepMembrane, energyDepCytoplasm, energyDepNucleus>*/

    //------------------–----------
    // Opening TTree file and creating TTreeReader
    std::unique_ptr<TFile> myFile(TFile::Open(filename, "READ"));
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

    std::vector<std::tuple<int,double,double,double>> storedInfoAllEvents;

    // Looping over every branch/event
    while(myReader.Next())
    {
        std::vector<std::tuple<int,double,double,double>> storedInfoEvent;

        // Looping over every leaf/step
        for(int i=0; i<energyDeps.GetSize(); i++)
        {
            // If first step store the energy
            if(storedInfoEvent.size()==0)
            {
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
        // Storing info from one event to tuple storing info from all events
        for(int k=0; k<storedInfoEvent.size(); k++)
        {
            storedInfoAllEvents.push_back(storedInfoEvent[k]);
        }
    }
    return storedInfoAllEvents;
}


void tuplesAnalysis()
{
    std::vector<std::tuple<int,double,double,double>> storedInfoDecaysSolution;
    // std::vector<std::tuple<int,double,double,double>> storedInfoDecaysMembrane;
    // std::vector<std::tuple<int,double,double,double>> storedInfoDecaysCytoplasm;


    storedInfoDecaysSolution = organizeDataFromTTree("../GEANT4Simulations/B4aSolution-build/B4.root");
    // storedInfoDecaysMembrane = organizeDataFromTTree("../GEANT4Simulations/B4aMembrane-build/B4.root");
    // storedInfoDecaysCytoplasm = organizeDataFromTTree("../GEANT4Simulations/B4aCytoplasm-build/B4.root");

     // Histogram for total energy deposited in one nuclei per decay
    TH1D *hEnergyDepsNucleus = new TH1D("hEnergyDepsNucleus", "Energy Deposition in Cell Nuclei / Decay", 10000, 0.0, 20.0);

    // Histogram for total energy deposited in one membrane per decay
    TH1D *hEnergyDepsMembrane = new TH1D("hEnergyDepsMembrane", "Energy Depsition in Cell Membrane / Decay", 10000, 0.0, 20.0);

    // Histogram for total energy deposited in one cytoplasm per decay
    TH1D *hEnergyDepsCytoplasm = new TH1D("hEnergyDepsCytoplasm", "Energy Depsition in Cell Cytoplasm / Decay", 10000, 0.0, 20.0);

    //------------------–----------
    // Making outputfile
    auto OutputTuplesAnalysis = new TFile("OutputTuplesAnalysis.root", "RECREATE");


    //Filling histograms
    for(int k=0; k<storedInfoDecaysSolution.size(); k++)
    {
       hEnergyDepsMembrane->Fill(get<1>(storedInfoDecaysSolution[k]));
       // hEnergyDepsMembrane->Fill(get<1>(storedInfoDecaysMembrane[k]));
       // hEnergyDepsMembrane->Fill(get<1>(storedInfoDecaysCytoplasm[k]));

       hEnergyDepsCytoplasm->Fill(get<2>(storedInfoDecaysSolution[k]));
       hEnergyDepsNucleus->Fill(get<3>(storedInfoDecaysSolution[k]));
    }

    hEnergyDepsMembrane->Scale(1.0/storedInfoDecaysSolution.size());
    hEnergyDepsCytoplasm->Scale(1.0/storedInfoDecaysSolution.size());
    hEnergyDepsNucleus->Scale(1.0/storedInfoDecaysSolution.size());


    hEnergyDepsMembrane->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    hEnergyDepsMembrane->GetYaxis()->SetTitle("Counts / Num. Decays");
    hEnergyDepsCytoplasm->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    hEnergyDepsCytoplasm->GetYaxis()->SetTitle("Counts / Num. Decays");
    hEnergyDepsNucleus->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    hEnergyDepsNucleus->GetYaxis()->SetTitle("Counts / Num. Decays");


    hEnergyDepsMembrane->Write();
    hEnergyDepsCytoplasm->Write();
    hEnergyDepsNucleus->Write();

    OutputTuplesAnalysis->Write();
    OutputTuplesAnalysis->Close();


}
