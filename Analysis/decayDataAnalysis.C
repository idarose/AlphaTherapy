#include <TFile.h>
#include <TH1.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <vector>

#include <random>
#include "TRandom.h"



std::vector<std::vector<double>> importDataDoubles(std::string filename)
{
    std::vector<std::vector<double>> input_data;

    std::fstream myfile;
    myfile.open(filename);

    if (myfile.is_open())
    {
        std::string line;
        double x,y;

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
                input_data.push_back({x,y});
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

void dataAnalysis()
{
    //------------------–----------
    // Opening TTree file and creating TTreeReader
    std::unique_ptr<TFile> myFile(TFile::Open("../B4a-build/B4.root", "READ"));
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
    // Importing data on decay dynamics from Mathematica calculations
    std::vector<std::vector<double>> N212PbData;
    std::vector<std::vector<double>> N212PbTotalData;
    std::vector<double> decayData;

    // Data on N212Pb in solution
    N212PbData = importDataDoubles("../Mathematica/Output/C4-2/Solution/Activity10_N212Pb.dat");

    // Data on N212Pb in total in the system
    N212PbTotalData = importDataDoubles("../Mathematica/Output/C4-2/Solution/Activity10_N212PbTotal.dat");

    decayData = importDataStringDouble("../Mathematica/Output/C4-2/Solution/Decays_Activity10_Solution.dat");


    //------------------–----------
    // Calculating total number of decays in solution in first hour

    // Total number of alpha decays in solution in first hour
    double totAlphaDecays1h = decayData[decayData.size()];
    double totDecays212PbSolution1h = decayData[0];

    // Sum variable to calculate total number of decays
    double totDecays1h = 0.0;

    // Summing
    for(int i=0; i<(decayData.size()-1); i++)
    {
        totDecays1h += decayData[i];
    }

    std::cout << "Total number of decays of 212Pb in solution in first 1h : " <<totDecays212PbSolution1h << std::endl;

    int totDecays1hInt = totDecays1h;
    int totDecays212PbSolution1hInt = totDecays212PbSolution1h;

    //------------------–----------
    // Calculating ratio N212PbSolution/N212PbTotal

    std::vector<double> ratio;
    std::vector<double> t;

    // Calculating ratio and storing ratio as well as time
    for(int i=0; i < N212PbData.size(); i++)
    {
        ratio.push_back(N212PbData[i][1]/N212PbTotalData[i][1]);
        t.push_back(N212PbData[i][0]);
    }


    //------------------–----------
    // Making an interpolating function for the ratio function
    auto *interP = new ROOT::Math::Interpolator(t,ratio);

    // int n = ratio.size();
    // double x[n], interY[n], dataY[n];
    // for (int i=0; i<n; i++) {
    //     x[i] = t[i];
    //     interY[i] = interP->Eval(t[i]);
    //     dataY[i] = ratio[i];
    // }

    // auto *grInterpolatedY = new TGraph(n, x, interY);
    // auto *grDataY = new TGraph(n, x, dataY);

    // grInterpolatedY->SetLineColor(4);
    // grInterpolatedY->Draw("AC*");

    // grDataY->SetLineColor(2);
    // grDataY->Draw("CP");



    //------------------–----------
    // Determining if decay happens in cell or in solution

    // Vectors to store energy depositions from decays in solution and decays in cytoplasm
    std::vector<double> energyDepDecaysSolution;
    std::vector<double> energyDepDecaysCytoplasm;


    // Creating random number generator
    TRandom3 rand3;


    // Interaction time
    double intTime_h; //hours
    // Ratio evaluated at interaction time
    double ratioIntTime;

    int s = 0;
    // Looping over every branch/event
    while(myReader.Next() && s <= totDecays212PbSolution1hInt)
    {
        // Looping over every leaf/step
        for(int i=0; i<energyDeps.GetSize(); i++)
        {
            intTime_h = interactionTime[i]/3600;

            // If interaction time happened in first hour check if happened in solution or in cytoplasm
            if(intTime_h <= 1.0)
            {
                if(i == 0)
                {
                    s++;
                }
                // Calculating ratio at the interaction time
                ratioIntTime = interP->Eval(intTime_h);

                // If random uniform (0,1) falles below ratio value, count as decay in solution. If not count as decay in cytoplasm
                if(rand3.Uniform(1) <= ratioIntTime)
                {
                    // Decay happened in solution
                    energyDepDecaysSolution.push_back(energyDeps[i]);
                }
                else
                {
                    // Decay happened in cytoplasm
                    energyDepDecaysCytoplasm.push_back(energyDeps[i]);
                }
            }

        }
    }
    std::cout << "Number of 212Pb decays in solution in first hour : " << s << std::endl;

}

