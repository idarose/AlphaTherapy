#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <cmath>

//------------------–----------
std::tuple<std::vector<int>,std::vector<double>, std::vector<double>, std::vector<double>> ImportData(std::string filename)
{
    //----------------------
    //  Imports data from file where both columns contains data
    //  of type in double. Returns a tuple with two vectors corresponding
    //  to the two columns

    std::vector<int> activity;
    std::vector<double> decaysSolution;
    std::vector<double> decaysMembrane;
    std::vector<double> decaysCytoplasm;

    std::tuple<std::vector<int>,std::vector<double>, std::vector<double>, std::vector<double>> data;

    std::fstream myfile(filename, ios_base::in);


    if (myfile.is_open())
    {
        std::string line;

        int c1;
        double c2;
        double c3;
        double c4;

        while(std::getline(myfile,line))
        {
            std::stringstream mystream(line);
            mystream >> c1 >> c2 >> c3 >> c4;
            // std::cout << c2 << std::endl;
            activity.push_back(c1);
            decaysSolution.push_back(c2);
            decaysMembrane.push_back(c3);
            decaysCytoplasm.push_back(c4);
        }
    }
    else
    {
        std::cout << "Unable to open file " << filename << std::endl;
    }
    myfile.close();

    data = std::make_tuple(activity,decaysSolution,decaysMembrane,decaysCytoplasm);

    return data;

}


void GenerateFiles(int initialVolume, int runs, int activity, std::string cellLineName, int numberIterations, int caseNumber)
{
    // Defining simulation region
    std:string region;
    if(initialVolume==0)
        {
            region = "Solution";
        }
        else if(initialVolume==1)
        {
            region = "Membrane";
        }
        else if(initialVolume==2)
        {
            region = "Cytoplasm";
        }

    for(int i=0; i<numberIterations;i++)
    {
         // Making output file names
        std::string prefix("Run_");
        std::string ext(".mac");

        std::stringstream ss;
        ss << prefix << cellLineName << "_" << region << "_Case_" << caseNumber + i << ext;

        std::cout << ss.str() << std::endl;

        int seedVal = runs + numberIterations;

        // open the file. If not c++11 use  ss.str().c_str()  instead
        std::ofstream file( ss.str() );
        if ( !file )
        {
            std::cerr << "Error: failed to create file " << ss.str() << '\n';
        }
        // write something to the newly created file
        file << "/run/initialize\n" << "/run/printProgress 10000\n" << "/Sim/setSeed " << seedVal << "\n/Sim/setOutputFileName Output_Pb212_" << cellLineName << "_Activity" << activity << "kBq_" << region << "_Thread_" << i << ".root\n" << "/Sim/SetInitialRadionuclide_Z 82\n" << "/Sim/SetInitialRadionuclide_A 212\n" << "/Sim/SetInitialRadionuclide_excitationEnergy 0.0\n" << "/Sim/DefineInitialRadionuclide\n" << "/Sim/SetInitialRadionuclide_location "<< initialVolume <<"\n" << "/Sim/SetSampleActivity " << activity << "\n/Sim/SetCellLineName " << cellLineName << "\n/Sim/SetDecayCurveRadionuclide\n" << "/run/beamOn " << runs << "\n";

        if ( !file )
        {
            std::cerr << "Error: failed to write to file " << ss.str() << '\n';
        }
    }

}

void makeRunFiles()
{
    std::tuple<std::vector<int>,std::vector<double>, std::vector<double>, std::vector<double>> data_C4_2 = ImportData("../Output_NumberDecays_C4_2.csv");
    std::tuple<std::vector<int>,std::vector<double>, std::vector<double>, std::vector<double>> data_PC3_PIP = ImportData("../Output_NumberDecays_PC3_PIP.csv");
    std::tuple<std::vector<int>,std::vector<double>, std::vector<double>, std::vector<double>> data_PC3_Flu = ImportData("../Output_NumberDecays_PC3_Flu.csv");

    int activity;
    double decaysSolution;
    double decaysMembrane;
    double decaysCytoplasm;

    int numInterations = 10;

    //-------------------------
    // Making files for C4-2

    int caseNum = 0;
    for(int i = 0; i < std::get<0>(data_C4_2).size(); i++)
    {
        activity = std::get<0>(data_C4_2)[i];
        decaysSolution = std::ceil(std::get<1>(data_C4_2)[i]);
        decaysMembrane = std::ceil(std::get<2>(data_C4_2)[i]);
        decaysCytoplasm = std::ceil(std::get<3>(data_C4_2)[i]);


        GenerateFiles(0, decaysSolution, activity, "C4_2", numInterations, caseNum);
        GenerateFiles(1, decaysMembrane, activity, "C4_2", numInterations, caseNum);
        GenerateFiles(2, decaysCytoplasm, activity, "C4_2", numInterations, caseNum);
        caseNum += numInterations;
    }


    caseNum = 0;
    for(int i = 0; i < std::get<0>(data_PC3_PIP).size(); i++)
    {
        //-------------------------
        // Making files for PC3 PIP
        activity = std::get<0>(data_PC3_PIP)[i];
        decaysSolution = std::ceil(std::get<1>(data_PC3_PIP)[i]);
        decaysMembrane = std::ceil(std::get<2>(data_PC3_PIP)[i]);
        decaysCytoplasm = std::ceil(std::get<3>(data_PC3_PIP)[i]);

        GenerateFiles(0, decaysSolution, activity, "PC3_PIP", numInterations, caseNum);
        GenerateFiles(1, decaysMembrane, activity, "PC3_PIP", numInterations, caseNum);
        GenerateFiles(2, decaysCytoplasm, activity, "PC3_PIP", numInterations, caseNum);

        caseNum += numInterations;
    }

    caseNum = 0;
    for(int i = 0; i < std::get<0>(data_PC3_Flu).size(); i++)
    {
        //-------------------------
        // Making files for PC3 Flu
        activity = std::get<0>(data_PC3_Flu)[i];
        decaysSolution = std::ceil(std::get<1>(data_PC3_Flu)[i]);


        GenerateFiles(0, decaysSolution, activity, "PC3_Flu", numInterations, caseNum);

        caseNum += numInterations;
    }
}


/*
run0.mac
  run1.mac
  run2.mac
  run3.mac
  run4.mac
  run5.mac
  run6.mac
  run7.mac
  run8.mac
  run9.mac
  run10.mac
  run11.mac
  run12.mac
  run13.mac
  run14.mac
  run15.mac
  run16.mac
  run17.mac
  run18.mac
  run19.mac
  run20.mac
  run21.mac
  run22.mac
  run23.mac
  run24.mac
  run25.mac
  run26.mac
  run27.mac
  run28.mac
  run29.mac
  */