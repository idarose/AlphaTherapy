#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TMath.h>
#include <TChain.h>
#include <Math/Interpolator.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>
#include <ctime>
#include <tuple>


//------------------–----------
std::tuple<std::vector<double>,std::vector<double>> ImportData(std::string filename)
{
    //----------------------
    //  Imports data from file where both columns contains data
    //  of type in double. Returns a tuple with two vectors corresponding
    //  to the two columns

    std::vector<double> CDF;
    std::vector<double> time;

    std::tuple<std::vector<double>,std::vector<double>> data;

    std::fstream myfile(filename, ios_base::in);


    if (myfile.is_open())
    {
        std::string line;
        double x;
        double y;

        while(std::getline(myfile,line))
        {
            std::stringstream mystream(line);
            mystream >> x >> y;
            time.push_back(x);
            CDF.push_back(y);
        }
    }
    else
    {
        std::cout << "Unable to open file " << filename << std::endl;
    }
    myfile.close();

    data = std::make_tuple(time,CDF);

    return data;

}

void WriteDecayCurvesToFile(std::string cellLine)
{
    std::vector<int> activities;

    if(cellLine=="C4_2")
    {
        activities = {5,10,25,50,75,100,150};
    }
    else
    {
        activities = {10,25,50,75,100,150};
    }

    std::string outputMathematicaCellLine = "../../Mathematica/Output/" + cellLine + "/";

    double tMin;
    double tMax;

    for(auto entry : activities)
    {

        // --------------------------------
        // Importing data from Mathematica into a tuple containing
        // vectors for time and number of decays
        std::string outputMathematicaSolution = outputMathematicaCellLine + "Solution/Activity_" + std::to_string(entry) + "kBq/Decays212Pb212Bi.dat";
        std::tuple<std::vector<double>,std::vector<double>> solutionData = ImportData(outputMathematicaSolution.c_str());

        tMin = *std::min_element(std::get<0>(solutionData).begin(), std::get<0>(solutionData).end());
        tMax = *std::max_element(std::get<0>(solutionData).begin(), std::get<0>(solutionData).end());


        // --------------------------------
        // Interpolating data
        ROOT::Math::Interpolator *intSolution = new ROOT::Math::Interpolator(std::get<0>(solutionData),std::get<1>(solutionData));

        // --------------------------------
        // Functions for PDF
        TF1* fDecaysSolution = new TF1("fDecaysSolution",
            [&](double*x, double *p)
            {
                return intSolution->Eval(x[0]);
            }
            , tMin, tMax, 0);

        fDecaysSolution->SetNpx(26000.0);

        // --------------------------------
        // Writing solution pdf to file
        std::string fileNameSolution = "Decays212Pb212Bi_" + cellLine + "_Solution_Activity_" + std::to_string(entry) + "kBq.root";
        auto* outputFileSolution = new TFile(fileNameSolution.c_str(), "RECREATE");

        fDecaysSolution->Write();

        outputFileSolution->Write();
        outputFileSolution->Close();

        std::cout << cellLine << " , integral Solution, " << entry << "kBq : " << fDecaysSolution->Integral(0.,26.) << std::endl;;

        if(cellLine=="PC3_PIP"||cellLine=="C4_2")
        {

            // --------------------------------
            // Importing data from Mathematica into a tuple containing
            // vectors for time and PDF value
            std::string outputMathematicaCells = outputMathematicaCellLine + "Cells/Activity_" + std::to_string(entry) + "kBq/Decays212Pb212Bi.dat";
            std::tuple<std::vector<double>,std::vector<double>> cellsData = ImportData(outputMathematicaCells.c_str());

            tMin = *std::min_element(std::get<0>(cellsData).begin(), std::get<0>(cellsData).end());
            tMax = *std::max_element(std::get<0>(cellsData).begin(), std::get<0>(cellsData).end());


            // --------------------------------
            // Interpolating data
            ROOT::Math::Interpolator *intCells = new ROOT::Math::Interpolator(std::get<0>(cellsData),std::get<1>(cellsData));


            // Function for interpolated PDF
            TF1* fDecaysCells = new TF1("fDecaysCells",
                [&](double*x, double *p)
                {
                    return intCells->Eval(x[0]);
                }
                , tMin, tMax, 0);

            fDecaysCells->SetNpx(26000.0);

            // --------------------------------
            // Writing cells pdf to file
            std::string fileNameCells = "Decays212Pb212Bi_" + cellLine + "_Cells_Activity_" + std::to_string(entry) + "kBq.root";
            auto* outputFileCells = new TFile(fileNameCells.c_str(), "RECREATE");

            fDecaysCells->Write();
            std::cout << cellLine << " , integral Cells, " << entry << "kBq : " << fDecaysCells->Integral(0.,26.) << std::endl;;

            outputFileCells->Write();
            outputFileCells->Close();
        }
    }
}

void LoadDecayCurves()
{
    WriteDecayCurvesToFile("C4_2");
    WriteDecayCurvesToFile("PC3_Flu");
    WriteDecayCurvesToFile("PC3_PIP");
}






