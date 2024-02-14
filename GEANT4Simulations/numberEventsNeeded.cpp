#include <TFile.h>
#include <TH1.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TMath.h>
#include <TChain.h>

#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>
#include <ctime>
#include <tuple>

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>



//------------------–----------
std::vector<double> ImportDataStringDouble(std::string filename)
{
    //----------------------
    //  Imports data from file where the first column contains data in string form
    //  and second column contains data in double form. Returns the columns of double
    //  values in the form of a vector

    std::vector<double> input_data;

    std::fstream myfile(filename, ios_base::in);


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
class DecayDynamics
{
    //----------------------
    //  Class to store calculated decay dynamics data from the Mathematica calculations.
    //  Assumes the output file from Mathematica is structured in a specific way

    public:
        DecayDynamics(int activitySample_in, double U0InternalizedPerCell_in, double U0SurfaceBoundPerCell_in, std::string cellLine_in);

        void LoadDataFromMathematicaCalculations(std::string filepathToMathematicaOutput, double ratio);

        double GetNumberDecaysInSolutionFirstHour(){return numberDecays212PbInSolutionFirstHour;};
        double GetNumberDecaysInMembraneTotalTime(){return numberDecays212PbInMembraneTotalTime;};
        double GetNumberDecaysInCytoplasmTotalTime(){return numberDecays212PbInCytoplasmTotalTime;};

        int GetActivity(){return activitySample;};
        std::string GetCellLine(){return cellLine;};

    private:

        // These number of decays are calculated in Mathematica for a sample of 0.2mL in volume
        // The simulation volume has to be taken into account
        double numberDecays212PbInSolutionFirstHour;
        double numberDecays212PbInMembraneTotalTime;
        double numberDecays212PbInCytoplasmTotalTime;

        double U0InternalizedPerCell;
        double U0SurfaceBoundPerCell;
        int activitySample;  // Given in kBq/1mL
        std::string cellLine;
};


DecayDynamics::DecayDynamics(int activitySample_in, double U0InternalizedPerCell_in, double U0SurfaceBoundPerCell_in, std::string cellLine_in)
{
    activitySample = activitySample_in;
    U0InternalizedPerCell = U0InternalizedPerCell_in;
    U0SurfaceBoundPerCell = U0SurfaceBoundPerCell_in;
    cellLine = cellLine_in;

}



void DecayDynamics::LoadDataFromMathematicaCalculations(std::string filepathToMathematicaOutput, double ratio)
{
    //----------------------
    // Loads data from calculations, assuming output files are structured as shown in the file "212PbDecayDynamics.nb"

    std::string filepathSolutionData = filepathToMathematicaOutput + "/" + cellLine + "/Solution/Activity_" + std::to_string(activitySample) + "kBq/NumberDecays.dat";

    // Importing data for decays occuring in solution
    std::vector<double> decayDataSolution = ImportDataStringDouble(filepathSolutionData);
    numberDecays212PbInSolutionFirstHour = decayDataSolution[0]*ratio;


    // If there are no radionuclides internalized there is no output file for "Cells"
    // So only read "Cells" files if there is uptake
    if(U0InternalizedPerCell>0.0)
    {
        // Importing data for decays occuring in cells
        std::string filepathCellData = filepathToMathematicaOutput + "/" + cellLine + "/Cells/Activity_" + std::to_string(activitySample) + "kBq/NumberDecays.dat";
        std::vector<double> decayDataCells = ImportDataStringDouble(filepathCellData);

        numberDecays212PbInMembraneTotalTime = decayDataCells[7]*ratio;
        numberDecays212PbInCytoplasmTotalTime = decayDataCells[14]*ratio;
    }
    else
    {
        numberDecays212PbInMembraneTotalTime = 0.0;
        numberDecays212PbInCytoplasmTotalTime = 0.0;
    }

    // std::cout << cellLine << ", " << activitySample << "kBq, Solution : " << numberDecays212PbInSolutionFirstHour*ratio << " Membrane : " << numberDecays212PbInMembraneTotalTime*ratio << " Cytoplasm : " << numberDecays212PbInCytoplasmTotalTime*ratio << std::endl;
}



void numberEventsNeeded()
{


    // Calculating volume ratio
    double VolumeSample = 0.2*1000; // mm^3
    double volumeCellTube = TMath::Pi()*std::pow(0.5,2.0)*1.0; // mm^3
    double volumeRatio = volumeCellTube/VolumeSample;
    int numberCells = 1000000;

    int numberIterations = 1;

    std::cout << "Volume ratio of simulation : " << volumeRatio << std::endl;

    double ratioInSolution = 0.98;
    double ratioInCytoplasm = 0.68;
    double ratioInMembrane = 0.99;

    //------------------–----------
    // Defining decay dynamics


    // C4-2 Cells

    DecayDynamics decays_A5kBq_C4_2 = DecayDynamics(5,1.14,1.16,"C4_2");
    DecayDynamics decays_A10kBq_C4_2 = DecayDynamics(10,1.97,1.98,"C4_2");
    DecayDynamics decays_A25kBq_C4_2 = DecayDynamics(25,4.78,5.91,"C4_2");
    DecayDynamics decays_A50kBq_C4_2 = DecayDynamics(50,8.94,11.07,"C4_2");
    DecayDynamics decays_A75kBq_C4_2 = DecayDynamics(75,10.79,13.12,"C4_2");
    DecayDynamics decays_A100kBq_C4_2 = DecayDynamics(100,13.16,22.72,"C4_2");
    DecayDynamics decays_A150kBq_C4_2 = DecayDynamics(150,16.40,23.56,"C4_2");


    // PC3 PIP Cells

    DecayDynamics decays_A10kBq_PC3_PIP = DecayDynamics(10, 47., 2., "PC3_PIP");
    DecayDynamics decays_A25kBq_PC3_PIP = DecayDynamics(25, 119., 229., "PC3_PIP");
    DecayDynamics decays_A50kBq_PC3_PIP = DecayDynamics(50, 229., 11., "PC3_PIP");
    DecayDynamics decays_A75kBq_PC3_PIP = DecayDynamics(75, 335., 17., "PC3_PIP");
    DecayDynamics decays_A100kBq_PC3_PIP = DecayDynamics(100, 448., 22., "PC3_PIP");
    DecayDynamics decays_A150kBq_PC3_PIP = DecayDynamics(150, 565., 28., "PC3_PIP");


    // PC3 Flu Cells
    DecayDynamics decays_A10kBq_PC3_Flu = DecayDynamics(10, 0., 0., "PC3_Flu");
    DecayDynamics decays_A25kBq_PC3_Flu = DecayDynamics(25, 0., 0., "PC3_Flu");
    DecayDynamics decays_A50kBq_PC3_Flu = DecayDynamics(50, 0., 0., "PC3_Flu");
    DecayDynamics decays_A75kBq_PC3_Flu = DecayDynamics(75, 0., 0., "PC3_Flu");
    DecayDynamics decays_A100kBq_PC3_Flu = DecayDynamics(100, 0., 0., "PC3_Flu");
    DecayDynamics decays_A150kBq_PC3_Flu = DecayDynamics(150, 0., 0., "PC3_Flu");
    DecayDynamics decays_A1000kBq_PC3_Flu = DecayDynamics(1000,0.,0., "PC3_Flu");
    DecayDynamics decays_A300kBq_PC3_Flu = DecayDynamics(300,0.,0., "PC3_Flu");
    DecayDynamics decays_A500kBq_PC3_Flu = DecayDynamics(500,0.,0., "PC3_Flu");
    DecayDynamics decays_A200kBq_PC3_Flu = DecayDynamics(200,0.,0., "PC3_Flu");
    DecayDynamics decays_A250kBq_PC3_Flu = DecayDynamics(250,0.,0., "PC3_Flu");


    //------------------–----------
    // Loading decay dynamics calculations

    std::string mathematicaOutput = "../Mathematica/Output";

    decays_A5kBq_C4_2.LoadDataFromMathematicaCalculations(mathematicaOutput.c_str(),volumeRatio);
    decays_A10kBq_C4_2.LoadDataFromMathematicaCalculations(mathematicaOutput.c_str(),volumeRatio);
    decays_A25kBq_C4_2.LoadDataFromMathematicaCalculations(mathematicaOutput.c_str(),volumeRatio);
    decays_A50kBq_C4_2.LoadDataFromMathematicaCalculations(mathematicaOutput.c_str(),volumeRatio);
    decays_A75kBq_C4_2.LoadDataFromMathematicaCalculations(mathematicaOutput.c_str(),volumeRatio);
    decays_A100kBq_C4_2.LoadDataFromMathematicaCalculations(mathematicaOutput.c_str(),volumeRatio);
    decays_A150kBq_C4_2.LoadDataFromMathematicaCalculations(mathematicaOutput.c_str(),volumeRatio);

    decays_A10kBq_PC3_PIP.LoadDataFromMathematicaCalculations(mathematicaOutput.c_str(),volumeRatio);
    decays_A25kBq_PC3_PIP.LoadDataFromMathematicaCalculations(mathematicaOutput.c_str(),volumeRatio);
    decays_A50kBq_PC3_PIP.LoadDataFromMathematicaCalculations(mathematicaOutput.c_str(),volumeRatio);
    decays_A75kBq_PC3_PIP.LoadDataFromMathematicaCalculations(mathematicaOutput.c_str(),volumeRatio);
    decays_A100kBq_PC3_PIP.LoadDataFromMathematicaCalculations(mathematicaOutput.c_str(),volumeRatio);
    decays_A150kBq_PC3_PIP.LoadDataFromMathematicaCalculations(mathematicaOutput.c_str(),volumeRatio);

    decays_A10kBq_PC3_Flu.LoadDataFromMathematicaCalculations(mathematicaOutput.c_str(),volumeRatio);
    decays_A25kBq_PC3_Flu.LoadDataFromMathematicaCalculations(mathematicaOutput.c_str(),volumeRatio);
    decays_A50kBq_PC3_Flu.LoadDataFromMathematicaCalculations(mathematicaOutput.c_str(),volumeRatio);
    decays_A75kBq_PC3_Flu.LoadDataFromMathematicaCalculations(mathematicaOutput.c_str(),volumeRatio);
    decays_A100kBq_PC3_Flu.LoadDataFromMathematicaCalculations(mathematicaOutput.c_str(),volumeRatio);
    decays_A150kBq_PC3_Flu.LoadDataFromMathematicaCalculations(mathematicaOutput.c_str(),volumeRatio);
    // decays_A1000kBq_PC3_Flu.LoadDataFromMathematicaCalculations(mathematicaOutput.c_str(),volumeRatio);
    // decays_A300kBq_PC3_Flu.LoadDataFromMathematicaCalculations(mathematicaOutput.c_str(),volumeRatio);
    // decays_A500kBq_PC3_Flu.LoadDataFromMathematicaCalculations(mathematicaOutput.c_str(),volumeRatio);
    // decays_A200kBq_PC3_Flu.LoadDataFromMathematicaCalculations(mathematicaOutput.c_str(),volumeRatio);
    // decays_A250kBq_PC3_Flu.LoadDataFromMathematicaCalculations(mathematicaOutput.c_str(),volumeRatio);



    std::stringstream filename_C4_2;
    std::stringstream filename_PC3_PIP;
    std::stringstream filename_PC3_Flu;

    filename_C4_2 << "Output_NumberDecays_C4_2.csv";
    filename_PC3_PIP << "Output_NumberDecays_PC3_PIP.csv";
    filename_PC3_Flu << "Output_NumberDecays_PC3_Flu.csv";

    //----------------------------
    // Writing C4-2 data to file

    // open the file. If not c++11 use  ss.str().c_str()  instead
    std::ofstream file_C4_2( filename_C4_2.str() );
    if ( !file_C4_2 )
    {
        std::cerr << "Error: failed to create file " << filename_C4_2.str() << '\n';
    }

    file_C4_2 << "5 " << decays_A5kBq_C4_2.GetNumberDecaysInSolutionFirstHour()/ratioInSolution << " " << decays_A5kBq_C4_2.GetNumberDecaysInMembraneTotalTime()/ratioInMembrane << " " << decays_A5kBq_C4_2.GetNumberDecaysInCytoplasmTotalTime()/ratioInCytoplasm <<

    "\n10 " << decays_A10kBq_C4_2.GetNumberDecaysInSolutionFirstHour()/ratioInSolution << " " << decays_A10kBq_C4_2.GetNumberDecaysInMembraneTotalTime()/ratioInMembrane << " " << decays_A10kBq_C4_2.GetNumberDecaysInCytoplasmTotalTime()/ratioInCytoplasm <<

    "\n25 " << decays_A25kBq_C4_2.GetNumberDecaysInSolutionFirstHour()/ratioInSolution << " " << decays_A25kBq_C4_2.GetNumberDecaysInMembraneTotalTime()/ratioInMembrane << " " << decays_A25kBq_C4_2.GetNumberDecaysInCytoplasmTotalTime()/ratioInCytoplasm <<

    "\n50 " << decays_A50kBq_C4_2.GetNumberDecaysInSolutionFirstHour()/ratioInSolution << " " << decays_A50kBq_C4_2.GetNumberDecaysInMembraneTotalTime()/ratioInMembrane << " " << decays_A50kBq_C4_2.GetNumberDecaysInCytoplasmTotalTime()/ratioInCytoplasm <<

    "\n75 " << decays_A75kBq_C4_2.GetNumberDecaysInSolutionFirstHour()/ratioInSolution << " " << decays_A75kBq_C4_2.GetNumberDecaysInMembraneTotalTime()/ratioInMembrane << " " << decays_A75kBq_C4_2.GetNumberDecaysInCytoplasmTotalTime()/ratioInCytoplasm <<

    "\n100 " << decays_A100kBq_C4_2.GetNumberDecaysInSolutionFirstHour()/ratioInSolution << " " << decays_A100kBq_C4_2.GetNumberDecaysInMembraneTotalTime()/ratioInMembrane << " " << decays_A100kBq_C4_2.GetNumberDecaysInCytoplasmTotalTime()/ratioInCytoplasm <<

    "\n150 " << decays_A150kBq_C4_2.GetNumberDecaysInSolutionFirstHour()/ratioInSolution << " " << decays_A150kBq_C4_2.GetNumberDecaysInMembraneTotalTime()/ratioInMembrane << " " << decays_A150kBq_C4_2.GetNumberDecaysInCytoplasmTotalTime()/ratioInCytoplasm << "\n";

    if ( !file_C4_2 )
    {
        std::cerr << "Error: failed to write to file " << filename_C4_2.str() << '\n';
    }

    //----------------------------
    // Writing PC3 PIP data to file

    // open the file. If not c++11 use  ss.str().c_str()  instead
    std::ofstream file_PC3_PIP( filename_PC3_PIP.str() );
    if ( !file_PC3_PIP )
    {
        std::cerr << "Error: failed to create file " << filename_PC3_PIP.str() << '\n';
    }



    file_PC3_PIP << "10 " << decays_A10kBq_PC3_PIP.GetNumberDecaysInSolutionFirstHour()/ratioInSolution << " " << decays_A10kBq_PC3_PIP.GetNumberDecaysInMembraneTotalTime()/ratioInMembrane << " " << decays_A10kBq_PC3_PIP.GetNumberDecaysInCytoplasmTotalTime()/ratioInCytoplasm <<

    "\n25 " << decays_A25kBq_PC3_PIP.GetNumberDecaysInSolutionFirstHour()/ratioInSolution << " " << decays_A25kBq_PC3_PIP.GetNumberDecaysInMembraneTotalTime()/ratioInMembrane << " " << decays_A25kBq_PC3_PIP.GetNumberDecaysInCytoplasmTotalTime()/ratioInCytoplasm <<


    "\n50 " << decays_A50kBq_PC3_PIP.GetNumberDecaysInSolutionFirstHour()/ratioInSolution << " " << decays_A50kBq_PC3_PIP.GetNumberDecaysInMembraneTotalTime()/ratioInMembrane << " " << decays_A50kBq_PC3_PIP.GetNumberDecaysInCytoplasmTotalTime()/ratioInCytoplasm <<


    "\n75 " << decays_A75kBq_PC3_PIP.GetNumberDecaysInSolutionFirstHour()/ratioInSolution << " " << decays_A75kBq_PC3_PIP.GetNumberDecaysInMembraneTotalTime()/ratioInMembrane << " " << decays_A75kBq_PC3_PIP.GetNumberDecaysInCytoplasmTotalTime()/ratioInCytoplasm <<


    "\n100 " << decays_A100kBq_PC3_PIP.GetNumberDecaysInSolutionFirstHour()/ratioInSolution << " " << decays_A100kBq_PC3_PIP.GetNumberDecaysInMembraneTotalTime()/ratioInMembrane << " " << decays_A100kBq_PC3_PIP.GetNumberDecaysInCytoplasmTotalTime()/ratioInCytoplasm <<


    "\n150 " << decays_A150kBq_PC3_PIP.GetNumberDecaysInSolutionFirstHour()/ratioInSolution << " " << decays_A150kBq_PC3_PIP.GetNumberDecaysInMembraneTotalTime()/ratioInMembrane << " " << decays_A150kBq_PC3_PIP.GetNumberDecaysInCytoplasmTotalTime()/ratioInCytoplasm << "\n";

    if ( !file_PC3_PIP )
    {
        std::cerr << "Error: failed to write to file " << filename_PC3_PIP.str() << '\n';
    }


    //----------------------------
    // Writing PC3 Flu data to file

    // open the file. If not c++11 use  ss.str().c_str()  instead
    std::ofstream file_PC3_Flu( filename_PC3_Flu.str() );
    if ( !file_PC3_Flu )
    {
        std::cerr << "Error: failed to create file " << filename_PC3_Flu.str() << '\n';
    }

    file_PC3_Flu << "10 " << decays_A10kBq_PC3_Flu.GetNumberDecaysInSolutionFirstHour()/ratioInSolution << " " << decays_A10kBq_PC3_Flu.GetNumberDecaysInMembraneTotalTime()/ratioInMembrane << " " << decays_A10kBq_PC3_Flu.GetNumberDecaysInCytoplasmTotalTime()/ratioInCytoplasm <<

    "\n25 " << decays_A25kBq_PC3_Flu.GetNumberDecaysInSolutionFirstHour()/ratioInSolution << " " << decays_A25kBq_PC3_Flu.GetNumberDecaysInMembraneTotalTime()/ratioInMembrane << " " << decays_A25kBq_PC3_Flu.GetNumberDecaysInCytoplasmTotalTime()/ratioInCytoplasm <<


    "\n50 " << decays_A50kBq_PC3_Flu.GetNumberDecaysInSolutionFirstHour()/ratioInSolution << " " << decays_A50kBq_PC3_Flu.GetNumberDecaysInMembraneTotalTime()/ratioInMembrane << " " << decays_A50kBq_PC3_Flu.GetNumberDecaysInCytoplasmTotalTime()/ratioInCytoplasm <<


    "\n75 " << decays_A75kBq_PC3_Flu.GetNumberDecaysInSolutionFirstHour()/ratioInSolution << " " << decays_A75kBq_PC3_Flu.GetNumberDecaysInMembraneTotalTime()/ratioInMembrane << " " << decays_A75kBq_PC3_Flu.GetNumberDecaysInCytoplasmTotalTime()/ratioInCytoplasm <<


    "\n100 " << decays_A100kBq_PC3_Flu.GetNumberDecaysInSolutionFirstHour()/ratioInSolution << " " << decays_A100kBq_PC3_Flu.GetNumberDecaysInMembraneTotalTime()/ratioInMembrane << " " << decays_A100kBq_PC3_Flu.GetNumberDecaysInCytoplasmTotalTime()/ratioInCytoplasm <<


    "\n150 " << decays_A150kBq_PC3_Flu.GetNumberDecaysInSolutionFirstHour()/ratioInSolution << " " << decays_A150kBq_PC3_Flu.GetNumberDecaysInMembraneTotalTime()/ratioInMembrane << " " << decays_A150kBq_PC3_Flu.GetNumberDecaysInCytoplasmTotalTime()/ratioInCytoplasm << "\n";


    if ( !file_PC3_Flu )
    {
        std::cerr << "Error: failed to write to file " << filename_PC3_Flu.str() << '\n';
    }

}