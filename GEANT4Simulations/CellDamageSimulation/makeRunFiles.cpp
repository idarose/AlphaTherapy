#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>


void GenerateFiles(int initialVolume, int runs)
{
    std::string prefix("run");

    std::string ext(".mac");

    for ( int i = 0; i <= 9; ++i )
    {
        std::stringstream ss;
        int runNumber = i + initialVolume*10;
        ss << prefix << runNumber << ext;

        int seedVal = runNumber + 1;

        // open the file. If not c++11 use  ss.str().c_str()  instead
        std::ofstream file( ss.str() );
        if ( !file )
        {
            std::cerr << "Error: failed to create file " << ss.str() << '\n';
            break;
        }

        std:string output;

        if(initialVolume==0)
        {
            output = "Solution";
        }
        else if(initialVolume==1)
        {
            output = "Membrane";
        }
        else if(initialVolume==2)
        {
            output = "Cytoplasm";
        }

        // write something to the newly created file
        file << "/run/initialize\n" << "/run/printProgress 10000\n" << "/Sim/setSeed " << runNumber << "\n/Sim/setOutputFileName Output_" << output << "_Thread_" << i << ".root\n" << "/Sim/SetInitialRadionuclide_Z 82\n" << "/Sim/SetInitialRadionuclide_A 212\n" << "/Sim/SetInitialRadionuclide_excitationEnergy 0.0\n" << "/Sim/DefineInitialRadionuclide\n" << "/Sim/SetInitialRadionuclide_location "<< initialVolume <<"\n" << "/run/beamOn " << runs << "\n";
        if ( !file )
        {
            std::cerr << "Error: failed to write to file " << ss.str() << '\n';
            break;
        }
    }
}

void makeRunFiles()
{
    /*
    solution : Submitted batch job 9239377
    membrane : Submitted batch job 9239390
    cytoplasm : Submitted batch job 9239395






    /run/initialize
    /run/printProgress 10000
    /Sim/setSeed 1
    /Sim/setOutputFileName Output_thread1.root

    /Sim/SetInitialRadionuclide_Z 82
    /Sim/SetInitialRadionuclide_A 212
    /Sim/SetInitialRadionuclide_excitationEnergy 0.0
    /Sim/DefineInitialRadionuclide

    /Sim/SetInitialRadionuclide_location 0

    /run/beamOn 10000


    # slurm-9224236_1.out. 10000 runs, 40 threads, 60.0 s


    /run/initialize
    /run/printProgress 10000
    /Sim/setSeed 2
    /Sim/setOutputFileName Output_thread2.root

    /Sim/SetInitialRadionuclide_Z 82
    /Sim/SetInitialRadionuclide_A 212
    /Sim/SetInitialRadionuclide_excitationEnergy 0.0
    /Sim/DefineInitialRadionuclide

    /Sim/SetInitialRadionuclide_location 1

    /run/beamOn 10000


    #slurm-9224531_1.out. 10000 runs, 40 threads, 52.0 s



    /run/initialize
    /run/printProgress 10000
    /Sim/setSeed 2
    /Sim/setOutputFileName Output_thread2.root

    /Sim/SetInitialRadionuclide_Z 82
    /Sim/SetInitialRadionuclide_A 212
    /Sim/SetInitialRadionuclide_excitationEnergy 0.0
    /Sim/DefineInitialRadionuclide

    /Sim/SetInitialRadionuclide_location 2

    /run/beamOn 10000


    #slurm-9224531_1.out. 10000 runs, 40 threads, 45.0 s

    */

    GenerateFiles(0, 385000);
    GenerateFiles(1, 71590);
    GenerateFiles(2, 52000);

}