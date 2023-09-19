#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>

void makeRunFiles()
{
    /*
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
  run30.mac
  run31.mac
  run32.mac
  run33.mac
  run34.mac
  run35.mac
  run36.mac
  run37.mac
  run38.mac
  run39.mac
  run40.mac
  run41.mac

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

    std::string prefix("run");

    std::string ext(".mac");


    // SOLUTION
    // int n = 60;

    int initialVolume = 1;
    int runs = 20000;

    for ( int i = 10; i <= 59; ++i )
    {
        std::stringstream ss;
        ss << prefix << i << ext;

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
            output = "solution";
        }
        else if(initialVolume==1)
        {
            output = "membrane";
        }
        else if(initialVolume==2)
        {
            output = "cytoplasm";
        }

        // write something to the newly created file
        file << "/run/initialize\n" << "/run/printProgress 10000\n" << "/Sim/setSeed " << i << "\n/Sim/setOutputFileName Output_" << output << "_Thread" << i << ".root\n" << "/Sim/SetInitialRadionuclide_Z 82\n" << "/Sim/SetInitialRadionuclide_A 212\n" << "/Sim/SetInitialRadionuclide_excitationEnergy 0.0\n" << "/Sim/DefineInitialRadionuclide\n" << "/Sim/SetInitialRadionuclide_location "<< initialVolume <<"\n" << "/run/beamOn " << runs << "\n";
        if ( !file )
        {
            std::cerr << "Error: failed to write to file " << ss.str() << '\n';
            break;
        }
    }


}