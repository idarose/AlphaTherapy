#include "include/DecayDynamics.hpp"
#include "include/CellHit.hpp"
#include "include/EnergyDepositionHistograms.hpp"
#include "include/AddEnergyDepositionHistograms.hpp"
#include <future>
#include <thread>
#include <TROOT.h>

//------------------–----------
int AnalyzeHistogramsFromSimulation(DecayDynamics decayDynamicsInstance, int numberIterations)
{
    std::cout << "---------------------------------" << std::endl;
    std::cout << "Analyzing Histograms for : " << decayDynamicsInstance.GetCellLine() << ", Activity: " << decayDynamicsInstance.GetActivity() << std::endl;
    //------------------–----------
    // Loading decay dynamics
    double numberDecays212PbInSolution1hTo2h = decayDynamicsInstance.GetNumberDecaysInSolutionFirstHour()*decayDynamicsInstance.GetVolumeRatio();
    double numberDecays212PbInMembrane1hTo26h = decayDynamicsInstance.GetNumberDecaysInMembraneTotalTime()*decayDynamicsInstance.GetVolumeRatio();
    double numberDecays212PbInCytoplasm1hTo26h = decayDynamicsInstance.GetNumberDecaysInCytoplasmTotalTime()*decayDynamicsInstance.GetVolumeRatio();

    double MaxEnergyCytoplasm = 100.;

    //---------------------------
    // Making threads
    std::vector<std::thread> threads;
    std::vector<std::promise<TFile*>> promises(numberIterations);
    std::vector<std::future<TFile*>> futures;

    for (auto& promise : promises) {
        futures.push_back(promise.get_future());
    }

    //------------------–----------
    //  Function for make one histogram for one iteration
    auto MakeFileCellHitsOneIteration = [&](TFile *inputFile_solution, TFile *inputFile_membrane, TFile *inputFile_cytoplasm, int itNum)
    {

        //------------------–----------
        // Generating empty histograms
        EnergyDepositionHistograms energyDepHistograms = EnergyDepositionHistograms(MaxEnergyCytoplasm, itNum);
        energyDepHistograms.GenerateEmptyHistograms(decayDynamicsInstance);

        //------------------–----------
        // Opening TTree files and creating TTreeReaders

        // Reader for solution simulation
        auto treeSolutionSim = inputFile_solution->Get<TTree>("B4");
        TTreeReader myReaderSolutionSim(treeSolutionSim);


        // Reader for membrane simulation
        auto treeMembraneSim = inputFile_membrane->Get<TTree>("B4");
        TTreeReader myReaderMembraneSim(treeMembraneSim);


        // Reader for cytoplasm simulation
        auto treeCytoplasmSim = inputFile_cytoplasm->Get<TTree>("B4");
        TTreeReader myReaderCytoplasmSim(treeCytoplasmSim);


        //------------------–----------
        // Accessing brances of tree

        // Solution
        TTreeReaderArray<double> energyDepsSolutionSim(myReaderSolutionSim, "EnergyDeps");
        TTreeReaderArray<int> volumeTypesSolutionSim(myReaderSolutionSim, "VolumeTypes");
        TTreeReaderArray<int> cellIDsSolutionSim(myReaderSolutionSim, "CellIDs");
        TTreeReaderArray<double> kineticEnergySolutionSim(myReaderSolutionSim, "KineticEnergy");
        TTreeReaderArray<int> particleTypeSolutionSim(myReaderSolutionSim, "ParticleType");
        TTreeReaderArray<double> interactionTimeSolutionSim(myReaderSolutionSim, "InteractionTime");
        TTreeReaderArray<double> firstInteractionTimeSolutionSim(myReaderSolutionSim, "FirstInteractionTime");
        TTreeReaderArray<int> firstInteractionVolumeSolutionSim(myReaderSolutionSim, "FirstInteractionVolume");
        TTreeReaderArray<int> trackIDSolutionSim(myReaderSolutionSim, "TrackID");
        TTreeReaderArray<int> parentIDSolutionSim(myReaderSolutionSim, "ParentID");

        // Membrane
        TTreeReaderArray<double> energyDepsMembraneSim(myReaderMembraneSim, "EnergyDeps");
        TTreeReaderArray<int> volumeTypesMembraneSim(myReaderMembraneSim, "VolumeTypes");
        TTreeReaderArray<int> cellIDsMembraneSim(myReaderMembraneSim, "CellIDs");
        TTreeReaderArray<double> kineticEnergyMembraneSim(myReaderMembraneSim, "KineticEnergy");
        TTreeReaderArray<int> particleTypeMembraneSim(myReaderMembraneSim, "ParticleType");
        TTreeReaderArray<double> interactionTimeMembraneSim(myReaderMembraneSim, "InteractionTime");
        TTreeReaderArray<double> firstInteractionTimeMembraneSim(myReaderMembraneSim, "FirstInteractionTime");
        TTreeReaderArray<int> firstInteractionVolumeMembraneSim(myReaderMembraneSim, "FirstInteractionVolume");
        TTreeReaderArray<int> trackIDMembraneSim(myReaderMembraneSim, "TrackID");
        TTreeReaderArray<int> parentIDMembraneSim(myReaderMembraneSim, "ParentID");

        // Cytoplasm
        TTreeReaderArray<double> energyDepsCytoplasmSim(myReaderCytoplasmSim, "EnergyDeps");
        TTreeReaderArray<int> volumeTypesCytoplasmSim(myReaderCytoplasmSim, "VolumeTypes");
        TTreeReaderArray<int> cellIDsCytoplasmSim(myReaderCytoplasmSim, "CellIDs");
        TTreeReaderArray<double> kineticEnergyCytoplasmSim(myReaderCytoplasmSim, "KineticEnergy");
        TTreeReaderArray<int> particleTypeCytoplasmSim(myReaderCytoplasmSim, "ParticleType");
        TTreeReaderArray<double> interactionTimeCytoplasmSim(myReaderCytoplasmSim, "InteractionTime");
        TTreeReaderArray<double> firstInteractionTimeCytoplasmSim(myReaderCytoplasmSim, "FirstInteractionTime");
        TTreeReaderArray<int> firstInteractionVolumeCytoplasmSim(myReaderCytoplasmSim, "FirstInteractionVolume");
        TTreeReaderArray<int> trackIDCytoplasmSim(myReaderCytoplasmSim, "TrackID");
        TTreeReaderArray<int> parentIDCytoplasmSim(myReaderCytoplasmSim, "ParentID");


        //------------------–----------
        // Vector to store cell hits
        std::vector<CellHit> storedCellHits;

        int decayOriginSolution = 0;
        int decayOriginMembrane = 1;
        int decayOriginCytoplasm = 2;

        //------------------–----------
        //Looping through data for decays occurring in solution from hour 1 to hour 2

        // Counter variable
        int numberDecays212PbInSolution1hTo2h_counter = 0;
        bool whileLoopSolutionSimWasBroken = false;

        while(myReaderSolutionSim.Next())
        {

            //--------------------------------------
            // bool for if the decay occurred in solution from hour 1 to hour 2
            bool firstInteractionInSolution1hTo2h = false;

            // Checking when and where first interaction occured
            if(firstInteractionVolumeSolutionSim[0]==0)
            {

                if(firstInteractionTimeSolutionSim[0]/3600. >= 1. && firstInteractionTimeSolutionSim[0]/3600. <= 2.)
                {

                    // Sett bool true
                    firstInteractionInSolution1hTo2h = true;

                    // Update counter
                    numberDecays212PbInSolution1hTo2h_counter ++;
                }
            }

            //------------------–----------
            // Break loop if number of decays have been reached
            if(numberDecays212PbInSolution1hTo2h_counter >= numberDecays212PbInSolution1hTo2h)
            {
                std::cout << "Solution finished" << std::endl;
                whileLoopSolutionSimWasBroken = true;
                break;
            }

            //--------------------------------------
            // Vector to store <cellID, trackID, volume hit> for particle hits to a cell for one event
            std::vector<std::tuple<int,int,int>> particleHitsInfoForEventVec;

            // Sorting through interactions
            if(firstInteractionInSolution1hTo2h)
            {
                //--------------------------------------
                // looping over all steps for one event/decay
                for(int i=0; i<energyDepsSolutionSim.GetSize(); i++)
                {

                    // Only add if actual energy deposition
                    if(energyDepsSolutionSim[i]>0.)
                    {

                        // Only add if interaction happened between hour 1 and hour 2
                        if(interactionTimeSolutionSim[i]/3600.0 >= 1. && interactionTimeSolutionSim[i]/3600.0 <= 2.)
                        {

                            // boolean: true if cell has not already been hit before
                            // false if already hit
                            bool cellAlreadyHit = false;


                            //--------------------------------------
                            // Loop over all previous cell hits stored
                            for(int ii=0; ii<storedCellHits.size(); ii++)
                            {

                                //--------------------------------------
                                // Add energy deposition to cell that has already been hit before
                                if(cellIDsSolutionSim[i]==storedCellHits[ii].GetCellID())
                                {
                                    // Add energy deposition
                                    // std::cout << "Interaction volume " << volumeTypesSolutionSim[i] << std::endl;
                                    storedCellHits[ii].AddEnergyDeposition(energyDepsSolutionSim[i],volumeTypesSolutionSim[i],decayOriginSolution, particleTypeSolutionSim[i]);

                                    // Update boolean
                                    cellAlreadyHit = true;


                                    //---------------------------------------------
                                    // boolean: true if the particle track has already been stored for a specific cell ID in this specific component of the cell
                                    bool trackAlreadyRegisteredInThisCellComponent = false;

                                    // boolean: true if it is the first time the track has been counted in a cell
                                    bool firstTimeCountingTrackInCell = true;

                                    // Looping over all previously stored track IDs and corresponding cell IDs for this event
                                    for(int iii=0; iii<particleHitsInfoForEventVec.size(); iii++)
                                    {

                                        // If track ID is already stored for this cell ID, do not update particle hit counts for this cell ID
                                        if(cellIDsSolutionSim[i] == std::get<0>(particleHitsInfoForEventVec[iii]))
                                        {
                                            if(trackIDSolutionSim[i]==std::get<1>(particleHitsInfoForEventVec[iii]))
                                            {
                                                // Track has been counted before in this cell, just in a different component of the cell
                                                firstTimeCountingTrackInCell = false;

                                                if(volumeTypesSolutionSim[i]==std::get<2>(particleHitsInfoForEventVec[iii]))
                                                {
                                                    // Track already counted in this cell and cell component
                                                    trackAlreadyRegisteredInThisCellComponent = true;
                                                }
                                            }
                                        }
                                    }
                                    // If track ID has not already been registered for this component in the same cell
                                    if(!trackAlreadyRegisteredInThisCellComponent)
                                    {

                                        // If alpha particle type
                                        if(particleTypeSolutionSim[i]==1000020040)
                                        {
                                            storedCellHits[ii].HitByAlphaParticle(volumeTypesSolutionSim[i],firstTimeCountingTrackInCell);
                                        }

                                        // Save the hit and cellID
                                        std::tuple<int,int,int> particleHit = std::make_tuple(cellIDsSolutionSim[i],trackIDSolutionSim[i],volumeTypesSolutionSim[i]);
                                        particleHitsInfoForEventVec.push_back(particleHit);

                                    }
                                }
                            }

                            //--------------------------------------
                            // Register a new cell hit, and add energy deposition
                            if(!cellAlreadyHit)
                            {
                                // std::cout << "Interaction volume " << volumeTypesSolutionSim[i] << std::endl;
                                CellHit aNewCellHit = CellHit(cellIDsSolutionSim[i]);
                                aNewCellHit.AddEnergyDeposition(energyDepsSolutionSim[i], volumeTypesSolutionSim[i],decayOriginSolution,particleTypeSolutionSim[i]);
                                storedCellHits.push_back(aNewCellHit);


                                // If alpha particle
                                if(particleTypeSolutionSim[i]==1000020040)
                                {
                                    aNewCellHit.HitByAlphaParticle(volumeTypesSolutionSim[i],true);
                                }

                                std::tuple<int,int,int> particleHit = std::make_tuple(cellIDsSolutionSim[i],trackIDSolutionSim[i],volumeTypesSolutionSim[i]);
                                particleHitsInfoForEventVec.push_back(particleHit);

                            }
                        }
                    }
                }
            }
        }

        //--------------------------------------
        // Check if enough decays were processed
        if(!whileLoopSolutionSimWasBroken)
        {
            std::cout << "Not enough events in solution simulation file! Need " << numberDecays212PbInSolution1hTo2h << " number of decays. Only reached " << numberDecays212PbInSolution1hTo2h_counter << " number of decays at entry number " << myReaderSolutionSim.GetCurrentEntry() << std::endl;
        }


        //------------------–----------
        // Looping through data for decays occuring in Membrane from hour 1 to hour 26

        // Counter variable
        int numberDecays212PbInMembrane1hTo26h_counter = 0;
        bool whileLoopMembraneSimWasBroken = false;

        if(inputFile_membrane)
        {
            while(myReaderMembraneSim.Next())
            {
                //--------------------------------------
                // Chacking if first interaction took place between 1 and 26 hours
                bool firstInteractionInMembrane1hTo26h = false;

                if(firstInteractionTimeMembraneSim[0]/3600. >= 1. && firstInteractionTimeMembraneSim[0]/3600. <= 26.)
                {
                    firstInteractionInMembrane1hTo26h = true;
                    numberDecays212PbInMembrane1hTo26h_counter ++;
                }

                //--------------------------------------
                // Break loop if number of decays have been reached
                if(numberDecays212PbInMembrane1hTo26h_counter >= numberDecays212PbInMembrane1hTo26h)
                {
                    whileLoopMembraneSimWasBroken = true;
                    std::cout << "Membrane finished" << std::endl;
                    break;
                }

                //--------------------------------------
                // Vector to store <cellID, trackID> for particle hits to a cell for one event
                std::vector<std::tuple<int,int,int>> particleHitsInfoForEventVec;

                if(firstInteractionInMembrane1hTo26h)
                {
                    //--------------------------------------
                    // looping over all steps for one event/decay
                    for(int i=0; i<energyDepsMembraneSim.GetSize(); i++)
                    {

                        // Only add if actual energy deposition
                        if(energyDepsMembraneSim[i]>0.)
                        {

                            // Only add if interaction happened between hour 1 and hour 26
                            if(interactionTimeMembraneSim[i]/3600.0 >= 1. && interactionTimeMembraneSim[i]/3600.0 <= 26.)
                            {

                                // boolean: true if cell has already been hit before
                                bool cellAlreadyHit = false;

                                //--------------------------------------
                                // Loop over all previous cell hits stored
                                for(int ii=0; ii<storedCellHits.size(); ii++)
                                {

                                    // Add energy deposition to cell that has already been hit before
                                    if(cellIDsMembraneSim[i]==storedCellHits[ii].GetCellID())
                                    {
                                        // std::cout << "Interaction volume " << volumeTypesSolutionSim[i] << std::endl;

                                        storedCellHits[ii].AddEnergyDeposition(energyDepsMembraneSim[i],volumeTypesMembraneSim[i],decayOriginMembrane, particleTypeMembraneSim[i]);

                                        // Update boolean
                                        cellAlreadyHit = true;

                                        //--------------------------------------
                                        // boolean: true if the particle track has already been stored for a specific cell ID
                                        bool trackAlreadyRegisteredInThisCellComponent = false;
                                        bool firstTimeCountingTrackInCell = true;

                                        // Looping over all previously stored track IDs and corresponding cell IDs for this event
                                        for(int iii=0; iii<particleHitsInfoForEventVec.size(); iii++)
                                        {

                                            // If track ID is already stored for this cell ID, do not update particle hit counts for this cell ID
                                            if(cellIDsMembraneSim[i] == std::get<0>(particleHitsInfoForEventVec[iii]))
                                            {
                                                if(trackIDMembraneSim[i]==std::get<1>(particleHitsInfoForEventVec[iii]))
                                                {
                                                    firstTimeCountingTrackInCell = false;
                                                    if(volumeTypesMembraneSim[i]==std::get<2>(particleHitsInfoForEventVec[iii]))
                                                    {
                                                        trackAlreadyRegisteredInThisCellComponent = true;
                                                    }
                                                }
                                            }
                                        }
                                        // If track ID has not already been registered or not been registered for this cell ID
                                        if(!trackAlreadyRegisteredInThisCellComponent)
                                        {

                                            // If alpha particle type
                                            if(particleTypeMembraneSim[i]==1000020040)
                                            {
                                                storedCellHits[ii].HitByAlphaParticle(volumeTypesMembraneSim[i],firstTimeCountingTrackInCell);
                                            }

                                            // Save the hit and cellID
                                            std::tuple<int,int,int> particleHit = std::make_tuple(cellIDsMembraneSim[i],trackIDMembraneSim[i],volumeTypesMembraneSim[i]);
                                            particleHitsInfoForEventVec.push_back(particleHit);

                                        }
                                    }
                                }

                                // Register a new cell hit, and add energy deposition
                                if(!cellAlreadyHit)
                                {
                                    // std::cout << "Interaction volume " << volumeTypesSolutionSim[i] << std::endl;

                                    CellHit aNewCellHit = CellHit(cellIDsMembraneSim[i]);
                                    aNewCellHit.AddEnergyDeposition(energyDepsMembraneSim[i], volumeTypesMembraneSim[i], decayOriginMembrane,particleTypeMembraneSim[i]);
                                    storedCellHits.push_back(aNewCellHit);

                                    // If alpha particle
                                    if(particleTypeMembraneSim[i]==1000020040)
                                    {
                                        aNewCellHit.HitByAlphaParticle(volumeTypesMembraneSim[i],true);
                                    }

                                    std::tuple<int,int,int> particleHit = std::make_tuple(cellIDsMembraneSim[i],trackIDMembraneSim[i],volumeTypesMembraneSim[i]);
                                    particleHitsInfoForEventVec.push_back(particleHit);
                                }
                            }
                        }
                    }
                }
            }

            //--------------------------------------
            // Checking if enough decays were processed
            if(!whileLoopMembraneSimWasBroken)
            {
                std::cout << "Not enough events in membrane simulation file! Need " << numberDecays212PbInMembrane1hTo26h << " number of decays. Only reached " << numberDecays212PbInMembrane1hTo26h_counter << " number of decays at entry number " << myReaderMembraneSim.GetCurrentEntry() << std::endl;
            }
        }



        // ------------------–----------
        // Looping through data for decays occuring in Cytoplasm between hour 1 and hour 26

        // Counter variable
        int numberDecays212PbInCytoplasm1hTo26h_counter = 0;
        bool whileLoopCytoplasmSimWasBroken = false;

        if(inputFile_cytoplasm)
        {
            while(myReaderCytoplasmSim.Next())
            {

                //--------------------------------------
                // Chacking if first interaction took place between 1 and 26 hours
                bool firstInteractionInCytoplasm1hTo26h = false;

                if(firstInteractionVolumeCytoplasmSim[0] == 2)
                {
                    if(firstInteractionTimeCytoplasmSim[0]/3600. >= 1. && firstInteractionTimeCytoplasmSim[0]/3600. <= 26.)
                    {
                        firstInteractionInCytoplasm1hTo26h = true;
                        numberDecays212PbInCytoplasm1hTo26h_counter ++;
                    }
                }

                //--------------------------------------
                // Break loop if number of decays have been reached
                if(numberDecays212PbInCytoplasm1hTo26h_counter >= numberDecays212PbInCytoplasm1hTo26h)
                {
                    std::cout << "Cytoplasm finished" << std::endl;
                    whileLoopCytoplasmSimWasBroken = true;
                    break;
                }

                //--------------------------------------
                // Vector to store <cellID, trackID> for particle hits to a cell for one event
                std::vector<std::tuple<int,int,int>> particleHitsInfoForEventVec;

                if(firstInteractionInCytoplasm1hTo26h)
                {
                    // looping over all steps for one event/decay
                    for(int i=0; i<energyDepsCytoplasmSim.GetSize(); i++)
                    {

                        // Only add if actual energy deposition
                        if(energyDepsCytoplasmSim[i]>0.)
                        {

                            // Only add if interaction happened within between hour 1 and hour 26
                            if(interactionTimeCytoplasmSim[i]/3600.0 >= 1. && interactionTimeCytoplasmSim[i]/3600.0 <= 26.)
                            {

                                //--------------------------------------
                                // boolean: true if cell has already been hit before in this event
                                bool cellAlreadyHit = false;

                                // Loop over all previous cell hits stored
                                for(int ii=0; ii<storedCellHits.size(); ii++)
                                {

                                    // Add energy deposition to cell that has already been hit before
                                    if(cellIDsCytoplasmSim[i]==storedCellHits[ii].GetCellID())
                                    {
                                        // std::cout << "Interaction volume " << volumeTypesSolutionSim[i] << std::endl;

                                        storedCellHits[ii].AddEnergyDeposition(energyDepsCytoplasmSim[i],volumeTypesCytoplasmSim[i],decayOriginCytoplasm,particleTypeCytoplasmSim[i]);

                                        // Update boolean
                                        cellAlreadyHit = true;

                                        // boolean: true if the particle track has already been stored for a specific cell ID
                                        bool trackAlreadyRegisteredInThisCellComponent = false;
                                        bool firstTimeCountingTrackInCell = true;

                                        // Looping over all previously stored track IDs and corresponding cell IDs for this event
                                        for(int iii=0; iii<particleHitsInfoForEventVec.size(); iii++)
                                        {

                                            // If track ID is already stored for this cell ID, do not update particle hit counts for this cell ID
                                            if(cellIDsCytoplasmSim[i] == std::get<0>(particleHitsInfoForEventVec[iii]))
                                            {
                                                if(trackIDCytoplasmSim[i]==std::get<1>(particleHitsInfoForEventVec[iii]))
                                                {
                                                    firstTimeCountingTrackInCell = false;
                                                    if(volumeTypesCytoplasmSim[i]==std::get<2>(particleHitsInfoForEventVec[iii]))
                                                    {
                                                        trackAlreadyRegisteredInThisCellComponent = true;
                                                    }
                                                }
                                            }
                                        }
                                        // If track ID has not already been registered or not been registered for this cell ID
                                        if(!trackAlreadyRegisteredInThisCellComponent)
                                        {

                                            // If alapha particle type
                                            if(particleTypeCytoplasmSim[i]==1000020040)
                                            {
                                                storedCellHits[ii].HitByAlphaParticle(volumeTypesCytoplasmSim[i],firstTimeCountingTrackInCell);
                                                // std::cout << "New track : "<< cellIDsCytoplasmSim[i] << " , " << trackIDCytoplasmSim[i] << " , " << particleTypeCytoplasmSim[i] << std::endl;
                                            }

                                            // Save the hit and cellID
                                            std::tuple<int,int,int> particleHit = std::make_tuple(cellIDsCytoplasmSim[i],trackIDCytoplasmSim[i],volumeTypesCytoplasmSim[i]);
                                            particleHitsInfoForEventVec.push_back(particleHit);

                                        }
                                    }
                                }
                                // Register a new cell hit, and add energy deposition
                                if(!cellAlreadyHit)
                                {
                                    // std::cout << "Interaction volume " << volumeTypesSolutionSim[i] << std::endl;

                                    CellHit aNewCellHit = CellHit(cellIDsCytoplasmSim[i]);
                                    aNewCellHit.AddEnergyDeposition(energyDepsCytoplasmSim[i], volumeTypesCytoplasmSim[i],decayOriginCytoplasm,particleTypeCytoplasmSim[i]);
                                    storedCellHits.push_back(aNewCellHit);

                                    // If alpha particle
                                    if(particleTypeCytoplasmSim[i]==1000020040)
                                    {
                                        aNewCellHit.HitByAlphaParticle(volumeTypesCytoplasmSim[i],true);
                                    }

                                    std::tuple<int,int,int> particleHit = std::make_tuple(cellIDsCytoplasmSim[i],trackIDCytoplasmSim[i],volumeTypesCytoplasmSim[i]);
                                    particleHitsInfoForEventVec.push_back(particleHit);
                                }
                            }
                        }
                    }

                }
            }

            //--------------------------------------
            // Checking if enough decays where processed
            if(!whileLoopCytoplasmSimWasBroken)
            {
                std::cout << "Not enough events in cytoplasm simulation file! Need " << numberDecays212PbInCytoplasm1hTo26h << " number of decays. Only reached " << numberDecays212PbInCytoplasm1hTo26h_counter << " number of decays at entry number " << myReaderCytoplasmSim.GetCurrentEntry() << std::endl;
            }
        }

        // Create a new ROOT file
        std::string fileName_i = "FileIteration" + std::to_string(itNum) + ".root";
        TFile *FileIteration_i = new TFile(fileName_i.c_str(), "RECREATE");

        // Create a TTree named "tree" in the file
        std::string treeName_i = "CellHitTree_" + std::to_string(itNum);
        TTree *tree_i = new TTree(treeName_i.c_str(), "Data Cell Hits Tree");

        // Creating tree to store cell hit info
        int cellID;
        double doseMem, doseCyt, doseNuc;
        double doseMem_FromSolution, doseMem_FromMembrane, doseMem_FromCytoplasm;
        double doseCyt_FromSolution, doseCyt_FromMembrane, doseCyt_FromCytoplasm;
        double doseNuc_FromSolution, doseNuc_FromMembrane, doseNuc_FromCytoplasm;
        int NHitsAlpha_TotalCell, NHitsAlpha_Membrane, NHitsAlpha_Cytoplasm, NHitsAlpha_Nucleus;

        // Create a branch for each column with the respective variable
        tree_i->Branch("CellID", &cellID, "cellID/I");

        tree_i->Branch("DoseMembrane", &doseMem, "doseMem/D");
        tree_i->Branch("DoseCytoplasm", &doseCyt, "doseCyt/D");
        tree_i->Branch("DoseNucleus", &doseNuc, "doseNuc/D");

        tree_i->Branch("DoseMem_FromSolution", &doseMem_FromSolution, "doseMem_FromSolution/D");
        tree_i->Branch("doseMem_FromMembrane", &doseMem_FromMembrane, "doseMem_FromMembrane/D");
        tree_i->Branch("doseMem_FromCytoplasm", &doseMem_FromCytoplasm, "doseMem_FromCytoplasm/D");

        tree_i->Branch("DoseCyt_FromSolution", &doseCyt_FromSolution, "DoseCyt_FromSolution/D");
        tree_i->Branch("DoseCyt_FromMembrane", &doseCyt_FromMembrane, "DoseCyt_FromMembrane/D");
        tree_i->Branch("DoseCyt_FromCytoplasm", &doseCyt_FromCytoplasm, "DoseCyt_FromCytoplasm/D");

        tree_i->Branch("DoseNuc_FromSolution", &doseNuc_FromSolution, "DoseNuc_FromSolution/D");
        tree_i->Branch("DoseNuc_FromMembrane", &doseNuc_FromMembrane, "DoseNuc_FromMembrane/D");
        tree_i->Branch("DoseNuc_FromCytoplasm", &doseNuc_FromCytoplasm, "DoseNuc_FromCytoplasm/D");

        tree_i->Branch("NumHitsAlpha_TotalCell", &NHitsAlpha_TotalCell, "NHitsAlpha_TotalCell/I");
        tree_i->Branch("NHitsAlpha_Membrane", &NHitsAlpha_Membrane, "NHitsAlpha_Membrane/I");
        tree_i->Branch("NHitsAlpha_Cytoplasm", &NHitsAlpha_Cytoplasm, "NHitsAlpha_Cytoplasm/I");
        tree_i->Branch("NHitsAlpha_Nucleus", &NHitsAlpha_Nucleus, "NHitsAlpha_Nucleus/I");


        //------------------–----------
        // Looping over all stored cell energy depositions
        for(int i=0; i<storedCellHits.size(); i++)
        {
            //-------------------------------
            // Finalixing information on cell hit
            storedCellHits[i].FinalizeCellHit();

            //---------------------------
            // Sorting information into tree
            cellID = storedCellHits[i].GetCellID();

            doseMem = storedCellHits[i].GetEnergyDepositionMembrane()/storedCellHits[i].GetMassMembrane();
            doseCyt = storedCellHits[i].GetEnergyDepositionCytoplasm()/storedCellHits[i].GetMassCytoplasm();
            doseNuc = storedCellHits[i].GetEnergyDepositionNucleus()/storedCellHits[i].GetMassNucleus();

            doseMem_FromSolution = storedCellHits[i].GetEnergyDepositionMembrane_FromSolution()/storedCellHits[i].GetMassMembrane();
            doseMem_FromMembrane = storedCellHits[i].GetEnergyDepositionMembrane_FromMembrane()/storedCellHits[i].GetMassMembrane();
            doseMem_FromCytoplasm = storedCellHits[i].GetEnergyDepositionMembrane_FromCytoplasm()/storedCellHits[i].GetMassMembrane();

            doseCyt_FromSolution = storedCellHits[i].GetEnergyDepositionCytoplasm_FromSolution()/storedCellHits[i].GetMassCytoplasm();
            doseCyt_FromMembrane = storedCellHits[i].GetEnergyDepositionCytoplasm_FromMembrane()/storedCellHits[i].GetMassCytoplasm();
            doseCyt_FromCytoplasm = storedCellHits[i].GetEnergyDepositionCytoplasm_FromCytoplasm()/storedCellHits[i].GetMassCytoplasm();

            doseNuc_FromSolution = storedCellHits[i].GetEnergyDepositionNucleus_FromSolution()/storedCellHits[i].GetMassNucleus();
            doseNuc_FromMembrane = storedCellHits[i].GetEnergyDepositionNucleus_FromMembrane()/storedCellHits[i].GetMassNucleus();
            doseNuc_FromCytoplasm = storedCellHits[i].GetEnergyDepositionNucleus_FromCytoplasm()/storedCellHits[i].GetMassNucleus();

            NHitsAlpha_TotalCell = storedCellHits[i].GetNumberHitsAlphas_TotalCell();
            NHitsAlpha_Membrane = storedCellHits[i].GetNumberHitsAlphas_Membrane();
            NHitsAlpha_Cytoplasm = storedCellHits[i].GetNumberHitsAlphas_Cytoplasm();
            NHitsAlpha_Nucleus = storedCellHits[i].GetNumberHitsAlphas_Nucleus();

            //--------------------------
            // Filling tree
            tree_i->Fill();
        }

        // Write the tree into the ROOT file
        tree_i->Write();

        // Close the file
        FileIteration_i->Close();

        promises[itNum-1].set_value(FileIteration_i);

    };

    //----------------------------
    std::string filepathSolutionIteration_i;
    std::string filepathMembraneIteration_i;
    std::string filepathCytoplasmIteration_i;

    //------------------–----------
    std::string filepathSimulationOutput = "../../GEANT4Simulations/OutputCellDamageSimulation/";

    std::vector<TFile*> inputFiles_solution;
    std::vector<TFile*> inputFiles_membrane;
    std::vector<TFile*> inputFiles_cytoplasm;

    for(int i=0; i<numberIterations; i++)
    {
        std::string filepathSolutionIteration_i = filepathSimulationOutput + "Output_Pb212_" + decayDynamicsInstance.GetCellLine() + "_Activity" + std::to_string(decayDynamicsInstance.GetActivity()) + "kBq_Solution_Thread_" + std::to_string(i) + ".root";
        std::string filepathMembraneIteration_i = filepathSimulationOutput + "Output_Pb212_" + decayDynamicsInstance.GetCellLine() + "_Activity" + std::to_string(decayDynamicsInstance.GetActivity()) + "kBq_Membrane_Thread_" + std::to_string(i) + ".root";
        std::string filepathCytoplasmIteration_i = filepathSimulationOutput + "Output_Pb212_" + decayDynamicsInstance.GetCellLine() + "_Activity" + std::to_string(decayDynamicsInstance.GetActivity()) + "kBq_Cytoplasm_Thread_" + std::to_string(i) + ".root";

        inputFiles_solution.push_back(new TFile(filepathSolutionIteration_i.c_str(), "READ"));
        inputFiles_membrane.push_back(new TFile(filepathMembraneIteration_i.c_str(), "READ"));
        inputFiles_cytoplasm.push_back(new TFile(filepathCytoplasmIteration_i.c_str(), "READ"));
    }


    //------------------–----------
    // Filling histgrams
    for(int i=0; i<numberIterations; i++)
    {
        threads.push_back(std::thread(MakeFileCellHitsOneIteration, inputFiles_solution[i], inputFiles_membrane[i], inputFiles_cytoplasm[i], i+1));
    }

    // Join threads and handle exceptions.
    for (int i=0; i < numberIterations; ++i) {
        try {
            if (threads[i].joinable()) {
                std::cout << "Joining thread " << i << std::endl;
                threads[i].join();
            }
        } catch (const std::system_error& e) {
            std::cerr << "System error while joining thread " << i << ": " << e.what() << std::endl;
        } catch (const std::exception& e) {
            std::cerr << "Exception while joining thread " << i << ": " << e.what() << std::endl;
        } catch (...) {
            std::cerr << "Unknown exception while joining thread " << i << std::endl;
        }
    }

    //---------------------------
    // Create TFile to store average values, i.e. average over all iterations
    std::string fileName_main = "FileIteration" + std::to_string(itNum) + ".root";
    TFile *FileIteration_i = new TFile(fileName_i.c_str(), "RECREATE");

    // Create a TTree named "tree" in the file
    std::string treeName_i = "CellHitTree_" + std::to_string(itNum);
    TTree *tree_i = new TTree(treeName_i.c_str(), "Data Cell Hits Tree");

    // Creating tree to store cell hit info
    int cellID;
    double doseMem, doseCyt, doseNuc;
    double doseMem_FromSolution, doseMem_FromMembrane, doseMem_FromCytoplasm;
    double doseCyt_FromSolution, doseCyt_FromMembrane, doseCyt_FromCytoplasm;
    double doseNuc_FromSolution, doseNuc_FromMembrane, doseNuc_FromCytoplasm;
    int NHitsAlpha_TotalCell, NHitsAlpha_Membrane, NHitsAlpha_Cytoplasm, NHitsAlpha_Nucleus;

    // Create a branch for each column with the respective variable
    tree_i->Branch("CellID", &cellID, "cellID/I");

    tree_i->Branch("DoseMembrane", &doseMem, "doseMem/D");
    tree_i->Branch("DoseCytoplasm", &doseCyt, "doseCyt/D");
    tree_i->Branch("DoseNucleus", &doseNuc, "doseNuc/D");

    tree_i->Branch("DoseMem_FromSolution", &doseMem_FromSolution, "doseMem_FromSolution/D");
    tree_i->Branch("doseMem_FromMembrane", &doseMem_FromMembrane, "doseMem_FromMembrane/D");
    tree_i->Branch("doseMem_FromCytoplasm", &doseMem_FromCytoplasm, "doseMem_FromCytoplasm/D");

    tree_i->Branch("DoseCyt_FromSolution", &doseCyt_FromSolution, "DoseCyt_FromSolution/D");
    tree_i->Branch("DoseCyt_FromMembrane", &doseCyt_FromMembrane, "DoseCyt_FromMembrane/D");
    tree_i->Branch("DoseCyt_FromCytoplasm", &doseCyt_FromCytoplasm, "DoseCyt_FromCytoplasm/D");

    tree_i->Branch("DoseNuc_FromSolution", &doseNuc_FromSolution, "DoseNuc_FromSolution/D");
    tree_i->Branch("DoseNuc_FromMembrane", &doseNuc_FromMembrane, "DoseNuc_FromMembrane/D");
    tree_i->Branch("DoseNuc_FromCytoplasm", &doseNuc_FromCytoplasm, "DoseNuc_FromCytoplasm/D");

    tree_i->Branch("NumHitsAlpha_TotalCell", &NHitsAlpha_TotalCell, "NHitsAlpha_TotalCell/I");
    tree_i->Branch("NHitsAlpha_Membrane", &NHitsAlpha_Membrane, "NHitsAlpha_Membrane/I");
    tree_i->Branch("NHitsAlpha_Cytoplasm", &NHitsAlpha_Cytoplasm, "NHitsAlpha_Cytoplasm/I");
    tree_i->Branch("NHitsAlpha_Nucleus", &NHitsAlpha_Nucleus, "NHitsAlpha_Nucleus/I");


    // Retrieve and print the results from the futures
    for (auto &future : futures) {
        // addHistograms.AddHistograms(histMain, future.get());
    }

    // //------------------–----------
    // // Scaling histograms
    // double scalingFactor = decayDynamicsInstance.GetNumberCells()*((double) numberIterations);
    // histMain.ScaleHistograms(1./scalingFactor);


    //------------------–----------
    return 2;

}

int main(int argc, char *argv[])
{
    //------------------–----------
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " string_argument integer_argument" << std::endl;
        return 1;
    }

    //------------------–----------
    std::string cellLine = argv[1];
    int activity;
    int numberIterations;

    //------------------–----------
    std::vector<int> validActivities = {1,3,5,10,25,50,75,100,150};

    //--------------------------
    try {
        activity = std::stoi(argv[2]);
        bool activityIsValid = false;
        for(auto & entry : validActivities){
            if(entry == activity)
            {
                activityIsValid = true;
            }
        }
        if(!activityIsValid)
        {
            throw std::invalid_argument("Activity case not valid.");
        }
    } catch (const std::invalid_argument& e) {
        std::cerr << "Error: second argument is not a valid integer! " << e.what() << std::endl;
        return 2;
    } catch (const std::out_of_range& e) {
        std::cerr << "Error: second argument is out of range for an integer! " << e.what() << std::endl;
        return 3;
    }

    //--------------------------
    try {
        numberIterations = std::stoi(argv[3]);
        if(numberIterations<=0){
            throw std::invalid_argument("Number of iterations cannot be smaller than one.");
        }
    } catch (const std::invalid_argument& e) {
        std::cerr << "Error: third argument is not a valid integer! " << e.what() << std::endl;
        return 4;
    } catch (const std::out_of_range& e) {
        std::cerr << "Error: third argument is out of range for an integer! " << e.what() << std::endl;
        return 5;
    }

    std::cout << "Cell line: " << cellLine << ", Activity: " << activity <<  ", Number Iterations: " << numberIterations  << std::endl;


    //------------------–----------
    ROOT::EnableThreadSafety();
    TH1::AddDirectory(false);

    //------------------–----------
    // Defining decay dynamics

    DecayDynamics decays = DecayDynamics(activity, cellLine);

    //------------------–----------
    // Loading decay dynamics calculations

    std::string mathematicaOutput = "../../Mathematica/Output";
    decays.LoadDataFromMathematicaCalculations(mathematicaOutput.c_str());

    int o = AnalyzeHistogramsFromSimulation(decays, numberIterations);
    std::cout << o << std::endl;

    // EnergyDepositionHistograms hists = AnalyzeHistogramsFromSimulation(decays, numberIterations);
    // std::string outputName = "Output_" + cellLine + "_" + std::to_string(activity) + "kBq.root";
    // auto output = new TFile(outputName.c_str(), "RECREATE");
    // hists.WriteHistogramsToFile();
    // output->Write();
    // output->Close();

    return 0;
}