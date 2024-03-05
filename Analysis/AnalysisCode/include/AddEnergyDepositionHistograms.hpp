#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
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
#include <future>
#include <thread>

#ifndef ADDENERGYDEPOSITIONHISTOGRAMS_HPP
#define ADDENERGYDEPOSITIONHISTOGRAMS_HPP


class AddEnergyDepositionHistograms{
    public:
        void AddHistograms(EnergyDepositionHistograms& h1, const EnergyDepositionHistograms& h2);
};


#endif // ADDENERGYDEPOSITIONHISTOGRAMS_HPP