#include "TFile.h"
#include "TCanvas.h"
#include "TStopwatch.h"
#include "TGraphAsymmErrors.h"
#include "TCutG.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH2Poly.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TArrow.h"
#include "TLatex.h"
#include "TF1Convolution.h"
#include "TFitResult.h"
#include "TParameter.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TMathText.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TAxis.h"
#include "TGaxis.h"
#include <iomanip>
#include <sstream>
#include <vector>
#include <tuple>

#ifndef DOSEANALYSIS_HPP
#define DOSEANALYSIS_HPP

class SurvivalFit;

//----------------------------
class DoseAnalysis
{
    public:
        DoseAnalysis(SurvivalFit survivalFitInstance);
        void MakeMeanDose_PerNHits_Graphs();
        void MakeDose_PerNHits_Average_Histograms();
        void MakeMeanDose_PerNHits_Average_Graph();
        void MakeDoseAnalysis();
        void WriteToFile(TFile *file);

    private:
        std::string regionName;

        //-----------------------------
        // Define number of bins and max number of hits
        int nBinsHits;
        int nBinsDose;
        double maxDose;

        double scalingFactorHistograms;

        std::vector<std::tuple<double,double>> parametersFit;

        //-------------------------------
        std::vector<std::tuple<double, TH1D*, TH1D*>> hDose_Activity_Vec;

        std::vector<std::tuple<double, TH2D*>> hDose_HitsCellComponent_Vec;

        std::vector<std::tuple<int, TH1D*>> hDose_HitsCellComponent_Average_Vec;

        //--------------------------
        std::vector<std::tuple<double, TH1D*>> hHitMultiplicity_vec;

        //---------------------------
        std::vector<std::tuple<double, TGraphAsymmErrors*>> grMeanDose_PerNHits_Vec;

        TGraphAsymmErrors* grMeanDose_PerNHits_Average = new TGraphAsymmErrors();
};


#endif // DOSEANALYSIS_HPP