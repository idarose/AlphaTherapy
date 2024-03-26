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

#ifndef HITANALYSIS_HPP
#define HITANALYSIS_HPP

class SurvivalFit;

//----------------------------
class HitAnalysis
{
    public:
        HitAnalysis(SurvivalFit survivalFitInstance);

        void MakeHitMulitplicityGraphs();
        void MakeHitMultiplicity_SurvivalFraction_Histograms();
        void MakeHitMultiplicity_SurvivalFraction_Graphs();
        void MakeHitMultiplicity_PercentKilled_Histograms();
        void MakeHitMultiplicity_PercentKilled_Graphs();
        void MakeGraph_ProbabilityDeath_PerNHits(int MaxNHits);
        void MakeGraph_ProbabilityDeath_UWA_ForNHits();
        void MakeHitAnalysis(int MaxNHits);
        void WriteToFile();

    private:
        std::string regionName;

        //-----------------------------
        // Define number of bins and max number of hits
        int nBins;
        double maxDose;

        double scalingFactorHistograms;

        std::vector<std::tuple<double,double>> parametersFit;

        //-------------------------------
        std::vector<std::tuple<double, TH1D*, TH1D*>> hDose_Activity_Vec;

        std::vector<std::tuple<double, TH2D*>> hDose_HitsCellComponent_Vec;

        //--------------------------
        std::vector<std::tuple<double, TH1D*>> hHitMultiplicity_vec;

        std::vector<std::tuple<double, TGraph*>> grHitMultiplicity_vec;

        //------------------------------
        std::vector<std::tuple<double, TH1D*>> hHitMultiplicity_SurvivalFraction_vec;

        std::vector<std::tuple<double, TGraphErrors*>> grHitMultiplicity_SurvivalFraction_vec;

        //--------------------------
        std::vector<std::tuple<double, TH1D*>> hHitMultiplicity_PercentKilled_Vec;

        std::vector<std::tuple<double, TGraphErrors*>> grHitMultiplicity_PercentKilled_Vec;

        //-----------------------------
        std::vector<std::tuple<int, TGraphErrors*>> grPercentKilled_PerOneHit_Vec;

        TGraphErrors* grProbabilityDeath_ForNHits = new TGraphErrors();
};


#endif // HITANALYSIS_HPP