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


#ifndef SURVIVALFIT_HPP
#define SURVIVALFIT_HPP


class CellSurvival;

//----------------------------
class SurvivalFit
{
    public:
        SurvivalFit();
        void LoadHistogramsFromAnalysis(CellSurvival cellSurvivalInstance, std::string regionName);
        void FitCellSurvival(CellSurvival cellSurvivalInstance, std::string modelName_in, int region);

        void MakeHitMultiplicitySurvivalHistogram();
        void WriteToFile(TFile *file);

        std::vector<std::tuple<double, TH1D*, TH1D*>> Get_hDose_Activity_Vec(){return hDose_Activity_Vec;};
        std::vector<std::tuple<double, TH1D*>> Get_hHitMultiplicity_vec(){return hHitMultiplicity_vec;};
        std::vector<std::tuple<double, TH2D*>> Get_hDose_HitsCellComponent_Vec(){return hDose_HitsCellComponent_Vec;};

        std::vector<std::tuple<double,double>> Get_ParametersAndUncertainties_Vec(){return parametersAndUncertainties_Vec;};
        std::string Get_RegionName(){return regionName;};

    private:
        std::string regionName;

        std::string cellLine;

        std::string modelName;

        std::string cellGeometryType;


        //-------------------------
        std::vector<std::tuple<double,double>> parametersAndUncertainties_Vec;


        //-------------------------------
        std::vector<std::tuple<double, TH1D*, TH1D*>> hDose_Activity_Vec;

        std::vector<std::tuple<double, TH2D*>> hDose_HitsCellComponent_Vec;

        //--------------------------
        std::vector<std::tuple<double, TH1D*>> hHitMultiplicity_vec;

        std::vector<std::tuple<double, TH1D*>> hHitMultiplicitySurvival_vec;


        //-----------------------
        TGraphErrors *gr_clonogenicSurvival;

        TGraphErrors gr_cellSurvivability_vs_activitykBqPerMl;

        TF1* f_cellSurvivalVsDose;

        TGraphErrors* grFitPoints = new TGraphErrors();
};


#endif // SURVIVALFIT_HPP
