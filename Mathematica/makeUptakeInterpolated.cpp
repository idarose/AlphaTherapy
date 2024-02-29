#include <iostream>
#include <vector>
#include <tuple>
#include <Math/Interpolator.h>

void makeUptakeInterpolated()
{

    auto WriteInterpolatedUptakeToFile = [&](std::vector<std::tuple<double,double>> U0TotalAbsorbed, std::vector<std::tuple<double,double>> U0AbsorbedCytoplasm, std::string cellLine)
    {
        std::vector<double> xVals;
        std::vector<double> yVals_total;
        std::vector<double> yVals_cytoplasm;

        for(int i=0; i<U0TotalAbsorbed.size(); i++)
        {
            double U0Cytoplasm = std::get<1>(U0AbsorbedCytoplasm[i]);
            double U0Total = std::get<1>(U0TotalAbsorbed[i]);
            double activity = std::get<0>(U0AbsorbedCytoplasm[i]);

            xVals.push_back(activity);
            yVals_total.push_back(U0Total);
            yVals_cytoplasm.push_back(U0Cytoplasm);
        }

        ROOT::Math::Interpolator interpolator_total(xVals, yVals_total, ROOT::Math::Interpolation::kLINEAR);
        ROOT::Math::Interpolator interpolator_cytoplasm(xVals, yVals_cytoplasm, ROOT::Math::Interpolation::kLINEAR);

        std::ofstream file;

        std::string fileName = "InterpolatedUptake_" + cellLine + ".csv";
        file.open(fileName.c_str(), std::ios::out | std::ios::trunc);

        std::vector<double> activities_needed = {1., 3., 5., 10., 25., 50., 75., 100., 150.};

        if(file.is_open()) {

            // write the column names
            file << "Activity,U0InternalizedTotal,U0Cytoplasm" << std::endl;

            for(auto & entry : activities_needed)
            {
                double interpolated_U0Total = interpolator_total.Eval(entry);
                double interpolated_U0Cytoplasm = interpolator_cytoplasm.Eval(entry);

                file << entry  << " " << interpolated_U0Total << " " << interpolated_U0Cytoplasm << std::endl;
            }
        }

        file.close();
    };

    std::vector<std::tuple<double,double>> U0InternalizedTotal_C4_2;
    std::vector<std::tuple<double,double>> U0Cytoplasm_C4_2;

    U0InternalizedTotal_C4_2.push_back(std::make_tuple(0.,0.));
    U0InternalizedTotal_C4_2.push_back(std::make_tuple(5.,2.31));
    U0InternalizedTotal_C4_2.push_back(std::make_tuple(10.,3.95));
    U0InternalizedTotal_C4_2.push_back(std::make_tuple(25.,10.69));
    U0InternalizedTotal_C4_2.push_back(std::make_tuple(50.,20.02));
    U0InternalizedTotal_C4_2.push_back(std::make_tuple(75.,23.92));
    U0InternalizedTotal_C4_2.push_back(std::make_tuple(100.,35.89));
    U0InternalizedTotal_C4_2.push_back(std::make_tuple(150.,39.96));

    U0Cytoplasm_C4_2.push_back(std::make_tuple(0.,0.));
    U0Cytoplasm_C4_2.push_back(std::make_tuple(5.,1.14));
    U0Cytoplasm_C4_2.push_back(std::make_tuple(10.,1.97));
    U0Cytoplasm_C4_2.push_back(std::make_tuple(25.,4.78));
    U0Cytoplasm_C4_2.push_back(std::make_tuple(50.,8.94));
    U0Cytoplasm_C4_2.push_back(std::make_tuple(75.,10.79));
    U0Cytoplasm_C4_2.push_back(std::make_tuple(100.,13.16));
    U0Cytoplasm_C4_2.push_back(std::make_tuple(150.,16.40));

    WriteInterpolatedUptakeToFile(U0InternalizedTotal_C4_2, U0Cytoplasm_C4_2, "C4_2");

    std::vector<std::tuple<double,double>> U0InternalizedTotal_PC3_PIP;
    std::vector<std::tuple<double,double>> U0Cytoplasm_PC3_PIP;

    U0InternalizedTotal_PC3_PIP.push_back(std::make_tuple(0.,0.));
    U0InternalizedTotal_PC3_PIP.push_back(std::make_tuple(10.,47.));
    U0InternalizedTotal_PC3_PIP.push_back(std::make_tuple(25.,119.));
    U0InternalizedTotal_PC3_PIP.push_back(std::make_tuple(50.,229.));
    U0InternalizedTotal_PC3_PIP.push_back(std::make_tuple(75.,335.));
    U0InternalizedTotal_PC3_PIP.push_back(std::make_tuple(100.,448.));
    U0InternalizedTotal_PC3_PIP.push_back(std::make_tuple(150.,565.));

    U0Cytoplasm_PC3_PIP.push_back(std::make_tuple(0.,0.));
    U0Cytoplasm_PC3_PIP.push_back(std::make_tuple(10.,2.));
    U0Cytoplasm_PC3_PIP.push_back(std::make_tuple(25.,6.));
    U0Cytoplasm_PC3_PIP.push_back(std::make_tuple(50.,11.));
    U0Cytoplasm_PC3_PIP.push_back(std::make_tuple(75.,17.));
    U0Cytoplasm_PC3_PIP.push_back(std::make_tuple(100.,22.));
    U0Cytoplasm_PC3_PIP.push_back(std::make_tuple(150.,28.));

    WriteInterpolatedUptakeToFile(U0InternalizedTotal_PC3_PIP, U0Cytoplasm_PC3_PIP, "PC3_PIP");

}