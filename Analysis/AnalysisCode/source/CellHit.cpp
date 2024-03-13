#include "../include/CellHit.hpp"

//------------------–----------
CellHit::CellHit(int cellID_in)
{
    cellID = cellID_in;
    nummberHitsAlphas_Nucleus = 0;
    numberHitsAlphas_Membrane = 0;
    numberHitsAlphas_Cytoplasm = 0;
    numberHitsAlphas_TotalCell = 0;

    energyDepMembrane = 0.0;
    energyDepCytoplasm = 0.0;
    energyDepNucleus = 0.0;

    energyDepMembrane_FromSolution = 0.0;
    energyDepMembrane_FromMembrane = 0.0;
    energyDepMembrane_FromCytoplasm = 0.0;

    energyDepCytoplasm_FromSolution = 0.0;
    energyDepCytoplasm_FromMembrane = 0.0;
    energyDepCytoplasm_FromCytoplasm = 0.0;

    energyDepNucleus_FromSolution = 0.0;
    energyDepNucleus_FromMembrane = 0.0;
    energyDepNucleus_FromCytoplasm = 0.0;

    energyDepTotalCell_FromSolution = 0.0;
    energyDepTotalCell_FromMembrane = 0.0;
    energyDepTotalCell_FromCytoplasm = 0.0;

    // energyDepNucleus_FromAlpha = 0.;
    // energyDepMembrane_FromAlpha = 0.;
    // energyDepCytoplasm_FromAlpha = 0.;
    // energyDepTotalCell_FromAlpha = 0.;

    // numberHitsAlphas = 0;
    // numberHitsBetas = 0;

    double densityWater = 1000. ; // kg/m^3
    double radiusCell = 9.0e-6; // m
    double radiusNucleus = 6.0e-6; // m
    double radiusCytoplasm = radiusCell - 4.0e-9; // m

    massNucleus = (4./3.)*TMath::Pi()*std::pow(radiusNucleus,3.)*densityWater; // kg
    massCytoplasm = (4./3.)*TMath::Pi()*std::pow(radiusCytoplasm,3.)*densityWater - massNucleus; // kg
    massCell = (4./3.)*TMath::Pi()*std::pow(radiusCell,3.)*densityWater; // kg
    massMembrane = massCell - (4./3.)*TMath::Pi()*std::pow(radiusCytoplasm,3.)*densityWater; // kg
}


//------------------–----------
void CellHit::AddEnergyDeposition(double energyDep_in, int volumeTypeInteraction_in, int volumeTypeOriginDecay_in, int particleType_in)
{
    std::tuple<int,double,int,int> hitTuple;
    hitTuple = std::make_tuple(volumeTypeInteraction_in, energyDep_in, volumeTypeOriginDecay_in, particleType_in);
    energyDepsVec.push_back(hitTuple);
}


//------------------–----------
void CellHit::HitByAlphaParticle(int volumeTypeHit, bool firstTimeCountingAlpha, int volumeTypeOriginDecay_in, double kineticEnergyAlpha)
{
    if(volumeTypeHit==1)
    {
        numberHitsAlphas_Membrane++;
    }
    if(volumeTypeHit==2)
    {
        numberHitsAlphas_Cytoplasm++;
    }
    if(volumeTypeHit==3)
    {
        nummberHitsAlphas_Nucleus++;
        if(volumeTypeOriginDecay_in==0)
        {
            kineticEnergyAlphaNucleus_FromSolution_Vec.push_back(kineticEnergyAlpha);
        }
        if(volumeTypeOriginDecay_in==1)
        {
            kineticEnergyAlphaNucleus_FromMembrane_Vec.push_back(kineticEnergyAlpha);
        }
        if(volumeTypeOriginDecay_in==2)
        {
            kineticEnergyAlphaNucleus_FromCytoplasm_Vec.push_back(kineticEnergyAlpha);
        }
    }

    if(firstTimeCountingAlpha)
    {
        numberHitsAlphas_TotalCell++;
        if(volumeTypeOriginDecay_in==0)
        {
            kineticEnergyAlphaTotalCell_FromSolution_Vec.push_back(kineticEnergyAlpha);
        }
        if(volumeTypeOriginDecay_in==1)
        {
            kineticEnergyAlphaTotalCell_FromMembrane_Vec.push_back(kineticEnergyAlpha);
        }
        if(volumeTypeOriginDecay_in==2)
        {
            kineticEnergyAlphaTotalCell_FromMembrane_Vec.push_back(kineticEnergyAlpha);
        }
    }
}

// //----------------------------
// void CellHit::HitByBetaParticle()
// {
//     numberHitsBetas ++;
// }


//------------------–----------
void CellHit::FinalizeCellHit()
{

    //----------------------
    // Sum the energy depositions per cell component

    int interactionVolume;
    int originVolume;
    double energyDepInteraction;
    int particleType;

    for(int i=0; i<energyDepsVec.size(); i++)
    {
        interactionVolume = std::get<0>(energyDepsVec[i]);
        energyDepInteraction = std::get<1>(energyDepsVec[i]);
        originVolume = std::get<2>(energyDepsVec[i]);
        particleType = std::get<3>(energyDepsVec[i]);


        if(interactionVolume==1)
        {
            // std::cout << "In membrane: " << energyDepInteraction << std::endl;
            energyDepMembrane += energyDepInteraction;
            if(originVolume == 0){energyDepMembrane_FromSolution += energyDepInteraction;}
            if(originVolume == 1){energyDepMembrane_FromMembrane += energyDepInteraction;}
            if(originVolume == 2){energyDepMembrane_FromCytoplasm += energyDepInteraction;}

            // if(particleType==1000020040)
            // {
            //     energyDepMembrane_FromAlpha += energyDepInteraction;
            // }
        }
        if(interactionVolume==2)
        {
            // std::cout << "In cytoplasm : " << energyDepInteraction << std::endl;
            energyDepCytoplasm += energyDepInteraction;
            if(originVolume == 0){energyDepCytoplasm_FromSolution += energyDepInteraction;}
            if(originVolume == 1){energyDepCytoplasm_FromMembrane += energyDepInteraction;}
            if(originVolume == 2){energyDepCytoplasm_FromCytoplasm += energyDepInteraction;}

            // if(particleType==1000020040)
            // {
            //     energyDepCytoplasm_FromAlpha += energyDepInteraction;
            // }
        }
        if(interactionVolume==3)
        {
            // std::cout << "In nucleus : " << energyDepInteraction << std::endl;
            energyDepNucleus += energyDepInteraction;
            if(originVolume == 0){energyDepNucleus_FromSolution += energyDepInteraction;}
            if(originVolume == 1){energyDepNucleus_FromMembrane += energyDepInteraction;}
            if(originVolume == 2){energyDepNucleus_FromCytoplasm += energyDepInteraction;}
        }

    }

    // //-----------------------------
    // // Finding fractions for origin of decay
    // fractionEnergyDepMembrane_FromSolution = energyDepMembrane_FromSolution/energyDepMembrane;
    // fractionEnergyDepMembrane_FromMembrane = energyDepMembrane_FromMembrane/energyDepMembrane;
    // fractionEnergyDepMembrane_FromCytoplasm = energyDepMembrane_FromCytoplasm/energyDepMembrane;

    // fractionEnergyDepCytoplasm_FromSolution = energyDepCytoplasm_FromSolution/energyDepCytoplasm;
    // fractionEnergyDepCytoplasm_FromMembrane = energyDepCytoplasm_FromMembrane/energyDepCytoplasm;
    // fractionEnergyDepCytoplasm_FromCytoplasm = energyDepCytoplasm_FromCytoplasm/energyDepCytoplasm;

    // fractionEnergyDepNucleus_FromSolution = energyDepNucleus_FromSolution/energyDepNucleus;
    // fractionEnergyDepNucleus_FromMembrane = energyDepNucleus_FromMembrane/energyDepNucleus;
    // fractionEnergyDepNucleus_FromCytoplasm = energyDepNucleus_FromCytoplasm/energyDepNucleus;

    energyDepTotalCell = energyDepMembrane + energyDepCytoplasm + energyDepNucleus;

    energyDepTotalCell_FromSolution = energyDepMembrane_FromSolution + energyDepCytoplasm_FromSolution + energyDepNucleus_FromSolution;
    energyDepTotalCell_FromMembrane = energyDepMembrane_FromMembrane + energyDepCytoplasm_FromMembrane + energyDepNucleus_FromMembrane;
    energyDepTotalCell_FromCytoplasm = energyDepMembrane_FromCytoplasm + energyDepCytoplasm_FromCytoplasm + energyDepNucleus_FromCytoplasm;

    // fractionEnergyDepTotalCell_FromSolution = energyDepTotalCell_FromSolution/energyDepTotalCell;
    // fractionEnergyDepTotalCell_FromMembrane = energyDepTotalCell_FromMembrane/energyDepTotalCell;
    // fractionEnergyDepTotalCell_FromCytoplasm = energyDepTotalCell_FromCytoplasm/energyDepTotalCell;

    // energyDepTotalCell_FromAlpha = energyDepMembrane_FromAlpha + energyDepNucleus_FromAlpha + energyDepCytoplasm_FromAlpha;

}
