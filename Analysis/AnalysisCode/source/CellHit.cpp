#include "../include/CellHit.hpp"

//------------------–----------
CellHit::CellHit(int cellID_in)
{
    cellID = cellID_in;
    numberHitsAlphas_Nucleus = 0;
    numberHitsAlphas_Membrane = 0;
    numberHitsAlphas_Cytoplasm = 0;
    numberHitsAlphas_TotalCell = 0;

    numberHitsAlphasTotalCell_FromSolution = 0;
    numberHitsAlphasTotalCell_FromMembrane = 0;
    numberHitsAlphasTotalCell_FromCytoplasm = 0;

    numberHitsAlphasNucleus_FromSolution = 0;
    numberHitsAlphasNucleus_FromMembrane = 0;
    numberHitsAlphasNucleus_FromCytoplasm = 0;

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
    // If membrane hit
    if(volumeTypeHit==1)
    {
        numberHitsAlphas_Membrane++;
    }
    // If cytoplasm hit
    if(volumeTypeHit==2)
    {
        numberHitsAlphas_Cytoplasm++;
    }
    // If nucleus hit
    if(volumeTypeHit==3)
    {
        numberHitsAlphas_Nucleus++;

        // If decay originated in solution
        if(volumeTypeOriginDecay_in==0)
        {
            kineticEnergyAlphaNucleus_FromSolution_Vec.push_back(kineticEnergyAlpha);
            numberHitsAlphasNucleus_FromSolution++;
        }

        // If decay originated in membrane
        if(volumeTypeOriginDecay_in==1)
        {
            kineticEnergyAlphaNucleus_FromMembrane_Vec.push_back(kineticEnergyAlpha);
            numberHitsAlphasNucleus_FromMembrane++;
        }

        // If decay originated in cytoplasm
        if(volumeTypeOriginDecay_in==2)
        {
            kineticEnergyAlphaNucleus_FromCytoplasm_Vec.push_back(kineticEnergyAlpha);
            numberHitsAlphasNucleus_FromCytoplasm++;
        }
    }

    // If it's the first time the alpha is registered in the cell
    if(firstTimeCountingAlpha)
    {
        // Count alpha hit for whole cell
        numberHitsAlphas_TotalCell++;

        // If decay originated in solution
        if(volumeTypeOriginDecay_in==0)
        {
            kineticEnergyAlphaTotalCell_FromSolution_Vec.push_back(kineticEnergyAlpha);
            // std::cout << "Filled alpha hit cell in solution" << std::endl;
            numberHitsAlphasTotalCell_FromSolution++;
        }

        // If decay originated in membrane
        if(volumeTypeOriginDecay_in==1)
        {
            kineticEnergyAlphaTotalCell_FromMembrane_Vec.push_back(kineticEnergyAlpha);
            numberHitsAlphasTotalCell_FromMembrane++;
        }

        // If decay originated in cytoplasm
        if(volumeTypeOriginDecay_in==2)
        {
            kineticEnergyAlphaTotalCell_FromCytoplasm_Vec.push_back(kineticEnergyAlpha);
            numberHitsAlphasTotalCell_FromCytoplasm++;
        }
    }
}

//------------------–----------
void CellHit::FinalizeCellHit(std::string cellGeometry)
{
    double densityWater = 1000. ; // kg/m^3
    double radiusCell = 9.0e-6; // m
    double radiusCytoplasm = radiusCell - 4.0e-9; // m

    double radiusNucleus;

    if(cellGeometry=="D12RP"||cellGeometry=="D12CP")
    {
        radiusNucleus = 6.0e-6; // m
    }

    if(cellGeometry=="D5RP"||cellGeometry=="D5CP")
    {
        radiusNucleus = 2.5e-6; // m
    }

    massNucleus = (4./3.)*TMath::Pi()*std::pow(radiusNucleus,3.)*densityWater; // kg
    massCytoplasm = (4./3.)*TMath::Pi()*std::pow(radiusCytoplasm,3.)*densityWater - massNucleus; // kg
    massCell = (4./3.)*TMath::Pi()*std::pow(radiusCell,3.)*densityWater; // kg
    massMembrane = massCell - (4./3.)*TMath::Pi()*std::pow(radiusCytoplasm,3.)*densityWater; // kg

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

        // If in membrane
        if(interactionVolume==1)
        {
            energyDepMembrane += energyDepInteraction;
            if(originVolume == 0){energyDepMembrane_FromSolution += energyDepInteraction;}
            if(originVolume == 1){energyDepMembrane_FromMembrane += energyDepInteraction;}
            if(originVolume == 2){energyDepMembrane_FromCytoplasm += energyDepInteraction;}
        }
        // If in cytoplasm
        if(interactionVolume==2)
        {
            energyDepCytoplasm += energyDepInteraction;
            if(originVolume == 0){energyDepCytoplasm_FromSolution += energyDepInteraction;}
            if(originVolume == 1){energyDepCytoplasm_FromMembrane += energyDepInteraction;}
            if(originVolume == 2){energyDepCytoplasm_FromCytoplasm += energyDepInteraction;}
        }
        // If in nucleus
        if(interactionVolume==3)
        {
            energyDepNucleus += energyDepInteraction;
            if(originVolume == 0){energyDepNucleus_FromSolution += energyDepInteraction;}
            if(originVolume == 1){energyDepNucleus_FromMembrane += energyDepInteraction;}
            if(originVolume == 2){energyDepNucleus_FromCytoplasm += energyDepInteraction;}
        }

    }

    //-----------------------------
    // Finding fractions for origin of decay
    fractionEnergyDepMembrane_FromSolution = energyDepMembrane_FromSolution/energyDepMembrane;
    fractionEnergyDepMembrane_FromMembrane = energyDepMembrane_FromMembrane/energyDepMembrane;
    fractionEnergyDepMembrane_FromCytoplasm = energyDepMembrane_FromCytoplasm/energyDepMembrane;

    fractionEnergyDepCytoplasm_FromSolution = energyDepCytoplasm_FromSolution/energyDepCytoplasm;
    fractionEnergyDepCytoplasm_FromMembrane = energyDepCytoplasm_FromMembrane/energyDepCytoplasm;
    fractionEnergyDepCytoplasm_FromCytoplasm = energyDepCytoplasm_FromCytoplasm/energyDepCytoplasm;

    fractionEnergyDepNucleus_FromSolution = energyDepNucleus_FromSolution/energyDepNucleus;
    fractionEnergyDepNucleus_FromMembrane = energyDepNucleus_FromMembrane/energyDepNucleus;
    fractionEnergyDepNucleus_FromCytoplasm = energyDepNucleus_FromCytoplasm/energyDepNucleus;

    energyDepTotalCell = energyDepMembrane + energyDepCytoplasm + energyDepNucleus;

    energyDepTotalCell_FromSolution = energyDepMembrane_FromSolution + energyDepCytoplasm_FromSolution + energyDepNucleus_FromSolution;
    energyDepTotalCell_FromMembrane = energyDepMembrane_FromMembrane + energyDepCytoplasm_FromMembrane + energyDepNucleus_FromMembrane;
    energyDepTotalCell_FromCytoplasm = energyDepMembrane_FromCytoplasm + energyDepCytoplasm_FromCytoplasm + energyDepNucleus_FromCytoplasm;

    fractionEnergyDepTotalCell_FromSolution = energyDepTotalCell_FromSolution/energyDepTotalCell;
    fractionEnergyDepTotalCell_FromMembrane = energyDepTotalCell_FromMembrane/energyDepTotalCell;
    fractionEnergyDepTotalCell_FromCytoplasm = energyDepTotalCell_FromCytoplasm/energyDepTotalCell;


}
