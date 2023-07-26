#ifndef mainAnalysisCode_h
#define mainAnalysisCode_h 1

class cellHit
{
    public:
        CellHit(int cellID_in);

        void AddEnergyDeposition(double energyDep_in, double volumeType_in);
        void AddInteractionTime(double interactionTime_in);

        void GetSumEnergyDepositions(){return energyDepMembrane + energyDepCytoplasm + energyDepNucleus;};

    private:
        int cellID;
        int volumeType
        double interactionTime;

        double energyDepMembrane;
        double energyDepCytoplasm;
        double energyDepNucleus;

        double sumEnergyDepositions;
};

#endif