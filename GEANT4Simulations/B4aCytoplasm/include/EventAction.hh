#ifndef B4aEventAction_h
#define B4aEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

#include <vector>


class EventAction : public G4UserEventAction
{
  public:
    EventAction();
    ~EventAction() override;

    void  BeginOfEventAction(const G4Event* event) override;
    void    EndOfEventAction(const G4Event* event) override;

    void AddAbs(G4double de, G4double dl);
    void AddGap(G4double de, G4double dl);

    void StoreInteractionInformation(G4double energyDep, G4int cellID, G4int volumeType, G4double kineticEnergy, G4int particleType, G4double interactionTime);


    // Functions for extracting vectores with interaction information
    std::vector<G4double>& GetEnergyDepVec() {return energyDepVec;};
    std::vector<G4int>& GetCellIDVec() {return cellIDVec;};
    std::vector<G4int>& GetVolumeTypeVec() {return volumeTypeVec;};
    std::vector<G4double>& GetKineticEnergyVec() {return kineticEnergyVec;};
    std::vector<G4int>& GetParticleTypeVec() {return particleTypeVec;};
    std::vector<G4double>& GetInteractionTime() {return interactionTimeVec;};

    // Methods to make sure information of events where first interaction takes place inside cell nucleus is not stored
    void BooleanFirstInteractionInCellNucleusMakeTrue() {booleanFirstInteractionInCellNucleus = true;};
    void BooleanFirstInteractionInCellNucleusMakeFalse() {booleanFirstInteractionInCellNucleus = false;};
    bool GetBooleanFirstInteractionInCellNucleus(){return booleanFirstInteractionInCellNucleus;};

    void AddStepNumber(){stepNumber++;};
    int GetStepNumber(){return stepNumber;};
    void ResetStepNumber(){stepNumber=0;};

  private:
    std::vector<G4double> energyDepVec;
    std::vector<G4int> cellIDVec;
    std::vector<G4int> volumeTypeVec;
    std::vector<G4double> kineticEnergyVec;
    std::vector<G4int> particleTypeVec;
    std::vector<G4double> interactionTimeVec;

    int stepNumber = 0;


    // Boolean for first interaction information
    bool booleanFirstInteractionInCellNucleus;
};




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


