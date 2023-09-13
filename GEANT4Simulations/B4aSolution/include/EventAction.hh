#ifndef B4aEventAction_h
#define B4aEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

#include <vector>

class RunAction;

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
    void SetFirstInteractionInfo(G4double firstInteractionTime, G4int firstInteractionVolume);

    // Functions for extracting vectores with interaction information
    std::vector<G4double>& GetEnergyDepVec() {return energyDepVec;};
    std::vector<G4int>& GetCellIDVec() {return cellIDVec;};
    std::vector<G4int>& GetVolumeTypeVec() {return volumeTypeVec;};
    std::vector<G4double>& GetKineticEnergyVec() {return kineticEnergyVec;};
    std::vector<G4int>& GetParticleTypeVec() {return particleTypeVec;};
    std::vector<G4double>& GetInteractionTime() {return interactionTimeVec;};
    std::vector<G4int>& GetParticleTypeVec() {return particleTypeVec;};
    std::vector<G4double>& GetInteractionTime() {return interactionTimeVec;};


    void SetFirstInteractionVolume(G4int volumeType) {firstInteractionVolume = volumeType;};
    G4int& GetFirstInteractionVolume() {return firstInteractionVolume;};

    // Methods to add, reset and extract step number
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
    std::vector<G4double> firstInteractionTimeVec;
    std::vector<G4int> firstInteractionVolumeVec;

    // If = 1 first interaction in solution, if = 0 first interaction not in solution
    G4int firstInteractionVolume;

    // Number to count steps
    int stepNumber = 0;

};




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


