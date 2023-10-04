#ifndef B4aEventAction_h
#define B4aEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

#include <vector>

class PrimaryGeneratorAction;


class EventAction : public G4UserEventAction
{
  public:
    EventAction(PrimaryGeneratorAction* primaryGeneratorAction);
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
    std::vector<G4double>& GetFirstInteractionTime() {return firstInteractionTimeVec;};
    std::vector<G4int>& GetFirstInteractionVolume() {return firstInteractionVolumeVec;};



    // Methods to keep track of step number
    void AddStepNumber() {stepNumber++;};
    G4int GetStepNumber() {return stepNumber;};
    void ResetStepNumber() {stepNumber=0;};

  private:
    std::vector<G4double> energyDepVec;
    std::vector<G4int> cellIDVec;
    std::vector<G4int> volumeTypeVec;
    std::vector<G4double> kineticEnergyVec;
    std::vector<G4int> particleTypeVec;
    std::vector<G4double> interactionTimeVec;

    std::vector<G4double> firstInteractionTimeVec;
    std::vector<G4int> firstInteractionVolumeVec;


    PrimaryGeneratorAction* fPrimaryGeneratorAction;


    G4int stepNumber = 0;

};




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


