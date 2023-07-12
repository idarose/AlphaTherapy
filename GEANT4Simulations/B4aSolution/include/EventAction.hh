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

    // Functions for extracting vectores with interaction information
    std::vector<G4double>& GetEnergyDepVec() {return energyDepVec;};
    std::vector<G4int>& GetCellIDVec() {return cellIDVec;};
    std::vector<G4int>& GetVolumeTypeVec() {return volumeTypeVec;};
    std::vector<G4double>& GetKineticEnergyVec() {return kineticEnergyVec;};
    std::vector<G4int>& GetParticleTypeVec() {return particleTypeVec;};
    std::vector<G4double>& GetInteractionTime() {return interactionTimeVec;};




    // Methods for determining if first interaction took place in the solution:

    // Function for extracting the value of the boolean
    bool& GetBooleanFirstInteractionInSolution() {return booleanFirstInteractionInSolution;};

    // Function for updating the value of the boolean
    void UpdateBooleanFirstInteractionInSolution(bool newValue) {booleanFirstInteractionInSolution = newValue;};

  private:
    std::vector<G4double> energyDepVec;
    std::vector<G4int> cellIDVec;
    std::vector<G4int> volumeTypeVec;
    std::vector<G4double> kineticEnergyVec;
    std::vector<G4int> particleTypeVec;
    std::vector<G4double> interactionTimeVec;


    // Boolean. True: First interaction took place in solution. False: First interaction did not take place in solution
    bool booleanFirstInteractionInSolution;
};




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


