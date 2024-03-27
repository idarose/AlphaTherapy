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

    void AddStepNumber(){stepNumber++;};
    void ResetStepNumber(){stepNumber=0;};

    void StoreInteractionInformation(G4double energyDep, G4double positionZ, G4double stepLength);

    void AddStep(G4double step){stepLengthCumulative+=step;};
    void ResetStep(){stepLengthCumulative=0.;};

    G4double GetStep(){return stepLengthCumulative;};

    // Functions for extracting vectores with interaction information
    std::vector<G4double>& GetEnergyDepVec() {return energyDepVec;};
    std::vector<G4double>& GetPositionZVec() {return positionZVec;};
    std::vector<G4double>& GetStepLengthVec() {return stepLengthVec;};

  private:
    std::vector<G4double> energyDepVec;
    std::vector<G4double> positionZVec;
    std::vector<G4double> stepLengthVec;

    G4int stepNumber = 0;

    G4double stepLengthCumulative = 0.;

};




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


