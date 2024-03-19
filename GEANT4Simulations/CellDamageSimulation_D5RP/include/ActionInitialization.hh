#ifndef B4aActionInitialization_h
#define B4aActionInitialization_h 1

#include "G4VUserActionInitialization.hh"

class DetectorConstruction;

/// Action initialization class.

class ActionInitialization : public G4VUserActionInitialization
{
  public:
    ActionInitialization(DetectorConstruction*);
    ~ActionInitialization() override;

    void BuildForMaster() const override;
    void Build() const override;

  private:
    DetectorConstruction* fDetConstruction = nullptr;
};


#endif


