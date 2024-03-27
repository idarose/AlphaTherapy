#include "RunAction.hh"
#include "EventAction.hh"

#include "G4AnalysisManager.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(EventAction* eventAction)
: fEventAction(eventAction),
fOutputFileName("SimulationOutput.root"),
fSeedN(1)
{
    // set printing event number per each event
    G4RunManager::GetRunManager()->SetPrintProgress(1);

    auto analysisManager = G4AnalysisManager::Instance();

    analysisManager->SetVerboseLevel(1);
    analysisManager->SetNtupleMerging(true);

    //Creating ntuple to hold interaction information
    analysisManager->CreateNtuple("B4", "Energydep Info");
    analysisManager->CreateNtupleDColumn(0, "EnergyDep", fEventAction->GetEnergyDepVec());
    analysisManager->CreateNtupleDColumn(0, "PositionZ", fEventAction->GetPositionZVec());
    analysisManager->CreateNtupleDColumn(0, "StepLength", fEventAction->GetStepLengthVec());
    analysisManager->FinishNtuple();

    fRunMessenger = new RunMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* /*run*/)
{
    auto analysisManager = G4AnalysisManager::Instance();

    // Open an output file
    analysisManager->OpenFile(fOutputFileName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* /*run*/)
{
  auto analysisManager = G4AnalysisManager::Instance();

  //Writing to and closing ROOT file
  analysisManager->Write();
  analysisManager->CloseFile();
}

//----------------------------------------------------------------------------

void RunAction::SetOutputFileName(G4String name)
{
  fOutputFileName = name;
}

//----------------------------------------------------------------------------

void RunAction::SetSeed(G4int seedN)
{
    fSeedN = seedN;
    G4cout << "fSeedN (OCLRunAction): " << fSeedN << G4endl;
    G4Random::setTheSeed(fSeedN);
}
