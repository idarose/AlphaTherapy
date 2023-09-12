#include "RunAction.hh"
#include "EventAction.hh"

#include "G4AnalysisManager.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(EventAction* eventAction)
: fEventAction(eventAction)
{
    // set printing event number per each event
    G4RunManager::GetRunManager()->SetPrintProgress(1);

    auto analysisManager = G4AnalysisManager::Instance();

    analysisManager->SetVerboseLevel(1);
    analysisManager->SetNtupleMerging(true);

    //Creating ntuple to hold interaction information
    analysisManager->CreateNtuple("B4", "Energydep Info");
    analysisManager->CreateNtupleDColumn(0, "EnergyDeps", fEventAction->GetEnergyDepVec());
    analysisManager->CreateNtupleIColumn(0, "CellIDs", fEventAction->GetCellIDVec());
    analysisManager->CreateNtupleIColumn(0, "VolumeTypes", fEventAction->GetVolumeTypeVec());
    analysisManager->CreateNtupleDColumn(0, "KineticEnergy", fEventAction->GetKineticEnergyVec());
    analysisManager->CreateNtupleIColumn(0, "ParticleType", fEventAction->GetParticleTypeVec());
    analysisManager->CreateNtupleDColumn(0, "InteractionTime", fEventAction->GetInteractionTime());
    analysisManager->FinishNtuple();
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
  G4String fileName = "Output_212Pb_C4_2_Membrane.root";
  analysisManager->OpenFile(fileName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* /*run*/)
{
  auto analysisManager = G4AnalysisManager::Instance();

  //Writing to and closing ROOT file
  analysisManager->Write();
  analysisManager->CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


