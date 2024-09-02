#include "PrimaryGeneratorAction.hh"
#include "PrimaryGeneratorMessenger.hh"
#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4IonTable.hh"


// https://stackoverflow.com/questions/6942273/how-to-get-a-random-element-from-a-c-container

#include  <random>
#include  <iterator>
#include  <TFile.h>
#include  <TROOT.h>

#include <TMath.h>
#include <Math/Interpolator.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>
#include <ctime>
#include <tuple>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* detConstruction):
fDetConstruction(detConstruction)
{
    TF1::DefaultAddToGlobalList(false);
    //-------------------------------
    G4int n_particle = 1;
    fParticleGun  = new G4ParticleGun(n_particle);

    fPrimaryGeneratorMessenger = new PrimaryGeneratorMessenger(this);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    auto particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("e-");

    fParticleGun->SetParticleDefinition(particleDefinition);


    fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,0.));

    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));

    fParticleGun->SetParticleEnergy(initialAlphaEnergy*MeV);

    fParticleGun->GeneratePrimaryVertex(anEvent);
}

//----------------------------------------------------------------------------

void PrimaryGeneratorAction::SetInitialAlphaEnergy(G4double value)
{
    initialAlphaEnergy = value;
}

