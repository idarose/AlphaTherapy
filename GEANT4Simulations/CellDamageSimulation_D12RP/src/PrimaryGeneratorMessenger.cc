//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
//

#include "PrimaryGeneratorMessenger.hh"

#include "PrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4Material.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(PrimaryGeneratorAction* pga)
 : G4UImessenger(),
   fPrimaryGeneratorAction(pga),
   fDirectory(0),
   fInitialRadionuclide_Z(0),
   fInitialRadionuclide_A(0),
   fInitialRadionuclide_excitationEnergy(0),
   fInitialRadionuclide_define(0),
   fSampleActivity(0),
   fCellLineName(0),
   fDecayCurve(0)
{
  fDirectory = new G4UIdirectory("/Sim/");
  fDirectory->SetGuidance("UI commands of Sim");

  fInitialRadionuclide_Z = new G4UIcmdWithAnInteger("/Sim/SetInitialRadionuclide_Z",this);
  fInitialRadionuclide_Z->SetGuidance("Set Z of initial radionuclide to be simulated.");
  fInitialRadionuclide_Z->SetParameterName("choice",false);
  fInitialRadionuclide_Z->AvailableForStates(G4State_Idle);

  fInitialRadionuclide_A = new G4UIcmdWithAnInteger("/Sim/SetInitialRadionuclide_A",this);
  fInitialRadionuclide_A->SetGuidance("Set A of initial radionuclide to be simulated.");
  fInitialRadionuclide_A->SetParameterName("choice",false);
  fInitialRadionuclide_A->AvailableForStates(G4State_Idle);

  fInitialRadionuclide_excitationEnergy = new G4UIcmdWithADouble("/Sim/SetInitialRadionuclide_excitationEnergy",this);
  fInitialRadionuclide_excitationEnergy->SetGuidance("Set excitation energy of initial radionuclide to be simulated.");
  fInitialRadionuclide_excitationEnergy->SetParameterName("choice",false);
  fInitialRadionuclide_excitationEnergy->AvailableForStates(G4State_Idle);

  fInitialRadionuclide_location = new G4UIcmdWithAnInteger("/Sim/SetInitialRadionuclide_location",this);
  fInitialRadionuclide_location->SetGuidance("Set location of initial radionuclide to be simulated.");
  fInitialRadionuclide_location->SetParameterName("choice",false);
  fInitialRadionuclide_location->AvailableForStates(G4State_Idle);

  fInitialRadionuclide_define = new G4UIcmdWithoutParameter("/Sim/DefineInitialRadionuclide",this);
  fInitialRadionuclide_define->SetGuidance("Define the initial radionuclide to be simulated (command must be executed after Z, A and excitationEnergy have been set).");
  fInitialRadionuclide_define->AvailableForStates(G4State_Idle);

  fSampleActivity = new G4UIcmdWithAnInteger("/Sim/SetSampleActivity", this);
  fSampleActivity->SetGuidance("Set the sample activity in kBq/mL.");
  fSampleActivity->AvailableForStates(G4State_Idle);

  fCellLineName = new G4UIcmdWithAString("/Sim/SetCellLineName",this);
  fCellLineName->SetGuidance("Set the name of the cell line.");
  fCellLineName->AvailableForStates(G4State_Idle);

  fDecayCurve = new G4UIcmdWithoutParameter("/Sim/SetDecayCurveRadionuclide",this);
  fDecayCurve->SetGuidance("Load number of decays as a function of time in the given analysisregion for the radionuclide chosen. (command must be executed after sample activity and cell line name have been set).");
  fDecayCurve->AvailableForStates(G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
{
  delete fDirectory;
  delete fInitialRadionuclide_Z;
  delete fInitialRadionuclide_A;
  delete fInitialRadionuclide_excitationEnergy;
  delete fInitialRadionuclide_define;
  delete fSampleActivity;
  delete fCellLineName;
  delete fDecayCurve;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{

  if( command == fInitialRadionuclide_Z ) {
    fPrimaryGeneratorAction->SetInitialRadionuclide_Z(fInitialRadionuclide_Z->GetNewIntValue(newValue));
  }
  if( command == fInitialRadionuclide_A ) {
    fPrimaryGeneratorAction->SetInitialRadionuclide_A(fInitialRadionuclide_A->GetNewIntValue(newValue));
  }
  if( command == fInitialRadionuclide_excitationEnergy ) {
    fPrimaryGeneratorAction->SetInitialRadionuclide_excitationEnergy(fInitialRadionuclide_excitationEnergy->GetNewDoubleValue(newValue));
  }
  if( command == fInitialRadionuclide_location ) {
    fPrimaryGeneratorAction->SetInitialRadionuclide_location(fInitialRadionuclide_location->GetNewIntValue(newValue));
  }
  if( command == fInitialRadionuclide_define ) {
    fPrimaryGeneratorAction->DefineInitialRadionuclide();
  }
  if( command == fSampleActivity ) {
    fPrimaryGeneratorAction->SetSampleActivity(fSampleActivity->GetNewIntValue(newValue));
  }
  if(command == fCellLineName ) {
    fPrimaryGeneratorAction->SetCellLineName(newValue);
  }
  if(command == fDecayCurve ) {
    fPrimaryGeneratorAction->SetDecayCurveRadionuclide();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

