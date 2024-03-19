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

#include "RunMessenger.hh"

#include "RunAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4Material.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunMessenger::RunMessenger(RunAction* run)
 : G4UImessenger(),
   fRun(run),
   fDirectory(0),
   fOutNameCmd(0)
{
  fDirectory = new G4UIdirectory("/Sim/");
  fDirectory->SetGuidance("UI commands of Sim");

  fOutNameCmd = new G4UIcmdWithAString("/Sim/setOutputFileName",this);
  fOutNameCmd->SetGuidance("Select output root (directory and) file name.");
  fOutNameCmd->SetParameterName("choice",false);
  fOutNameCmd->AvailableForStates(G4State_Idle);

  fSetSeed = new G4UIcmdWithAnInteger("/Sim/setSeed",this);
  fSetSeed->SetGuidance("Set seed of run");
  fSetSeed->SetParameterName("choice",false);
  fSetSeed->AvailableForStates(G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunMessenger::~RunMessenger()
{
  delete fDirectory;
  delete fOutNameCmd;
  delete fSetSeed;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if( command == fOutNameCmd ) {
    fRun->SetOutputFileName(newValue);
  }
  else if( command == fSetSeed ) {
    G4cout << "HERE IN RUN MESSENGER, SetSeed: " << newValue << G4endl;
    fRun->SetSeed(fSetSeed->GetNewIntValue(newValue));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4String RunMessenger::GetCurrentValue(G4UIcommand * command)
{
  G4String ans;
  if( command == fOutNameCmd ){
    // ans=fRun->GetOutName();
  }
  return ans;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
