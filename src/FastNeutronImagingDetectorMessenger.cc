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
/// \file FastNeutronImaging/src/FastNeutronImagingDetectorMessenger.cc
/// \brief Implementation of the FastNeutronImagingDetectorMessenger class
//
//

#include "FastNeutronImagingDetectorMessenger.hh"

#include "FastNeutronImagingDetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4Material.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FastNeutronImagingDetectorMessenger::FastNeutronImagingDetectorMessenger(FastNeutronImagingDetectorConstruction* det)
 : G4UImessenger(),
   fDetector(det),
   fDirectory(nullptr),
   fSampleMaterialCmd(nullptr),
   fSampleShapeCmd(nullptr),
   fLayerThicknessCmd(nullptr),
   fSampleThicknessCmd(nullptr)
{ 
  fDirectory = new G4UIdirectory("/FastNeutronImaging/");
  fDirectory->SetGuidance("UI commands of this example");
/*
  G4String matList;
  const G4MaterialTable* matTbl = G4Material::GetMaterialTable();
  for(size_t i=0;i<G4Material::GetNumberOfMaterials();i++)
  {
    matList += (*matTbl)[i]->GetName();
    matList += " ";
  }
*/
  fSampleMaterialCmd = new G4UIcmdWithAnInteger("/FastNeutronImaging/Sample/SetSampleMaterial",this);
  fSampleMaterialCmd->SetGuidance("Select material number.");
  fSampleMaterialCmd->SetParameterName("material",false);
  fSampleMaterialCmd->AvailableForStates(G4State_Idle);
  fSampleMaterialCmd->SetRange("material>=0");
  
  fSampleShapeCmd = new G4UIcmdWithAnInteger("/FastNeutronImaging/Sample/SetSampleShape",this);
  fSampleShapeCmd->SetGuidance("Select shape number.");
  fSampleShapeCmd->SetParameterName("shape",false);
  fSampleShapeCmd->AvailableForStates(G4State_Idle);
  fSampleShapeCmd->SetRange("shape>=0");

  fSampleThicknessCmd = new G4UIcmdWithADoubleAndUnit("/FastNeutronImaging/Sample/SampleThickness",this);
  fSampleThicknessCmd->SetGuidance("set Sample thickness.");
  fSampleThicknessCmd->SetParameterName("SampleThickness",true,true);
  fSampleThicknessCmd->SetDefaultUnit("mm");
  fSampleThicknessCmd->SetUnitCandidates("mm cm m");
  fSampleThicknessCmd->SetDefaultValue(20.);

  fLayerThicknessCmd = new G4UIcmdWithADoubleAndUnit("/FastNeutronImaging/detector/LayerThickness",this);
  fLayerThicknessCmd->SetGuidance("set Layer thickness.");
  fLayerThicknessCmd->SetParameterName("LayerThickness",true,true);
  fLayerThicknessCmd->SetDefaultUnit("mm");
  fLayerThicknessCmd->SetUnitCandidates("mm cm m");
  fLayerThicknessCmd->SetDefaultValue(20.);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FastNeutronImagingDetectorMessenger::~FastNeutronImagingDetectorMessenger()
{
  delete fSampleMaterialCmd;
  delete fSampleShapeCmd;
  delete fSampleThicknessCmd;
  delete fLayerThicknessCmd;
  delete fDirectory;  
}

void FastNeutronImagingDetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == fSampleMaterialCmd ) {
    fDetector->SetSampleMaterial(fSampleMaterialCmd->GetNewIntValue(newValue));
  
  } else if( command == fSampleShapeCmd ) {
    fDetector->SetSampleShape(fSampleShapeCmd->GetNewIntValue(newValue));

  } else if(command == fSampleThicknessCmd){
      fDetector->SetSampleThickness(fSampleThicknessCmd->GetNewDoubleValue(newValue));
  }
  else if (command == fLayerThicknessCmd) {
      fDetector->SetLayerThickness(fLayerThicknessCmd->GetNewDoubleValue(newValue));
  }

}

G4String FastNeutronImagingDetectorMessenger::GetCurrentValue(G4UIcommand * command)
{
  G4String ans;
  if( command == fSampleMaterialCmd ){
    ans=fSampleMaterialCmd->ConvertToString(fDetector->GetSampleMaterial());


  } else if( command == fSampleShapeCmd ) {
    ans=fSampleShapeCmd->ConvertToString(fDetector->GetSampleShape());
  }

  else if(command == fSampleThicknessCmd){
      ans = fSampleThicknessCmd->ConvertToStringWithBestUnit(fDetector->GetSampleThickness());
  }
  else if (command == fLayerThicknessCmd) {
      ans = fLayerThicknessCmd->ConvertToStringWithBestUnit(fDetector->GetLayerThickness());
  }

  return ans;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/*
void    FastNeutronImagingDetectorMessenger::UpdateMaterialList()
{
  G4String matList;
  const G4MaterialTable* matTbl = G4Material::GetMaterialTable();
  for(size_t i=0;i<G4Material::GetNumberOfMaterials();i++)
  {
    matList += (*matTbl)[i]->GetName();
    matList += " ";
  }

  if(fSampleMaterialCmd != nullptr) {
    fSampleMaterialCmd->SetCandidates(matList);
  }
}
*/
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
