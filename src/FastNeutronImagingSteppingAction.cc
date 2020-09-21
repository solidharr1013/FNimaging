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
/// \file FastNeutronImagingSteppingAction.cc
/// \brief Implementation of the FastNeutronImagingSteppingAction class

#include "FastNeutronImagingSteppingAction.hh"
#include "FastNeutronImagingEventAction.hh"
#include "FastNeutronImagingDetectorConstruction.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FastNeutronImagingSteppingAction::FastNeutronImagingSteppingAction(FastNeutronImagingEventAction* eventAction)
: G4UserSteppingAction(),
  fEventAction(eventAction),
  fLayerUserLog(nullptr),
  fEventNumber(-1),
  fGammaNumber(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FastNeutronImagingSteppingAction::~FastNeutronImagingSteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FastNeutronImagingSteppingAction::UserSteppingAction(const G4Step* step)
{
  if ( !fLayerUserLog) {
   const FastNeutronImagingDetectorConstruction* detectorConstruction
      = static_cast<const FastNeutronImagingDetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    fLayerUserLog = detectorConstruction->GetLayerLogicalVolume();

  }


  // get volume of the current step
  G4LogicalVolume* volume 
    = step->GetPreStepPoint()->GetTouchableHandle()
      ->GetVolume()->GetLogicalVolume();
/*
  if(step->GetTrack()->GetParticleDefinition()->GetParticleName()
          == "gamma"){
      G4String GammaCreateProcess = step->GetTrack()->GetCreatorProcess()
              ->GetProcessName();
      G4double GammaKineticEnergy = step->GetTrack()->GetKineticEnergy();
      G4double GammaTotalEnergy = step->GetTrack()->GetTotalEnergy();
      G4cout<< "the gamma create by " << GammaCreateProcess << G4endl
            << "the Kinetic energy is " << GammaKineticEnergy << G4endl
            << "the Total energy is " << GammaTotalEnergy << G4endl;

  }
*/
  G4int eventNumber = G4RunManager::GetRunManager()->GetCurrentEvent()
          ->GetEventID();
  if(eventNumber!=fEventNumber){
      fEventNumber = eventNumber;
      fGammaNumber = 0;
  }

  const std::vector<const G4Track*>* secondary = step->GetSecondaryInCurrentStep();

  G4int GammaNumber = 0;
  if(secondary->size()>0){
      for(unsigned int i = 0; i<secondary->size(); ++i){
          if(secondary->at(i)->GetParentID() > 0){
              if(secondary->at(i)->GetDynamicParticle()->GetParticleDefinition()
                      ->GetParticleName() == "gamma"){
                  GammaNumber++;

              }
          }
      }
  }
  fEventAction->AddGammaNumber(GammaNumber);

  // check if we are in scoring volume
  if ( volume != fLayerUserLog ) return;

  // get information of the step
  G4StepPoint* preStep = step->GetPreStepPoint();
  G4StepPoint* postStep = step->GetPostStepPoint();
  G4Track* atrack = step->GetTrack();
  G4String particleName = atrack->GetDynamicParticle()->GetParticleDefinition()->GetParticleName();

  G4String PresentVolumeName = atrack->GetTouchable()->GetVolume()->GetName();
  G4String NextVolumeName = atrack->GetNextTouchable()->GetVolume()->GetName();

  G4String PreStepVolumeName = preStep->GetTouchable()->GetVolume()->GetName();
  G4String PostStepVolumeName = postStep->GetTouchable()->GetVolume()->GetName();

 /* if using the prestepPoint without this code, it will get Segmentation Fault(because it
  might point to a empty pointer).  */
  G4String ProcessName = "";
  if(postStep!= nullptr)
  {
    const G4VProcess* proc = postStep->GetProcessDefinedStep();
    if(proc != nullptr)
    {
      ProcessName = proc->GetProcessName();
    }
  }


  if(volume == fLayerUserLog){

 G4double edepStep = step->GetTotalEnergyDeposit();

 const std::vector<const G4Track*>* secondaryInLayer = step->GetSecondaryInCurrentStep();
 G4double edepStepGamma = 0.,
         edepStepProton = 0.,
         edepStepElectron = 0.,
         edepStepBe9 = 0.,
         edepStepC12 = 0.,
         edepStepAlpha = 0.;

 if(secondaryInLayer->size()>0){
     for(unsigned int i = 0; i<secondaryInLayer->size(); ++i){
         auto secondaryPresent = secondaryInLayer->at(i);
         if(secondaryPresent->GetParentID()>0){
             if(secondaryPresent->GetDynamicParticle()->
                     GetParticleDefinition()->GetParticleName() == "gamma"){
                 edepStepGamma = secondaryInLayer->at(i)->GetKineticEnergy();
                 fEventAction->AddEdepGamma(edepStepGamma);
             }
             if(secondaryPresent->GetDynamicParticle()->
                     GetParticleDefinition()->GetParticleName() == "proton"){
                 edepStepProton = secondaryInLayer->at(i)->GetKineticEnergy();
                 fEventAction->AddEdepProton(edepStepProton);
             }
             if(secondaryPresent->GetDynamicParticle()->
                     GetParticleDefinition()->GetParticleName() == "e-"){
                 edepStepElectron = secondaryInLayer->at(i)->GetKineticEnergy();
                 fEventAction->AddEdepElectron(edepStepElectron);
             }
             if(secondaryPresent->GetDynamicParticle()->
                     GetParticleDefinition()->GetParticleName() == "Be9"){
                 edepStepBe9 = secondaryInLayer->at(i)->GetKineticEnergy();
                 fEventAction->AddEdepBe9(edepStepBe9);
             }
             if(secondaryPresent->GetDynamicParticle()->
                     GetParticleDefinition()->GetParticleName() == "C12"){
                 edepStepC12 = secondaryInLayer->at(i)->GetKineticEnergy();
                 fEventAction->AddEdepC12(edepStepC12);
             }
             if(secondaryPresent->GetDynamicParticle()->
                     GetParticleDefinition()->GetParticleName() == "alpha"){
                 edepStepAlpha = secondaryInLayer->at(i)->GetKineticEnergy();
                 fEventAction->AddEdepAlpha(edepStepAlpha);
             }
         }
     }
 }

 fEventAction->AddEdep(edepStep);


}

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

