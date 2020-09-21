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
/// \file FastNeutronImagingEventAction.cc
/// \brief Implementation of the FastNeutronImagingEventAction class

#include "FastNeutronImagingEventAction.hh"
#include "FastNeutronImagingRunAction.hh"
#include "FastNeutronImagingSteppingAction.hh"
#include "FastNeutronImagingAnalysis.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4Track.hh"
#include "G4RunManager.hh"
#include "G4SteppingManager.hh"

#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4THitsMap.hh"

#include "FastNeutronImagingLayerHit.hh"
#include "FastNeutronImagingLayerSD.hh"

#include "G4SystemOfUnits.hh"
#include "random"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace {

// Utility function which finds a hit collection with the given Id
// and print warnings if not found
G4VHitsCollection* GetHC(const G4Event* event, G4int collId) {
  auto hce = event->GetHCofThisEvent();
  if (!hce) {
      G4ExceptionDescription msg;
      msg << "No hits collection of this event found." << G4endl;
      G4Exception("FastNeutronImagingEventAction::EndOfEventAction()",
                  "FastNeutronImagingCode001", JustWarning, msg);
      return nullptr;
  }

  auto hc = hce->GetHC(collId);
  if ( ! hc) {
    G4ExceptionDescription msg;
    msg << "Hits collection " << collId << " of this event not found." << G4endl;
    G4Exception("FastNeutronImagingEventAction::EndOfEventAction()",
                "FastNeutronImagingCode001", JustWarning, msg);
  }
  return hc;
}
}

FastNeutronImagingEventAction::FastNeutronImagingEventAction(FastNeutronImagingRunAction* runAction)
: G4UserEventAction(),
  fRunAction(runAction),
  fEdep(0.),
  fEdepGamma(0.),
  fEdepProton(0.),
  fEdepElectron(0.),
  fEdepBe9(0.),
  fEdepC12(0.),
  fEdepAlpha(0.),
  fFastNeutronImagingNumber(0),
  fTotalNumber(0),
  fGammaNumber(0),
  fLayerID(-1),
  fLayerUserID(-1),
  fCollID_detectorEnergy(-1)
{
     G4RunManager::GetRunManager()->SetPrintProgress(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FastNeutronImagingEventAction::~FastNeutronImagingEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FastNeutronImagingEventAction::BeginOfEventAction(const G4Event*)
{    
  fEdep = 0.;
  fEdepGamma = 0.;
  fEdepProton = 0.;
  fEdepElectron = 0.;
  fEdepBe9 = 0.;
  fEdepC12 = 0.;
  fEdepAlpha = 0.;

  fFastNeutronImagingNumber = 0;
  fTotalNumber = 0;
  fGammaNumber = 0;

  if(fLayerUserID == -1){
      auto sdManager = G4SDManager::GetSDMpointer();

      fLayerUserID = sdManager->GetCollectionID("LayerUser/LayerUserColl");

  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FastNeutronImagingEventAction::EndOfEventAction(const G4Event* event)
{   
    auto analysisManager = G4AnalysisManager::Instance();

    auto HCE = event->GetHCofThisEvent();
    if(!HCE) return;

    if(fCollID_detectorEnergy<0){
        auto SDMan = G4SDManager::GetSDMpointer();
        fCollID_detectorEnergy = SDMan->GetCollectionID("detectorEnergy/energyDeposition");
    }

    G4double eThreshold = 300*eV;

    G4THitsMap<G4double>* evtMap =
                       static_cast<G4THitsMap<G4double>*>(HCE->GetHC(fCollID_detectorEnergy));

    std::map<G4int,G4double*>::iterator itr;
    for(itr = evtMap->GetMap()->begin(); itr != evtMap->GetMap()->end(); ++itr){

        G4double CopyNo = static_cast<G4double>(itr->first);
        G4double eDep = *(itr->second);

        if(eDep>eThreshold){
            analysisManager->FillH1(1,CopyNo);
           // G4cout << "CopyNo is " << CopyNo << G4endl;
        }
    }

    auto hcLayerUser = GetHC(event,fLayerUserID);
    if(!hcLayerUser) return;

    if(hcLayerUser->GetSize() != 0){
        ++fTotalNumber;
        fRunAction->AddTotalNumber(fTotalNumber);
    }


        //G4int columnID = CopyNo/16;
        //G4int rowID = CopyNo%16;

        // the threshold used in experiment is set for total signal,
        // so set ethreshold for total energy deposited.


  // accumulate statistics in run action
  fRunAction->AddEdep(fEdep);
  fRunAction->AddEdepGamma(fEdepGamma);
  fRunAction->AddEdepProton(fEdepProton);
  fRunAction->AddEdepElectron(fEdepElectron);
  fRunAction->AddEdepBe9(fEdepBe9);
  fRunAction->AddEdepC12(fEdepC12);
  fRunAction->AddEdepAlpha(fEdepAlpha);

  if(fEdep > 0.2*MeV){
      ++fFastNeutronImagingNumber;
      fRunAction->AddFastNeutronImagingNumber(fFastNeutronImagingNumber);
     // G4cout << "the eDep of this event is " << fEdep << G4endl;

      analysisManager->FillH1(0,addFluction(fEdep));
  }

  if(fEdepGamma>0.){
      analysisManager->FillH1(2,fEdepGamma);
  }
  if(fEdepProton>0.){
      analysisManager->FillH1(3,addFluction(fEdepProton));
  }
  if(fEdepElectron>0.){
      analysisManager->FillH1(4,addFluction(fEdepElectron));
  }
  if(fEdepBe9>0.){
      analysisManager->FillH1(5,addFluction(fEdepBe9+fEdepAlpha));
  }
  if(fEdepC12>0.){
      analysisManager->FillH1(6,addFluction(fEdepC12));
  }
  if(fEdepAlpha>0.){
      analysisManager->FillH1(7,addFluction(fEdepAlpha));
  }

  fRunAction->AddGammaNumber(fGammaNumber);


}

G4double FastNeutronImagingEventAction::addFluction(G4double val){
    //calculate the coefficient of the resolution, based on 662keV
    G4double Resolution662 =0.4,
             coefficient = Resolution662*std::sqrt(661.6*keV);

    G4double valVar = coefficient*std::sqrt(val)/(2*std::sqrt(2*std::log(2)));

    std::random_device rd;
    std::mt19937 generator(rd());

    std::normal_distribution<G4double> distribution(val,valVar);
    G4double var_new = distribution(generator);

    return var_new;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
