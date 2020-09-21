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
/// \file FastNeutronImagingEventAction.hh
/// \brief Definition of the FastNeutronImagingEventAction class

#ifndef FastNeutronImagingEventAction_h
#define FastNeutronImagingEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class FastNeutronImagingRunAction;


/// Event action class
///

class FastNeutronImagingEventAction : public G4UserEventAction
{
  public:
    FastNeutronImagingEventAction(FastNeutronImagingRunAction* runAction);
    virtual ~FastNeutronImagingEventAction();

    virtual void BeginOfEventAction(const G4Event* event);
    virtual void EndOfEventAction(const G4Event*);

    void AddEdep(G4double edep) { fEdep += edep; }
    void AddEdepGamma(G4double edepGamma){ fEdepGamma += edepGamma; }
    void AddEdepProton(G4double edepProton){fEdepProton += edepProton;}
    void AddEdepElectron(G4double edepElectron){fEdepElectron += edepElectron;}
    void AddEdepBe9(G4double edepBe9){fEdepBe9 += edepBe9;}
    void AddEdepC12(G4double edepC12){fEdepC12 += edepC12;}
    void AddEdepAlpha(G4double edepAlpha){fEdepAlpha += edepAlpha;}

    void AddFastNeutronImagingNumber(G4int FastNeutronImagingnumber) { fFastNeutronImagingNumber += FastNeutronImagingnumber; }
    void AddTotalNumber(G4int totalnumber) {fTotalNumber += totalnumber;}
    void AddGammaNumber(G4int GammaNumber) {fGammaNumber += GammaNumber;}

    G4double addFluction(G4double val);

    G4int GetTotalNumber() { return fTotalNumber; }

  private:
    FastNeutronImagingRunAction* fRunAction;
    G4double     fEdep;
    G4double fEdepGamma;
    G4double fEdepProton;
    G4double fEdepElectron;
    G4double fEdepBe9;
    G4double fEdepC12;
    G4double fEdepAlpha;

    G4int        fFastNeutronImagingNumber;
    G4int        fTotalNumber;
    G4int fGammaNumber;

    G4int fLayerID;
    G4int fLayerUserID;

    G4int fCollID_detectorEnergy;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
