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
/// \file FastNeutronImagingRunAction.hh
/// \brief Definition of the FastNeutronImagingRunAction class

#ifndef FastNeutronImagingRunAction_h
#define FastNeutronImagingRunAction_h 1

#include "G4UserRunAction.hh"
#include "G4Accumulable.hh"
#include "globals.hh"

//txt format output
#include "iostream"

class G4Run;

/// Run action class
///
/// In EndOfRunAction(), it calculates the dose in the selected volume 
/// from the energy deposit accumulated via stepping and event actions.
/// The computed dose is then printed on the screen.

class FastNeutronImagingRunAction : public G4UserRunAction
{
  public:
    FastNeutronImagingRunAction();
    virtual ~FastNeutronImagingRunAction();

    // virtual G4Run* GenerateRun();
    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);
    void AddEdep (G4double edep);
    void AddEdepGamma(G4double edepGamma);
    void AddEdepProton(G4double edepProton);
    void AddEdepElectron(G4double edepElectron);
    void AddEdepBe9(G4double edepBe9);
    void AddEdepC12(G4double edepC12);
    void AddEdepAlpha(G4double edepAlpha);

    G4double GetEdep(){ return fEdep.GetValue();}
    G4double GetEdepGamma(){return fEdepGamma.GetValue();}

    void AddFastNeutronImagingNumber (G4int FastNeutronImagingnumber);
    void AddTotalNumber (G4int totalnumber);
    void AddGammaNumber(G4int GammaNumber);


    //txt format output
    void OutputFile(std::string outputFileName) {fOutput = outputFileName;}

  private:
    G4Accumulable<G4double> fEdep;
    G4Accumulable<G4double> fEdep2;
    G4Accumulable<G4double> fEdepGamma;
    G4Accumulable<G4double> fEdepProton;
    G4Accumulable<G4double> fEdepElectron;
    G4Accumulable<G4double> fEdepBe9;
    G4Accumulable<G4double> fEdepC12;
    G4Accumulable<G4double> fEdepAlpha;

    G4Accumulable<G4int> fFastNeutronImagingNumber;
    G4Accumulable<G4int> fTotalNumber;
    G4Accumulable<G4int> fGammaNumber;

    //txt format output
    std::string fOutput;
};

#endif

