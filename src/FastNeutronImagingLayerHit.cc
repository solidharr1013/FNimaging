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
/// \file FastNeutronImagingLayerHit.cc
/// \brief Implementation of the FastNeutronImagingLayerHit class

#include "FastNeutronImagingLayerHit.hh"

#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4AttDefStore.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4UIcommand.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal G4Allocator<FastNeutronImagingLayerHit>* FastNeutronImagingLayerHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FastNeutronImagingLayerHit::FastNeutronImagingLayerHit()
: G4VHit(), 
   fEdep(0.), fLocalPos(0), fWorldPos(0), fProcessName(""), fCellID(-1)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//FastNeutronImagingLayerHit::FastNeutronImagingLayerHit(G4int layerID)
//: G4VHit(),
 // fLayerID(layerID), fTime(0.), fLocalPos(0), fWorldPos(0), fProcessName("")
//{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FastNeutronImagingLayerHit::~FastNeutronImagingLayerHit()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FastNeutronImagingLayerHit::FastNeutronImagingLayerHit(const FastNeutronImagingLayerHit &right)
: G4VHit(),
  //fLayerID(right.fLayerID),
  fEdep(right.fEdep),
  fLocalPos(right.fLocalPos),
  fWorldPos(right.fWorldPos),
  fProcessName(right.fProcessName),
  fCellID(right.fCellID)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const FastNeutronImagingLayerHit& FastNeutronImagingLayerHit::operator=(const FastNeutronImagingLayerHit &right)
{
  //fLayerID = right.fLayerID;
  fEdep = right.fEdep;
  fLocalPos = right.fLocalPos;
  fWorldPos = right.fWorldPos;
  fProcessName = right.fProcessName;
  fCellID = right.fCellID;
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool FastNeutronImagingLayerHit::operator==(const FastNeutronImagingLayerHit &/*right*/) const
{
  return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FastNeutronImagingLayerHit::Draw()
{
  auto visManager = G4VVisManager::GetConcreteInstance();
  if (! visManager) return;

  G4Circle circle(fWorldPos);
  circle.SetScreenSize(2);
  circle.SetFillStyle(G4Circle::filled);
  G4Colour colour(0.4,0.4,0.);
  G4VisAttributes attribs(colour);
  circle.SetVisAttributes(attribs);
  visManager->Draw(circle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const std::map<G4String,G4AttDef>* FastNeutronImagingLayerHit::GetAttDefs() const
{
  G4bool isNew;
  auto store = G4AttDefStore::GetInstance("FastNeutronImagingLayerHit",isNew);

  if (isNew) {
      (*store)["HitType"] 
        = G4AttDef("HitType","Hit Type","Physics","","G4String");
      
      (*store)["ID"]
        = G4AttDef("ID","ID","Physics","","G4int");
      
      (*store)["Energy"]
        = G4AttDef("Energy","Energy","Physics","G4BestUnit","G4double");
      
      (*store)["Pos"] 
        = G4AttDef("Pos", "Position", "Physics","G4BestUnit","G4ThreeVector");

      (*store)["Process"]
              =G4AttDef("Process", "ProcessName", "Physics" ,"", "G4String");
  }
  
  return store;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::vector<G4AttValue>* FastNeutronImagingLayerHit::CreateAttValues() const
{
  auto values = new std::vector<G4AttValue>;
  
  values
    ->push_back(G4AttValue("HitType","DriftChamberHit",""));
  values
    ->push_back(G4AttValue("ID",G4UIcommand::ConvertToString(fCellID),""));
  values
    ->push_back(G4AttValue("EnergyDeposition",G4BestUnit(fEdep,"Energy"),""));
  values
    ->push_back(G4AttValue("Pos",G4BestUnit(fWorldPos,"Length"),""));
  values
    ->push_back(G4AttValue("Process", fProcessName,""));
  
  return values;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FastNeutronImagingLayerHit::Print()
{
  G4cout << "] : time " << fEdep/keV
  << " (nsec) --- local (x,y) " << fLocalPos.x()
  << ", " << fLocalPos.y() << G4endl
  << "process name " << fProcessName << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
