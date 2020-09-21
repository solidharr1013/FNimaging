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
/// \file FastNeutronImagingLayerHit.hh
/// \brief Definition of the FastNeutronImagingLayerHit class

#ifndef FastNeutronImagingLayerHit_h
#define FastNeutronImagingLayerHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

class G4AttDef;
class G4AttValue;

/// Drift chamber hit
///
/// It records:
/// - the layer ID
/// - the particle time
/// - the particle local and global positions

class FastNeutronImagingLayerHit : public G4VHit
{
  public:
    FastNeutronImagingLayerHit();
    //FastNeutronImagingLayerHit(G4int layerID);
    FastNeutronImagingLayerHit(const FastNeutronImagingLayerHit &right);
    virtual ~FastNeutronImagingLayerHit();

    const FastNeutronImagingLayerHit& operator=(const FastNeutronImagingLayerHit &right);
    G4bool operator==(const FastNeutronImagingLayerHit &right) const;
    
    inline void *operator new(size_t);
    inline void operator delete(void *aHit);
    
    virtual void Draw();
    virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
    virtual std::vector<G4AttValue>* CreateAttValues() const;
    virtual void Print();

    //void SetLayerID(G4int z) { fLayerID = z; }
    //G4int GetLayerID() const { return fLayerID; }

    void SetEnergyDepositon(G4double edep) { fEdep = edep; }
    G4double GetEnergyDepositon() const { return fEdep; }

    void SetLocalPos(G4ThreeVector xyz) { fLocalPos = xyz; }
    G4ThreeVector GetLocalPos() const { return fLocalPos; }

    void SetWorldPos(G4ThreeVector xyz) { fWorldPos = xyz; }
    G4ThreeVector GetWorldPos() const { return fWorldPos; }

    void SetProcessName(G4String procName) { fProcessName = procName; }
    G4String GetProcessName() const { return fProcessName; }

    void SetCellID(G4int z) {fCellID = z;}
    G4int GetCellID() const {return fCellID;}

    
  private:
    //G4int fLayerID;
    G4int fCellID;
    G4double fEdep;
    G4ThreeVector fLocalPos;
    G4ThreeVector fWorldPos;
    G4String fProcessName;
};

using FastNeutronImagingLayerHitsCollection = G4THitsCollection<FastNeutronImagingLayerHit>;

extern G4ThreadLocal G4Allocator<FastNeutronImagingLayerHit>* FastNeutronImagingLayerHitAllocator;

inline void* FastNeutronImagingLayerHit::operator new(size_t)
{
  if (!FastNeutronImagingLayerHitAllocator) {
       FastNeutronImagingLayerHitAllocator = new G4Allocator<FastNeutronImagingLayerHit>;
  }
  return (void*)FastNeutronImagingLayerHitAllocator->MallocSingle();
}

inline void FastNeutronImagingLayerHit::operator delete(void* aHit)
{
  FastNeutronImagingLayerHitAllocator->FreeSingle((FastNeutronImagingLayerHit*) aHit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
