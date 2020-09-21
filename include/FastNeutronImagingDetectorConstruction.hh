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
/// \file FastNeutronImagingDetectorConstruction.hh
/// \brief Definition of the FastNeutronImagingDetectorConstruction class

#ifndef FastNeutronImagingDetectorConstruction_h
#define FastNeutronImagingDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4Box;
class G4Tubs;
class G4UnionSolid;
class G4VSolid;
class G4Material;
class G4VPhysicalVolume;
class G4LogicalVolume;
class G4VSensitiveDetector;

class G4GenericMessenger;
class FastNeutronImagingDetectorMessenger;
/// Detector construction class to define materials and geometry.

class FastNeutronImagingDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    FastNeutronImagingDetectorConstruction();
    virtual ~FastNeutronImagingDetectorConstruction();

    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();

    void SetLayerThickness(G4double val);
    G4double GetLayerThickness() const { return fLayerThickness; }
    
    G4LogicalVolume* GetLayerLogicalVolume() const { return fLayerUserLog; }
    G4LogicalVolume* GetSampleLogicalVolume() const { return fSampleLog; }

    void SetSampleThickness(G4double val);
    G4double GetSampleThickness() const { return fSampleThickness; }

    void SetSampleMaterial(G4int material);
    G4int GetSampleMaterial() const {return fSampleMaterialNo;}

    void SetSampleShape(G4int shape);
    G4int GetSampleShape() const {return fSampleShapeNo;}

    void ConstructMaterial();

  private:
    //void DefineCommands();
    //G4GenericMessenger* fMessenger;
    G4Box* fSolidLayer;
    //G4Box* fSolidCell;
    G4Box* fSolidCell;
    G4VPhysicalVolume* fShieldPhy;
    G4VPhysicalVolume* fSamplePhys;

    G4Box* Sample_box_origin;
    G4Box* Sample_Lshape_part1;
    G4Box* Sample_Lshape_part2;
    G4Tubs* Sample_tube_origin;
    G4UnionSolid* Sample_Lshape_origin;

    G4VSolid* fSampleShape;

    G4Material* fSampleMaterial;

    G4Material* Air_mat;
    G4Material* PE_mat;
    G4Material* lead_mat;
    G4Material* Al_mat;

    G4int fSampleMaterialNo;

    G4int fSampleShapeNo;

  protected:
    G4LogicalVolume*  fLayerLog;
    G4LogicalVolume* fLayerUserLog;
    G4LogicalVolume* fSampleLog;

    G4double fLayerThickness;
    G4double fSampleThickness;

    FastNeutronImagingDetectorMessenger* fMessenger;


};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

