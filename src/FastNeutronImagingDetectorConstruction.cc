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
/// \file FastNeutronImagingDetectorConstruction.cc
/// \brief Implementation of the FastNeutronImagingDetectorConstruction class

#include "FastNeutronImagingDetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Orb.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"

#include "G4RunManager.hh"

#include "G4VisAttributes.hh"

#include "G4GenericMessenger.hh"
#include "G4PhysicalConstants.hh"

#include "FastNeutronImagingLayerSD.hh"

#include "FastNeutronImagingDetectorMessenger.hh"

#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FastNeutronImagingDetectorConstruction::FastNeutronImagingDetectorConstruction()
: G4VUserDetectorConstruction(),  
  fSolidLayer(nullptr),
  fSolidCell(nullptr),
  fShieldPhy(nullptr),
  fSamplePhys(nullptr),
  fLayerLog(nullptr),
  fLayerUserLog(nullptr),
  fLayerThickness(20.),
  fSampleThickness(20.),
  fMessenger(nullptr)
{
   // DefineCommands();
    fMessenger = new FastNeutronImagingDetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FastNeutronImagingDetectorConstruction::~FastNeutronImagingDetectorConstruction()
{
    delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* FastNeutronImagingDetectorConstruction::Construct()
{  
   ConstructMaterial();
   Air_mat = G4Material::GetMaterial("G4_AIR");
   //G4Material* Water_mat = G4Material::GetMaterial("G4_WATER");
   G4Material* layer_mat = G4Material::GetMaterial("Mylar");
   //G4Material* cell_mat = G4Material::GetMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
   G4Material* cell_mat = G4Material::GetMaterial("BC408");
   PE_mat = G4Material::GetMaterial("PE_bulk");
   lead_mat = G4Material::GetMaterial("Lead");
   Al_mat = G4Material::GetMaterial("aluminium");


  // Envelope parameters
  //

   
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //     
  // World
  //
  G4double world_sizeXY = 1.5*m;
  G4double world_sizeZ  = 7.5*m;
  
  G4Box* SolidWorld =
    new G4Box("World",                       //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size
      
  G4LogicalVolume* LogicWorld =
    new G4LogicalVolume(SolidWorld,          //its solid
                        Air_mat,           //its material
                        "World");            //its name
                                   
  G4VPhysicalVolume* PhysWorld =
    new G4PVPlacement(nullptr,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      LogicWorld,            //its logical volume
                      "World",               //its name
                      nullptr,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking

 
  //     
  // Layer
  //

  /*******************************************/
  G4int cell_num_x = 16, cell_num_y = 16;
  G4double cell_dx = 2*mm, cell_dy = 2*mm, cell_dz = fLayerThickness;
  G4double gap_dx = 0.02*mm, gap_dy = 0.02*mm;
  G4double layer_dx = cell_num_x*cell_dx+(cell_num_x+1)*gap_dx,
           layer_dy = cell_num_y*cell_dy+(cell_num_y+1)*gap_dy,
           layer_dz = fLayerThickness+2*gap_dx;

  G4ThreeVector pos_layer = G4ThreeVector(0, 0, 0);


  fSolidLayer =
    new G4Box("Layer",
              0.5*layer_dx, 0.5*layer_dy,
              0.5*layer_dz);
                      
  G4LogicalVolume* LogicLayer =
    new G4LogicalVolume(fSolidLayer,         //its solid
                        layer_mat,          //its material
                        "Layer");           //its name
               
  new G4PVPlacement(nullptr,                       //no rotation
                    pos_layer,                    //at position
                    LogicLayer,             //its logical volume
                    "Layer",                //its name
                    LogicWorld,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
  //
  // cells
  //

  fSolidCell = new G4Box("Cell", 0.5*cell_dx, 0.5*cell_dy, 0.5*cell_dz);

  G4LogicalVolume* LogicalCell = new G4LogicalVolume(fSolidCell, cell_mat, "cell");

  //
  // the arrangement of the cells needs adjustment if using the method below,
  // the adjustment causes by the origin is the center of mother volumn,
  //
  G4int NOCopy = 0;
  G4double adjust_displacement_x = (cell_num_x*cell_dx+(cell_num_x+1)*gap_dx)*0.5,
           adjust_displacement_y = (cell_num_y*cell_dy+(cell_num_y+1)*gap_dy)*0.5;


  for (G4int i=0; i!=cell_num_x; ++i) {

      for (G4int j=0; j!=cell_num_y; ++j) {


          new G4PVPlacement(nullptr,
                            G4ThreeVector((i+0.5)*cell_dx+(i+1)*gap_dx
                                         -adjust_displacement_x,
                                          ((j+0.5)*cell_dy+(j+1)*gap_dy
                                         -adjust_displacement_y),
                                          gap_dx),
                            LogicalCell,
                            "cell",
                            LogicLayer,
                            false,
                            NOCopy,
                            checkOverlaps);
            ++NOCopy;
      }

  }

  /********************/
  // optical boundry

  // layer and cell surface
  G4OpticalSurface* opScintMylarSurface = new G4OpticalSurface("ScintMylarSurface");
  opScintMylarSurface->SetType(dielectric_dielectric);
  opScintMylarSurface->SetFinish(polishedfrontpainted);
  opScintMylarSurface->SetModel(unified);

  G4LogicalSkinSurface* ScintMylarSurface =
          new G4LogicalSkinSurface("ScintMylarSurface", LogicalCell, opScintMylarSurface);

  G4OpticalSurface* opticalSurface = dynamic_cast <G4OpticalSurface*>
        (ScintMylarSurface->GetSurface(LogicalCell)->GetSurfaceProperty());
  if (opticalSurface) opticalSurface->DumpInfo();

  const G4int num = 2;
  G4double ephoton[num] = {2.034*eV, 4.136*eV};

  G4double specularlobe[num] = {0.9, 0.9};
  G4double specularspike[num] = {0.9, 0.9};
  G4double backscatter[num] = {0.9, 0.9};
  G4double reflectivity[num] = {1., 1.};
  G4double efficiency[num]   = {0.8, 1.0};

  G4MaterialPropertiesTable *myST1 = new G4MaterialPropertiesTable();

  //myST1->AddProperty("SPECULARSPIKECONSTANT", ephoton, specularspike, num);
  //myST1->AddProperty("SPECULARLOBECONSTANT", ephoton, specularlobe, num);
  //myST1->AddProperty("BACKSCATTERCONSTANT", ephoton, backscatter, num);
  myST1->AddProperty("REFLECTIVITY", ephoton, reflectivity, num);
  myST1->AddProperty("EFFICIENCY",   ephoton, efficiency,   num);

  opScintMylarSurface->SetMaterialPropertiesTable(myST1);

  /************************/
  // layer and air surface
  G4OpticalSurface* opLayerAirSurface = new G4OpticalSurface("LayerAirSurface");
  opLayerAirSurface->SetType(dielectric_dielectric);
  opLayerAirSurface->SetFinish(polished);
  opLayerAirSurface->SetModel(unified);

  G4LogicalSkinSurface* LayerAirSurface =
          new G4LogicalSkinSurface("LayerAirSurface", LogicLayer, opScintMylarSurface);

  opticalSurface = dynamic_cast <G4OpticalSurface*>
        (LayerAirSurface->GetSurface(LogicLayer)->GetSurfaceProperty());
  if (opticalSurface) opticalSurface->DumpInfo();

  G4double reflectivityAirLayer[num] = {0.1, 0.1};
  G4double efficiencyAirLayer[num]   = {0.8, 1.0};

  G4MaterialPropertiesTable *myST2 = new G4MaterialPropertiesTable();

  myST2->AddProperty("REFLECTIVITY", ephoton, reflectivityAirLayer, num);
  myST2->AddProperty("EFFICIENCY",   ephoton, efficiencyAirLayer,   num);

  opLayerAirSurface->SetMaterialPropertiesTable(myST2);


      /****************/

  /******************************
  G4double cell_radius = 0.5*mm, cell_thickness = 1*mm;
  G4double layer_dx = 50*mm,
           layer_dy = 50*mm,
           layer_dz = fLayerThickness;
  G4ThreeVector pos_layer = G4ThreeVector(0, 0, 0);


  fSolidLayer =
    new G4Box("Layer",
              0.5*layer_dx, 0.5*layer_dy,
              0.5*layer_dz);

  G4LogicalVolume* LogicLayer =
    new G4LogicalVolume(fSolidLayer,         //its solid
                        layer_mat,          //its material
                        "Layer");           //its name

  new G4PVPlacement(nullptr,                       //no rotation
                    pos_layer,                    //at position
                    LogicLayer,             //its logical volume
                    "Layer",                //its name
                    LogicWorld,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

  fSolidCell = new G4Tubs("cell",0,cell_radius,cell_thickness,0.,2*pi*rad);

  G4LogicalVolume* LogicalCell = new G4LogicalVolume(fSolidCell,cell_mat,"cell");

  G4int NOCopy = 0;
  G4double adjust_displacement_x = 49*(cell_radius+0.1*mm),
           adjust_displacement_y = 49*std::sqrt(3)*(cell_radius+0.1*mm);


  for (G4int i=0; i!=45; ++i) {

      for (G4int j=0; j!=45; ++j) {

          switch (j%2) {
            case 0: new G4PVPlacement(nullptr,
                        G4ThreeVector((i-1)*2*(cell_radius+0.1*mm)
                                       -adjust_displacement_x,
                                       ((j-1)*std::sqrt(3)*(cell_radius+0.1*mm)
                                       -adjust_displacement_y),
                                        0),
                                     LogicalCell,
                                     "cell",
                                     LogicLayer,
                                     false,
                                     NOCopy,
                                     checkOverlaps);
                     ++NOCopy;
              break;

          case 1: new G4PVPlacement(nullptr,
                                    G4ThreeVector((i-1)*2*(cell_radius+0.1*mm)
                                                   -adjust_displacement_x-(cell_radius+0.1*mm),
                                                   ((j-1)*std::sqrt(3)*(cell_radius+0.1*mm)
                                                   -adjust_displacement_y),
                                                    0),
                                                 LogicalCell,
                                                 "cell",
                                                 LogicLayer,
                                                 false,
                                                 NOCopy,
                                                 checkOverlaps);
                                 ++NOCopy;
                          break;


          }

      }

  }
***************/


  //     
  // Sample
  //

  G4double sample_dx = 2*cm, sample_dy = 2*cm, sample_dz = fSampleThickness;
  G4double sample_r_outer = 1.5*cm, sample_r_inner = 0.75*cm, sample_thickness = fSampleThickness;
  G4ThreeVector default_sample_pos = G4ThreeVector(0, 0, -25*cm);

  // box sample
  Sample_box_origin =
    new G4Box("Sample_box",                      //its name
              0.5*sample_dx, 0.5*sample_dy, 0.5*sample_dz); //its size

  // tub sample
  Sample_tube_origin =
    new G4Tubs("Sample_tube",sample_r_inner,sample_r_outer,sample_thickness,
               0.,2*pi*rad);

  //Sample Lshape
  G4double L1_dx =1*cm, L1_dy = 2*cm, L1_dz = fSampleThickness;
  G4double L2_dx =3*cm, L2_dy = 1*cm, L2_dz = fSampleThickness;
  Sample_Lshape_part1 =
          new G4Box("L_part1",
                    0.5*L1_dx,0.5*L1_dy,0.5*L1_dz);
  Sample_Lshape_part2 =
          new G4Box("L_part2",
                    0.5*L2_dx,0.5*L2_dy,0.5*L2_dz);

  Sample_Lshape_origin =
          new G4UnionSolid("Sample_Lshape",Sample_Lshape_part1,Sample_Lshape_part2,
                           nullptr,G4ThreeVector(-1.0*cm,-1.5*cm,0));

  // default is tube, PE material
  fSampleShape = Sample_tube_origin;
  fSampleMaterial = PE_mat;

  fSampleLog =
    new G4LogicalVolume(fSampleShape,         //its solid
                        fSampleMaterial,          //its material
                        "Sample");           //its name

  G4RotationMatrix* sample_rotation = new G4RotationMatrix();
  sample_rotation->rotateX(0.*deg);
               
  fSamplePhys = new G4PVPlacement(sample_rotation,                       //no rotation
                    default_sample_pos,                    //at position
                    fSampleLog,             //its logical volume
                    "Sample",                //its name
                    LogicWorld,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
/*
  //
  // this is the outter stainless tube
  // make sure the sample is set to Sample_tube and the material is PE_bulk
  //


  G4Tubs* metalTube = new G4Tubs("Sample_tube",
                                 sample_r_outer+1*mm,
                                 sample_r_outer+3.5*mm,
                                 1.5*cm,
                                 0.,2*pi*rad);
  G4LogicalVolume* metalTubeLog = new G4LogicalVolume(metalTube,
                                                      Al_mat,
                                                      "metalTube");

  new G4PVPlacement(sample_rotation,
                    default_sample_pos,
                    metalTubeLog,
                    "metalTubePhy",
                    LogicWorld,
                    false,
                    checkOverlaps);
*/

  /***************************
  //
  // the PE shield
  //
  G4Material* shield_mat = PE_mat;
  G4double outterSize_x = 25*cm, outterSize_y = 25*cm, outterSize_z = 10*cm;
  G4double innerSize_x = layer_dx+5*mm, innerSize_y = layer_dy+5*mm,
           innerSize_z = outterSize_z;


  G4Box* outter = new G4Box("outter",
                            0.5*outterSize_x,
                            0.5*outterSize_y,
                            0.5*outterSize_z);
  G4Box* inner = new G4Box("inner",
                           0.5*innerSize_x,
                           0.5*innerSize_y,
                           0.5*innerSize_z);

  G4SubtractionSolid* solidShield = new G4SubtractionSolid("outter-inner",
                                                           outter,inner);
  G4LogicalVolume* logicShield = new G4LogicalVolume(solidShield,shield_mat,"Shield");

  fShieldPhy = new G4PVPlacement(nullptr,
                    pos_layer,
                    logicShield,
                    "Shield",
                    LogicWorld,
                    false,
                    0,
                    checkOverlaps);
  ***************************/

  // Set Layer as scoring volume
  //
  fLayerUserLog = LogicalCell;

  LogicLayer->SetVisAttributes(new G4VisAttributes(G4Colour(1.,1.,1.)));
  LogicalCell->SetVisAttributes(new G4VisAttributes(G4Colour(1.,1.,0)));
  fSampleLog->SetVisAttributes(new G4VisAttributes(G4Colour(0,1.,0)));


  //
  //always return the physical World
  //
  return PhysWorld;
}

void FastNeutronImagingDetectorConstruction::ConstructMaterial(){
    // Get NistManager
    G4NistManager* man = G4NistManager::Instance();

    G4Element* H = man->FindOrBuildElement("H",false);
    G4Element* C = man->FindOrBuildElement("C",false);
    G4Element* Ba = man->FindOrBuildElement("Ba",false);
    G4Element* S  = man->FindOrBuildElement("S",false);
    G4Element* O = man->FindOrBuildElement("O",false);
    G4Element* Al = man->FindOrBuildElement("Al",false);
    G4Element* Pb = man->FindOrBuildElement("Pb",false);

    man->FindOrBuildMaterial("G4_AIR");

    man->FindOrBuildMaterial("G4_WATER");

    G4Material* EJ200 = man->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
    EJ200->GetIonisation()->SetBirksConstant(0.16*mm/MeV);
    man->FindOrBuildMaterial("G4_POLYETHYLENE");

    G4Material* BC408 = new G4Material("BC408", 1.032*g/cm3, 2);
    BC408->AddElement(C,27);
    BC408->AddElement(H,30);

    // define BaSO4

    G4Material* Barium_Sulfate = new G4Material("BaSO4",4.49*g/cm3,3);
    Barium_Sulfate->AddElement(Ba,1);
    Barium_Sulfate->AddElement(S,1);
    Barium_Sulfate->AddElement(O,4);

    // define Polyenthylene bulk material

    G4Material* PE_bulk = new G4Material("PE_bulk",0.96*g/cm3,2);
    PE_bulk->AddElement(C,2);
    PE_bulk->AddElement(H,4);

    // define Mylar membrane
    G4Material* Mylar = new G4Material("Mylar",1.39*g/cm3,3);
    Mylar->AddElement(C,10);
    Mylar->AddElement(H,8);
    Mylar->AddElement(O,4);

    // define lead
    G4Material* Lead = new G4Material("Lead",11.34*g/cm3,1);
    Lead->AddElement(Pb,1);

    // define aluminium
    G4Material* aluminium = new G4Material("aluminium",2.70*g/cm3,1);
    aluminium->AddElement(Al,1);

    const G4int nEntries = 12;
    G4double PhotonEnergy[nEntries] =
    { 2.08*eV, 2.38*eV, 2.58*eV, 2.7*eV,
      2.76*eV, 2.82*eV, 2.92*eV, 2.95*eV,
      3.02*eV, 3.1*eV, 3.26*eV, 3.44*eV};
    G4double Scintillation1[nEntries] =
    { 0.00, 0.03, 0.17, 0.40,
      0.55, 0.83, 1.00, 0.84,
      0.49, 0.20, 0.07, 0.04 };
    G4double RefractiveIndex1[nEntries] =
    { 1.58, 1.58, 1.58, 1.58,
      1.58, 1.58, 1.58, 1.58,
      1.58, 1.58, 1.58, 1.58 };
    G4double Absorption1[nEntries] =
    { 380*cm, 380*cm, 380*cm, 380*cm,
      380*cm, 380*cm, 380*cm, 380*cm,
      380*cm, 380*cm, 380*cm, 380*cm };

    G4MaterialPropertiesTable* myMPT1 = new G4MaterialPropertiesTable();
    myMPT1->AddProperty("RINDEX", PhotonEnergy, RefractiveIndex1, nEntries);
    myMPT1->AddProperty("ABSLENGTH", PhotonEnergy, Absorption1, nEntries);
    myMPT1->AddProperty("FASTCOMPONENT", PhotonEnergy, Scintillation1, nEntries);
    myMPT1->AddConstProperty("SCINTILLATIONYIELD",100./MeV);
    myMPT1->AddConstProperty("RESOLUTIONSCALE",1.0);
    myMPT1->AddConstProperty("FASTTIMECONSTANT", 0.9*ns);
    myMPT1->AddConstProperty("SLOWTIMECONSTANT",2.1*ns);
    myMPT1->AddConstProperty("YIELDRATIO",1.);

    BC408->SetMaterialPropertiesTable(myMPT1);

    BC408->GetIonisation()->SetBirksConstant(0.154*mm/MeV);

    G4double photonEnergy[] =
              { 2.034*eV, 2.068*eV, 2.103*eV, 2.139*eV,
                2.177*eV, 2.216*eV, 2.256*eV, 2.298*eV,
                2.341*eV, 2.386*eV, 2.433*eV, 2.481*eV,
                2.532*eV, 2.585*eV, 2.640*eV, 2.697*eV,
                2.757*eV, 2.820*eV, 2.885*eV, 2.954*eV,
                3.026*eV, 3.102*eV, 3.181*eV, 3.265*eV,
                3.353*eV, 3.446*eV, 3.545*eV, 3.649*eV,
                3.760*eV, 3.877*eV, 4.002*eV, 4.136*eV };

    const G4int nEntriesAir = sizeof(photonEnergy)/sizeof(G4double);

    G4double refractiveIndex2[] =
              { 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
                1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
                1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
                1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
                1.00, 1.00, 1.00, 1.00 };

    G4MaterialPropertiesTable* myMPT2 = new G4MaterialPropertiesTable();
    myMPT2->AddProperty("RINDEX", photonEnergy, refractiveIndex2, nEntriesAir);


}

void FastNeutronImagingDetectorConstruction::ConstructSDandField(){

    auto sdManger = G4SDManager::GetSDMpointer();
    sdManger->SetVerboseLevel(1);
    G4String SDname;

    auto LayerUser = new FastNeutronImagingLayerSD(SDname = "/LayerUser");
    sdManger->AddNewDetector(LayerUser);
    fLayerUserLog->SetSensitiveDetector(LayerUser);

    auto detectorEnergy = new G4MultiFunctionalDetector("detectorEnergy");
    sdManger->AddNewDetector(detectorEnergy);
    G4VPrimitiveScorer* primitive = new G4PSEnergyDeposit("energyDeposition");
    detectorEnergy->RegisterPrimitive(primitive);
    SetSensitiveDetector(fLayerUserLog,detectorEnergy);

}

void FastNeutronImagingDetectorConstruction::SetLayerThickness(G4double val){
    fLayerThickness = val;
    //G4double shield_outter_z = 10*cm;
    fSolidLayer->SetZHalfLength(0.5*fLayerThickness);
    fSolidCell->SetZHalfLength(0.5*fLayerThickness);
    //fShieldPhy->
      //SetTranslation(G4ThreeVector(0,0,0.5*fLayerThickness-0.5*shield_outter_z));



    // tell G4RunManager that we change the geometry
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

void FastNeutronImagingDetectorConstruction::SetSampleThickness(G4double val){
    fSampleThickness = val;
    Sample_box_origin->SetZHalfLength(0.5*fSampleThickness);
    Sample_tube_origin->SetZHalfLength(0.5*fSampleThickness);
    Sample_Lshape_part1->SetZHalfLength(0.5*fSampleThickness);
    Sample_Lshape_part2->SetZHalfLength(0.5*fSampleThickness);

    // tell G4RunManager that we change the geometry
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

void FastNeutronImagingDetectorConstruction::SetSampleMaterial(G4int material){
    fSampleMaterialNo = material;
    switch (fSampleMaterialNo) {
    case 0:
        fSampleLog->SetMaterial(Air_mat);
        break;
    case 1:
        fSampleLog->SetMaterial(PE_mat);
        break;
    case 2:
        fSampleLog->SetMaterial(lead_mat);
    }
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

void FastNeutronImagingDetectorConstruction::SetSampleShape(G4int shape){
    fSampleShapeNo = shape;
    switch (fSampleShapeNo) {
    case 1:
        fSampleLog->SetSolid(Sample_box_origin);
        break;
    case 2:
        fSampleLog->SetSolid(Sample_tube_origin);
        break;
    case 3:
        fSampleLog->SetSolid(Sample_Lshape_origin);
        break; }
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

/*
void FastNeutronImagingDetectorConstruction::DefineCommands(){

    fMessenger = new G4GenericMessenger(this,
                                        "/FastNeutronImagingCamera/detector/",
                                        "Detector control");

    //  command layerthickness
    auto& LayerThicknessCmd
      = fMessenger->DeclareMethodWithUnit("LayerThickness","mm",
                                  &FastNeutronImagingDetectorConstruction::SetLayerThickness,
                                  "Set thickness of the Layer1.");
    LayerThicknessCmd.SetParameterName("Thickness", true);
    LayerThicknessCmd.SetRange("Thickness>=0. && Thickness<100.");
    LayerThicknessCmd.SetDefaultValue("20.");

    //  command samplethickness
    auto& SampleThicknessCmd
      = fMessenger->DeclareMethodWithUnit("SampleThickness","mm",
                                  &FastNeutronImagingDetectorConstruction::SetSampleThickness,
                                  "Set thickness of the Sample.");
    SampleThicknessCmd.SetParameterName("Thickness", true);
    SampleThicknessCmd.SetRange("Thickness>=0. && Thickness<100.");
    SampleThicknessCmd.SetDefaultValue("20.");


}*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
