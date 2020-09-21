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
/// \file FastNeutronImagingRunAction.cc
/// \brief Implementation of the FastNeutronImagingRunAction class

#include "FastNeutronImagingRunAction.hh"
#include "FastNeutronImagingPrimaryGeneratorAction.hh"
#include "FastNeutronImagingDetectorConstruction.hh"
#include "FastNeutronImagingEventAction.hh"
#include "FastNeutronImagingAnalysis.hh"
// #include "FastNeutronImagingRun.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FastNeutronImagingRunAction::FastNeutronImagingRunAction()
: G4UserRunAction(),
  fEdep(0.),
  fEdep2(0.),
  fEdepGamma(0.),
  fEdepProton(0.),
  fEdepElectron(0.),
  fEdepBe9(0.),
  fEdepC12(0.),
  fEdepAlpha(0.),
  fFastNeutronImagingNumber(0),
  fTotalNumber(0),
  fGammaNumber(0)
{ 
  // add new units for dose
  // 
  const G4double milligray = 1.e-3*gray;
  const G4double microgray = 1.e-6*gray;
  const G4double nanogray  = 1.e-9*gray;  
  const G4double picogray  = 1.e-12*gray;
   
  new G4UnitDefinition("milligray", "milliGy" , "Dose", milligray);
  new G4UnitDefinition("microgray", "microGy" , "Dose", microgray);
  new G4UnitDefinition("nanogray" , "nanoGy"  , "Dose", nanogray);
  new G4UnitDefinition("picogray" , "picoGy"  , "Dose", picogray); 

  // Register accumulable to the accumulable manager
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->RegisterAccumulable(fEdep);
  accumulableManager->RegisterAccumulable(fEdep2); 
  accumulableManager->RegisterAccumulable(fEdepGamma);
  accumulableManager->RegisterAccumulable(fEdepProton);
  accumulableManager->RegisterAccumulable(fEdepElectron);
  accumulableManager->RegisterAccumulable(fEdepBe9);
  accumulableManager->RegisterAccumulable(fEdepC12);
  accumulableManager->RegisterAccumulable(fEdepAlpha);

  accumulableManager->RegisterAccumulable(fFastNeutronImagingNumber);
  accumulableManager->RegisterAccumulable(fTotalNumber);
  accumulableManager->RegisterAccumulable(fGammaNumber);

  auto analysisManager = G4AnalysisManager::Instance();

  G4double ticMax = 16.*MeV , binNum = 2000;

  analysisManager->SetVerboseLevel(1);
  analysisManager->SetFileName("FNImagingHisto");

  analysisManager->CreateH1("eDep","Deposited Energy",binNum,0.,ticMax,"MeV"); //id=0
  analysisManager->CreateH1("countMap","countMap",
                            16*16,0,16*16);
  analysisManager->CreateH1("eDep of Gamma","Deposited Energy",binNum,0.,ticMax*MeV,"MeV");
  analysisManager->CreateH1("eDep of Proton","Deposited Energy",binNum,0.,ticMax*MeV,"MeV");
  analysisManager->CreateH1("eDep of Electron","Deposited Energy",binNum,0.,ticMax*MeV,"MeV");
  analysisManager->CreateH1("eDep of Be9","Deposited Energy",binNum,0.,ticMax*MeV,"MeV");
  analysisManager->CreateH1("eDep of C12","Deposited Energy",binNum,0.,ticMax*MeV,"MeV");
  analysisManager->CreateH1("eDep of Alpha","Deposited Energy",binNum,0.,ticMax*MeV,"MeV");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FastNeutronImagingRunAction::~FastNeutronImagingRunAction()
{delete G4AnalysisManager::Instance();}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FastNeutronImagingRunAction::BeginOfRunAction(const G4Run*)
{ 
  // inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);

  // reset accumulables to their initial values
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Reset();

  auto analysisManager = G4AnalysisManager::Instance();

  analysisManager->OpenFile();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FastNeutronImagingRunAction::EndOfRunAction(const G4Run* run)
{
  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;

  G4double edepTest = fEdep.GetValue();
  G4cout << "the test energy value is " << edepTest << G4endl;

  // Merge accumulables 
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Merge();

  // Compute dose = total energy deposit in a run and its variance
  //
  G4double edep  = fEdep.GetValue();
  G4double edep2 = fEdep2.GetValue();
  // cast int to double type
  G4double NOTotalFastNeutronImagingEvent = fFastNeutronImagingNumber.GetValue();
  G4double NOTotalEvent = fTotalNumber.GetValue();
  G4double NOGammaEvent = fGammaNumber.GetValue();

  G4double rms = edep2 - edep*edep/nofEvents;
  if (rms > 0.) rms = std::sqrt(rms); else rms = 0.;  

  const FastNeutronImagingDetectorConstruction* detectorConstruction
   = static_cast<const FastNeutronImagingDetectorConstruction*>
     (G4RunManager::GetRunManager()->GetUserDetectorConstruction());


  G4double FastNeutronImagingEventRatio = NOTotalFastNeutronImagingEvent/NOTotalEvent;

  // Run conditions
  //  note: There is no primary generator action object for "master"
  //        run manager for multi-threaded mode.
  const FastNeutronImagingPrimaryGeneratorAction* generatorAction
   = static_cast<const FastNeutronImagingPrimaryGeneratorAction*>
     (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  G4String runCondition;
  if (generatorAction)
  {
    const G4ParticleGun* particleGun = generatorAction->GetParticleGun();
    runCondition += particleGun->GetParticleDefinition()->GetParticleName();
    runCondition += " of ";
    G4double particleEnergy = particleGun->GetParticleEnergy();
    runCondition += G4BestUnit(particleEnergy,"Energy");
  }

  G4double thicknessLayer1 = detectorConstruction->GetLayerThickness();
  G4double thicknessSample = detectorConstruction->GetSampleThickness();
  G4String SampleName = detectorConstruction->GetSampleLogicalVolume()->GetSolid()->GetName();
  G4String SampleMaterial = detectorConstruction->GetSampleLogicalVolume()->GetMaterial()->GetName();

  // Print
  //  
  if (IsMaster()) {
    G4cout
     << G4endl
     << "--------------------End of Global Run-----------------------";
  }
  else {
    G4cout
     << G4endl
     << "--------------------End of Local Run------------------------";
  }
  
  G4cout
     << G4endl
     << " The thickness of layer1 is " << G4BestUnit(thicknessLayer1, "Length") << G4endl
     << " The Sample Shape is " << SampleName << G4endl
     << " The Sample Material is " << SampleMaterial << G4endl
     << " The thickness of sample is " << G4BestUnit(thicknessSample,"Length") << G4endl
     << " The run consists of " << nofEvents << " "<< runCondition
     << G4endl
     << " Cumulated energy per run, in scoring volume : "
     << G4BestUnit(edep,"Energy") << " energy despositon per event = " << G4BestUnit(edep/NOTotalFastNeutronImagingEvent,"Energy")
     << G4endl
     << "The number of total event is " << NOTotalEvent << G4endl
     << "The number of FastNeutronImaging event is " << NOTotalFastNeutronImagingEvent << G4endl
     << "The FastNeutronImaging event ratio is " << std::setw(3) << FastNeutronImagingEventRatio
     << G4endl
     << "the number of Secondary Gamma is " << NOGammaEvent << G4endl
     << "------------------------------------------------------------"
     << G4endl
     << G4endl;


  //txt format output
  std::ofstream SingleFastNeutronImagingEfficiency;


    SingleFastNeutronImagingEfficiency.open(fOutput + "SingleFastNeutronImagingEfficiency.txt", std::ios_base::app);
    if (SingleFastNeutronImagingEfficiency.is_open()){

      SingleFastNeutronImagingEfficiency
              << "The run consists of " << nofEvents << " "<< "; "
              << "The layer thickness is " << G4BestUnit(thicknessLayer1, "Length") << "; "
              << "The Sample Shape is " << SampleName << "; "
              << "The Sample Material is " << SampleMaterial << "; "
              << "The thickness of sample is " << G4BestUnit(thicknessSample,"Length") << "; "
              << "the total event number is " << NOTotalEvent << "; "
              << "the effective event number is " << NOTotalFastNeutronImagingEvent << "; "
              << "the efficient is " << std::setw(4) << FastNeutronImagingEventRatio << "; "
              << std::endl;
      }
    else SingleFastNeutronImagingEfficiency << "Unable to open file" << std::endl;

  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FastNeutronImagingRunAction::AddEdep(G4double edep)
{
  fEdep  += edep;
  fEdep2 += edep*edep;
}

void FastNeutronImagingRunAction::AddEdepGamma(G4double edepGamma){
    fEdepGamma += edepGamma;
}

void FastNeutronImagingRunAction::AddEdepProton(G4double edepProton){
    fEdepProton += edepProton;
}

void FastNeutronImagingRunAction::AddEdepElectron(G4double edepElectron){
    fEdepElectron += edepElectron;
}

void FastNeutronImagingRunAction::AddEdepBe9(G4double edepBe9){
    fEdepBe9 += edepBe9;
}

void FastNeutronImagingRunAction::AddEdepC12(G4double edepC12){
    fEdepC12 += edepC12;
}

void FastNeutronImagingRunAction::AddEdepAlpha(G4double edepAlpha){
    fEdepAlpha += edepAlpha;
}

void FastNeutronImagingRunAction::AddFastNeutronImagingNumber(G4int FastNeutronImagingnumber){
    fFastNeutronImagingNumber += FastNeutronImagingnumber;
}

void FastNeutronImagingRunAction::AddTotalNumber(G4int totalnumber){
    fTotalNumber += totalnumber;
}

void FastNeutronImagingRunAction::AddGammaNumber(G4int GammaNumber){
    fGammaNumber += GammaNumber;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

