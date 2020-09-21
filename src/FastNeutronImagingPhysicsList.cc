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
/// \file B3PhysicsList.cc
/// \brief Implementation of the B3PhysicsList class

#include "FastNeutronImagingPhysicsList.hh"

#include "G4DecayPhysics.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"
//#include "G4RadioactiveDecayPhysics.hh"
#include "G4OpticalPhysics.hh"
//#include "QGSP_BERT_HP.hh"
#include "QGSP_BIC_HP.hh"
#include "G4HadronElasticPhysicsHP.hh"
//#include "G4HadronPhysicsQGSP_BERT_HP.hh"
#include "G4HadronPhysicsQGSP_BIC_HP.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/**********************************************/

FastNeutronImagingPhysicsList::FastNeutronImagingPhysicsList()
: G4VModularPhysicsList(){
  SetVerboseLevel(1);

  // Default physics
  RegisterPhysics(new G4DecayPhysics());
  
  // Hadronic physics
  //RegisterPhysics(new G4HadronPhysicsQGSP_BERT_HP());
  RegisterPhysics(new G4HadronPhysicsQGSP_BIC_HP());
  RegisterPhysics(new G4HadronElasticPhysicsHP());

  // EM physics
  RegisterPhysics(new G4EmLivermorePhysics());

  // Optical physics
  G4OpticalPhysics* optical = new G4OpticalPhysics();
  //optical->Configure(kCerenkov,true);
  optical->Configure(kScintillation,true);
  optical->SetScintillationYieldFactor(1.);
  optical->SetTrackSecondariesFirst(kScintillation,true);
  optical->SetMaxNumPhotonsPerStep(1000);
  optical->Configure(kAbsorption,true);
  RegisterPhysics(optical);

  // Radioactive decay
  //RegisterPhysics(new G4RadioactiveDecayPhysics());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FastNeutronImagingPhysicsList::~FastNeutronImagingPhysicsList()
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FastNeutronImagingPhysicsList::SetCuts()
{
  G4VUserPhysicsList::SetCuts();
}  
/**************************************
#include "globals.hh"
#include "FastNeutronImagingPhysicsList.hh"
#include "FastNeutronImagingPhysicsListMessenger.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

#include "G4ProcessManager.hh"

#include "G4Cerenkov.hh"
#include "G4Scintillation.hh"
#include "G4OpAbsorption.hh"
#include "G4OpRayleigh.hh"
#include "G4OpMieHG.hh"
#include "G4OpBoundaryProcess.hh"

#include "G4LossTableManager.hh"
#include "G4EmSaturation.hh"

G4ThreadLocal G4int FastNeutronImagingPhysicsList::fVerboseLevel = 1;
G4ThreadLocal G4int FastNeutronImagingPhysicsList::fMaxNumPhotonStep = 20;
G4ThreadLocal G4Cerenkov* FastNeutronImagingPhysicsList::fCerenkovProcess = 0;
G4ThreadLocal G4Scintillation* FastNeutronImagingPhysicsList::fScintillationProcess = 0;
G4ThreadLocal G4OpAbsorption* FastNeutronImagingPhysicsList::fAbsorptionProcess = 0;
G4ThreadLocal G4OpRayleigh* FastNeutronImagingPhysicsList::fRayleighScatteringProcess = 0;
G4ThreadLocal G4OpMieHG* FastNeutronImagingPhysicsList::fMieHGScatteringProcess = 0;
G4ThreadLocal G4OpBoundaryProcess* FastNeutronImagingPhysicsList::fBoundaryProcess = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FastNeutronImagingPhysicsList::FastNeutronImagingPhysicsList()
 : G4VUserPhysicsList()
{
  fMessenger = new FastNeutronImagingPhysicsListMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FastNeutronImagingPhysicsList::~FastNeutronImagingPhysicsList() { delete fMessenger; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FastNeutronImagingPhysicsList::ConstructParticle()
{
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program.

  G4BosonConstructor bConstructor;
  bConstructor.ConstructParticle();

  G4LeptonConstructor lConstructor;
  lConstructor.ConstructParticle();

  G4MesonConstructor mConstructor;
  mConstructor.ConstructParticle();

  G4BaryonConstructor rConstructor;
  rConstructor.ConstructParticle();

  G4IonConstructor iConstructor;
  iConstructor.ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FastNeutronImagingPhysicsList::ConstructProcess()
{
  AddTransportation();
  ConstructDecay();
  ConstructEM();
  ConstructOp();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4Decay.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FastNeutronImagingPhysicsList::ConstructDecay()
{
  // Add Decay Process
  G4Decay* theDecayProcess = new G4Decay();
  auto particleIterator=GetParticleIterator();
  particleIterator->reset();
  while( (*particleIterator)() ){
    G4ParticleDefinition* particle = particleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    if (theDecayProcess->IsApplicable(*particle)) {
      pmanager ->AddProcess(theDecayProcess);
      // set ordering for PostStepDoIt and AtRestDoIt
      pmanager ->SetProcessOrdering(theDecayProcess, idxPostStep);
      pmanager ->SetProcessOrdering(theDecayProcess, idxAtRest);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4eMultipleScattering.hh"
#include "G4MuMultipleScattering.hh"
#include "G4hMultipleScattering.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4hIonisation.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FastNeutronImagingPhysicsList::ConstructEM()
{
  auto particleIterator=GetParticleIterator();
  particleIterator->reset();
  while( (*particleIterator)() ){
    G4ParticleDefinition* particle = particleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    if (particleName == "gamma") {
    // gamma
      // Construct processes for gamma
      pmanager->AddDiscreteProcess(new G4GammaConversion());
      pmanager->AddDiscreteProcess(new G4ComptonScattering());
      pmanager->AddDiscreteProcess(new G4PhotoElectricEffect());

    } else if (particleName == "e-") {
    //electron
      // Construct processes for electron
      pmanager->AddProcess(new G4eMultipleScattering(),-1, 1, 1);
      pmanager->AddProcess(new G4eIonisation(),       -1, 2, 2);
      pmanager->AddProcess(new G4eBremsstrahlung(),   -1, 3, 3);

    } else if (particleName == "e+") {
    //positron
      // Construct processes for positron
      pmanager->AddProcess(new G4eMultipleScattering(),-1, 1, 1);
      pmanager->AddProcess(new G4eIonisation(),       -1, 2, 2);
      pmanager->AddProcess(new G4eBremsstrahlung(),   -1, 3, 3);
      pmanager->AddProcess(new G4eplusAnnihilation(),  0,-1, 4);

    } else if( particleName == "mu+" ||
               particleName == "mu-"    ) {
    //muon
     // Construct processes for muon
     pmanager->AddProcess(new G4MuMultipleScattering(),-1, 1, 1);
     pmanager->AddProcess(new G4MuIonisation(),      -1, 2, 2);
     pmanager->AddProcess(new G4MuBremsstrahlung(),  -1, 3, 3);
     pmanager->AddProcess(new G4MuPairProduction(),  -1, 4, 4);

    } else {
      if ((particle->GetPDGCharge() != 0.0) &&
          (particle->GetParticleName() != "chargedgeantino") &&
          !particle->IsShortLived()) {
       // all others charged particles except geantino
       pmanager->AddProcess(new G4hMultipleScattering(),-1,1,1);
       pmanager->AddProcess(new G4hIonisation(),       -1,2,2);
     }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#include "G4Threading.hh"

void FastNeutronImagingPhysicsList::ConstructOp()
{
  fCerenkovProcess = new G4Cerenkov("Cerenkov");
  fCerenkovProcess->SetMaxNumPhotonsPerStep(fMaxNumPhotonStep);
  fCerenkovProcess->SetMaxBetaChangePerStep(10.0);
  fCerenkovProcess->SetTrackSecondariesFirst(true);
  fScintillationProcess = new G4Scintillation("Scintillation");
  fScintillationProcess->SetScintillationYieldFactor(1.);
  fScintillationProcess->SetTrackSecondariesFirst(true);
  fAbsorptionProcess = new G4OpAbsorption();
  fRayleighScatteringProcess = new G4OpRayleigh();
  fMieHGScatteringProcess = new G4OpMieHG();
  fBoundaryProcess = new G4OpBoundaryProcess();

  fCerenkovProcess->SetVerboseLevel(fVerboseLevel);
  fScintillationProcess->SetVerboseLevel(fVerboseLevel);
  fAbsorptionProcess->SetVerboseLevel(fVerboseLevel);
  fRayleighScatteringProcess->SetVerboseLevel(fVerboseLevel);
  fMieHGScatteringProcess->SetVerboseLevel(fVerboseLevel);
  fBoundaryProcess->SetVerboseLevel(fVerboseLevel);

  // Use Birks Correction in the Scintillation process
  if(G4Threading::IsMasterThread())
  {
    G4EmSaturation* emSaturation =
              G4LossTableManager::Instance()->EmSaturation();
      fScintillationProcess->AddSaturation(emSaturation);
  }

  auto particleIterator=GetParticleIterator();
  particleIterator->reset();
  while( (*particleIterator)() ){
    G4ParticleDefinition* particle = particleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
    if (fCerenkovProcess->IsApplicable(*particle)) {
      pmanager->AddProcess(fCerenkovProcess);
      pmanager->SetProcessOrdering(fCerenkovProcess,idxPostStep);
    }
    if (fScintillationProcess->IsApplicable(*particle)) {
      pmanager->AddProcess(fScintillationProcess);
      pmanager->SetProcessOrderingToLast(fScintillationProcess, idxAtRest);
      pmanager->SetProcessOrderingToLast(fScintillationProcess, idxPostStep);
    }
    if (particleName == "opticalphoton") {
      G4cout << " AddDiscreteProcess to OpticalPhoton " << G4endl;
      pmanager->AddDiscreteProcess(fAbsorptionProcess);
      pmanager->AddDiscreteProcess(fRayleighScatteringProcess);
      pmanager->AddDiscreteProcess(fMieHGScatteringProcess);
      pmanager->AddDiscreteProcess(fBoundaryProcess);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FastNeutronImagingPhysicsList::SetVerbose(G4int verbose)
{
  fVerboseLevel = verbose;

  fCerenkovProcess->SetVerboseLevel(fVerboseLevel);
  fScintillationProcess->SetVerboseLevel(fVerboseLevel);
  fAbsorptionProcess->SetVerboseLevel(fVerboseLevel);
  fRayleighScatteringProcess->SetVerboseLevel(fVerboseLevel);
  fMieHGScatteringProcess->SetVerboseLevel(fVerboseLevel);
  fBoundaryProcess->SetVerboseLevel(fVerboseLevel);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FastNeutronImagingPhysicsList::SetNbOfPhotonsCerenkov(G4int MaxNumber)
{
  fMaxNumPhotonStep = MaxNumber;

  fCerenkovProcess->SetMaxNumPhotonsPerStep(fMaxNumPhotonStep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FastNeutronImagingPhysicsList::SetCuts()
{
  //  " G4VUserPhysicsList::SetCutsWithDefault" method sets
  //   the default cut value for all particle types
  //
  SetCutsWithDefault();

  if (verboseLevel>0) DumpCutValuesTable();
}
******************************/
