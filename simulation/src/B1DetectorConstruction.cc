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
/// \file B1DetectorConstruction.cc
/// \brief Implementation of the B1DetectorConstruction class

#include "B1DetectorConstruction.hh"
#include "B1SteppingAction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::B1DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::~B1DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B1DetectorConstruction::Construct()
{  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
  
  // Envelope parameters
  //
  G4double env_sizeXY = 10*m, env_sizeZ = 10*m;
  G4Material* water = nist->FindOrBuildMaterial("G4_WATER");
   
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;
  G4Material* vacuum =
    new G4Material("Vacuum",      //Name as String
                   1,             //Atomic Number,  in this case we use 1 for hydrogen  
                   1.008*g/mole,  //Mass per Mole "Atomic Weight"  1.008*g/mole for Hydoren 
                   1.e-25*g/cm3,  //Density of Vaccuum  *Cant be Zero, Must be small insted 
                   kStateGas,     //kStateGas for Gas
                   2.73*kelvin,   //Temperature for gas
                   1.e-25*g/cm3); //Pressure for Vaccum

  //     
  // World
  //
  G4double world_sizeXY = 1.2*env_sizeXY;
  G4double world_sizeZ  = 1.2*env_sizeZ;
  G4Material* air = nist->FindOrBuildMaterial("G4_AIR");
//  G4NistManager* man = G4NistManager::Instance();
//  man->FindOrBuildMaterial("G4_CARBON_DIOXIDE");

  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size
      
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        air,           //its material
                        "World");            //its name


                                   
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking
                     
  //     
  // Water bath
  //  
   G4RotationMatrix *flipmatrix = new G4RotationMatrix();
   flipmatrix->rotateX(90*deg);

  G4Tubs* solidBath =    
    new G4Tubs("WaterBath",                    //its name
        0*m,2*m,0.5*1*m,0*deg,360*deg); //its size
      
  G4LogicalVolume* logicBath =                         
    new G4LogicalVolume(solidBath,            //its solid
                        water,             //its material
                        "WaterBath");         //its name
               
  new G4PVPlacement(flipmatrix,                       //rotation
                    G4ThreeVector(),         //at (0,0,0)
                    logicBath,                //its logical volume
                    "WaterBath",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking


  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
