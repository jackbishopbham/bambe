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
/// \file B1PrimaryGeneratorAction.hh
/// \brief Definition of the B1PrimaryGeneratorAction class

#ifndef B1PrimaryGeneratorAction_h
#define B1PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"
#include "TGraph.h"
#include "TGraph2D.h"

class G4ParticleGun;
class G4Event;
class G4Box;

/// The primary generator action class with particle gun.
///
/// The default kinematic is a 6 MeV gamma, randomly distribued 
/// in front of the phantom across 80% of the (X,Y) phantom size.

class B1PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    B1PrimaryGeneratorAction();    
    virtual ~B1PrimaryGeneratorAction();
    G4double GetNeutron();
    double InteractionE(double);
    double Eloss(double,double);
    double InteractionTheta(double);
    double ConvertCMtoLab(double,double);
    double NeutronEnergy(double,double,double);
    // method from the base class
    virtual void GeneratePrimaries(G4Event*);         
        TGraph *gXS_0;
        TGraph *gXS_1;
        TGraph *gXS_2;
        TGraph *gstopping;
        TGraph *gXS_t;
        TGraph2D *g_gs;
        TGraph2D *gp;
        TGraph2D *gpp;
        double d2r=atan(1.)/45.;
        double pi=4.*atan(1.);
        double mb=4.002603;//Mass 4He
        double mt=9.012182;//Mass 9Be
        double ml=1.0008665;//Mass neutron
        double mh=12.;//Mass 12C
        double amu=931.4941;
        double Q=+5.702;
        double max_depth = 0.025;//mm
        double Einitial=5.48556;//MeV
        double Eth=0. ;//Threshold energy
        double thetaCM,thetalab;
        bool secondks=false;
        bool excited=false;
        bool excited2=false;
  
    // method to access particle gun
    const G4ParticleGun* GetParticleGun() const { return fParticleGun; }
  
  private:
    G4ParticleGun*  fParticleGun; // pointer a to G4 gun class
    G4Box* fEnvelopeBox;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
