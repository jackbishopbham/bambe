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
/// \file B1PrimaryGeneratorAction.cc
/// \brief Implementation of the B1PrimaryGeneratorAction class

#include "B1PrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
	double InteractionE(double);
	double InteractionTheta(double);
	void makesumgraph();
	double Eloss(double,double);
	double ConvertCMtoLab(double Eint,double thetaCM);
	double NeutronEnergy(double Eint, double thetalab, double thetaCM);

B1PrimaryGeneratorAction::B1PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(0), 
  fEnvelopeBox(0)
{
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);

  // default particle kinematic
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle
    = particleTable->FindParticle(particleName="neutron");
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleEnergy(9.*MeV);
  gXS_0 = new TGraph("XS_0.txt");
  gXS_1 = new TGraph("XS_1.txt");
  gXS_2 = new TGraph("XS_2.txt");
  gstopping = new TGraph("stopping_power.txt");
  gXS_t = new TGraph("XS_t.txt"); //Atomic Data and Nuclear Data Tables, Vol. 15, No. 1, January 1975
  g_gs = new TGraph2D("n0_ang.dat");
  gp  = new TGraph2D("n1_ang.dat");
  gpp = new TGraph2D("n2_ang.dat");
  d2r=atan(1.)/45.;
  pi=4.*atan(1.);
  mb=4.002603;//Mass 4He
  mt=9.012182;//Mass 9Be
  ml=1.0008665;//Mass neutron
  mh=12.;//Mass 12C
  amu=931.4941;
  Q=+5.702;
  max_depth = 0.025;//mm
  Einitial=5.48556;//MeV
  Eth=0. ;//Threshold energy
  secondks=false;
  excited=false;
  excited2=false;


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1PrimaryGeneratorAction::~B1PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of ecah event
  //

  // In order to avoid dependence of PrimaryGeneratorAction
  // on DetectorConstruction class we get Envelope volume
  // from G4LogicalVolumeStore.
  fParticleGun->SetParticlePosition(G4ThreeVector(0,0*cm,0*cm));
  fParticleGun->SetParticleEnergy(GetNeutron());

  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));

  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double B1PrimaryGeneratorAction::GetNeutron() {
	double Einitial = 5.48556;//Main peak
	if(G4UniformRand()<0.131) Einitial = 5.4428;//Satellite peak
	//Sample Ep_interaction
	double Eint = InteractionE(Einitial);
	//Sample theta CM
	thetaCM = InteractionTheta(Eint);
	if(thetaCM==0) return 0;//Skip this guy?
	thetalab = ConvertCMtoLab(Eint,thetaCM);
	double En = NeutronEnergy(Eint,thetalab,thetaCM);
	if(En==0) return 0;//Skip this guy
	return En;
}

double B1PrimaryGeneratorAction::InteractionE(double Ep) {
	//Rejection sampling for z_interaction based on XS(Ep) and stopping power
	//Then return Ep for z
	bool trials=true;
	double XSmax = 1000;
	int loopcounter=0;
	while(trials) {//Randomly choose a depth into the target (rather than energy to keep same target thickness)
		double trial_z = max_depth*G4UniformRand();
		//Calc XS at Z
		double energy = Eloss(Ep,trial_z);
		if(energy<=Eth) continue;
		double XS = gXS_t->Eval(energy);
//		cout<<gXS_1->Eval(energy)<<"\t"<<gXS_0->Eval(energy)<<endl;
		double choose_chan = G4UniformRand();
		double BRn0=gXS_0->Eval(energy)/(gXS_0->Eval(energy)+gXS_1->Eval(energy)+gXS_2->Eval(energy));
		double BRn1=gXS_1->Eval(energy)/(gXS_0->Eval(energy)+gXS_1->Eval(energy)+gXS_2->Eval(energy));
		double BRn2=gXS_2->Eval(energy)/(gXS_0->Eval(energy)+gXS_1->Eval(energy)+gXS_2->Eval(energy));
		excited=false;
		excited2=false;
		if(choose_chan<BRn0) {
			excited=false;
			excited2=false;
		}
		else if(choose_chan<BRn0+BRn1) {
			excited=true;
			excited2=false;
		}
		else {
			excited=false;
			excited2=true;
		}
		if(XS>XSmax) {G4cout<<"Max XS is exceeded!!"<<G4endl;}//
		if(XS>XSmax*G4UniformRand()) {
			return energy;
		}
		loopcounter++;
		if(loopcounter>1000) {
			G4cout<<"Rejection sampling not working properly - 1000 samples taken"<<G4endl;
			return 0;
		}
	}
	return 0;
}

double B1PrimaryGeneratorAction::Eloss(double Ea, double z) {//Energy of alpha of initial energy Ea after z mm in Be
	double range=0.000467*Ea*Ea+0.00238*Ea+0.000481;//mm
	if(z>range) return 0;//Fully stopped before Z
	range-=z;
	double Eout=-3232.85*range*range+282.754*range+0.0932;//MeV
	return Eout;
}

double B1PrimaryGeneratorAction::InteractionTheta(double Ep) {//COM theta
//	return (100+10*G4UniformRand())*d2r;
	double random_theta = acos(-1+2*G4UniformRand());
	random_theta = 2.*pi*G4UniformRand();
	double A[3]={1,0,0};
	bool sample=true;
	double XSmax=1.2;
	int counter=0;
	TGraph2D *select;
	if(!excited && !excited2) {
		select=g_gs;
	}
	if(excited && !excited2) {
		select=gp;
	}
	if(!excited && excited2) {
		select=gpp;
	}
	while(sample) {
		double XS_sample = XSmax*G4UniformRand();
		random_theta = acos(-1+2.*G4UniformRand());
		double XS=select->Interpolate(Ep,random_theta/d2r);
		if(XS>XSmax) G4cout<<"Sampling error for angle"<<G4endl;
		if(XS==0 || excited2) XS=1;//Outside of the energy range for this data set or Hoyle (poor data)
		if(XS_sample<XS) {
			return random_theta;//Isotropic for now
		}
		counter++;
		if(counter>100) {
			G4cout<<"Over 100 samples taken - no solution for Ealpha = "<<Ep<<"\t"<<excited<<"\t"<<excited2<<G4endl;
			return 0;
		}
	}
}

double B1PrimaryGeneratorAction::ConvertCMtoLab(double Ep, double thetaCM) {//Convert thetaCM to lab
	double Ex=0;
	if(excited) Ex=4.44;
	if(excited2) Ex=7.654;

	double gamma = sqrt(((mb*ml)/(mt*mh))*(Ep/(Ep+(Q-Ex)*(1.+mb/mt))));//Calculate gamma for the conversion to lab COM
	gamma = gamma*((mb+mt)/(ml+mh));//Account for change in the CM velocity
	double thetalab = atan2(sin(thetaCM),(cos(thetaCM)+gamma));//
	secondks = false;
	if(gamma*cos(thetaCM)<=-1.) secondks=true;//Neutron is going bac
//	thetalab=(12*G4UniformRand())*d2r;//TEST
//	if(G4UniformRand()>0.5) secondks=true;//TEST
	return thetalab;
}

double B1PrimaryGeneratorAction::NeutronEnergy(double Eb, double theta, double thetaCM) {//Get neutron energy for a given Ep and theta
	double En;
	double Ex=0;
	if(excited) Ex=4.44;//MeV
	if(excited2) Ex=7.654;

	double ECM = mt*Eb/(mt+mb)+Q-Ex;
	if(ECM<0) {
		G4cout<<"Below threshold! Something is wrong!"<<G4endl;
		return 0;
	}
	double p_nT = sqrt(2.*ECM*ml*(mh/(ml+mh)));
	double phi=2*pi*G4UniformRand();
	double pn[3]={p_nT*sin(thetaCM)*cos(phi),p_nT*sin(thetaCM)*sin(phi),p_nT*cos(thetaCM)};//COM neutron momentum
	pn[2]+=(ml/(mh+ml))*sqrt(2.*mb*Eb);//Add the momentum boost for the lab frame
	double En_new = 0;
	for(int i=0;i<3;i++) En_new+=pn[i]*pn[i];
	En_new /= (2.*ml);
	En = En_new;
	theta = acos(pn[2]/p_nT);
	return En;
}
