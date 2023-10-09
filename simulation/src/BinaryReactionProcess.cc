#include "BinaryReactionProcess.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"

BinaryReactionProcess::BinaryReactionProcess(const G4String& processName)
  : G4VDiscreteProcess(processName,fHadronic), fScatteringEnergy(1e6) {
  SetProcessSubType(111);
  fQValue = 0.;
}

BinaryReactionProcess::~BinaryReactionProcess() {
}

G4double BinaryReactionProcess::GetMeanFreePath(const G4Track& aTrack,
					    G4double /*previousStepSize*/,
					    G4ForceCondition* condition) {

  G4double energy = aTrack.GetKineticEnergy()/MeV;
  G4double mfp = DBL_MAX;//By default the MFP is 'infinite'
//  G4DynamicParticle* excitedstate = new G4DynamicParticle;
  //excitedstate = aTrack.GetDynamicParticle();

  G4String excitedname= aTrack.GetDynamicParticle()->GetDefinition()->GetParticleName();//Name of the particle, i.e. 'alpha','gamma' etc.
//  if(excitedname=="alpha" && energy<0.001) G4cout<<aTrack.GetPosition()<<G4endl;

//  G4cout<<"Excitedname: "<<excitedname<<"\t"<<mfp<<"\t"<<aTrack.GetTrackID()<<"\t"<<energy<<G4endl;
  *condition = NotForced;
/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Energy of the alpha is in variable energy, one can therefore now set the mfp (mean free path)      ///
// by pulling up a XS for the given energy and then converting to mean free path (via the density etc)///
/////////////////////////////////////////////////////////////////////////////////////////////////////////

// For now, I set the mfp = 20 cm
  if(excitedname=="alpha" && aTrack.GetTrackID()==1) mfp=20.*cm;  
//  G4cout<<aTrack.GetNextVolume()->GetLogicalVolume()->GetName()<<"\t"<<excitedname<<"\t"<<aTrack.GetParentID()<<G4endl;
  if(excitedname=="neutron"  && aTrack.GetParentID()==0 && aTrack.GetNextVolume()->GetLogicalVolume()->GetName()=="World") {
	mfp=50.*cm; 
  }
  if(excitedname=="deuteron" && aTrack.GetTrackID()==1 && aTrack.GetNextVolume()->GetLogicalVolume()->GetName()=="fGasCellLogical") {
	mfp=100.*cm; 
//        G4cout<<"Ready"<<G4endl;
  }
//  if(excitedname=="proton") G4cout<<aTrack.GetNextVolume()->GetLogicalVolume()->GetName()<<"\t"<<energy<<G4endl;
//  G4cout<<neutronbeam<<"\t"<<mfp<<G4endl;
  return mfp;
}

G4VParticleChange* BinaryReactionProcess::PostStepDoIt( const G4Track& aTrack,
						    const G4Step& aStep) {
  
  G4StepPoint* postStepPoint = aStep.GetPostStepPoint();
  if (postStepPoint->GetStepStatus()==fGeomBoundary) {
    return G4VDiscreteProcess::PostStepDoIt(aTrack,aStep);
  }
  // Do the TwoBody decay -- set to a+12C->n+15O at the moment
  aParticleChange.Initialize(aTrack);
  if(aTrack.GetDynamicParticle()->GetDefinition()->GetParticleName()=="alpha") TwoBody(aTrack,aStep,/*Target charge*/6,/*Target mass*/12,/*Light product charge*/0,/*Light product mass*/1,/*Heavy product charge*/8,/*Heavy product mass*/15,/*Excitation energy of light particle*/0*MeV,/*Excitation of heavy product*/0*MeV);
  if(aTrack.GetDynamicParticle()->GetDefinition()->GetParticleName()=="deuteron") TwoBody(aTrack,aStep,1,2,0,1,2,3,0*MeV,0*MeV);
  if(aTrack.GetDynamicParticle()->GetDefinition()->GetParticleName()=="neutron") {
//	G4cout<<"Neutron elastic"<<G4endl;
//	TwoBody(aTrack,aStep,1,1,0,1,1,1,0*MeV,0*MeV);
	TwoBody(aTrack,aStep,6,12,2,4,4,9,0*MeV,0*MeV);
  }
  return &aParticleChange;
}

void BinaryReactionProcess::StartTracking(G4Track* track) {
  G4VProcess::StartTracking(track);	// Apply base class actions
  fScatteringEnergy = track->GetKineticEnergy()*G4UniformRand()/MeV;
}

void BinaryReactionProcess::ParseParams(std::map<std::string,double> &params) {
    double lightProductMass = -1;
    double lightProductCharge = -1;
    double heavyProductMass = -1;
    double heavyProductCharge = -1;
    double targetMass = -1;
    double targetCharge = -1;

    if(lightProductCharge>0 && lightProductMass>0) {
      SetLightProduct(lightProductCharge,lightProductMass);
    }
    
    if(lightProductCharge==0 && lightProductMass==1) {//added for neutron
      SetLightProduct(lightProductCharge,lightProductMass);
    }
    if(heavyProductCharge>0 && heavyProductMass>0) {
      SetHeavyProduct(heavyProductCharge,heavyProductMass);
    }
    if(targetCharge>0 && targetMass>0) {
      SetTarget(targetCharge,targetMass);
    }
}


G4VParticleChange* BinaryReactionProcess::TwoBody( const G4Track& aTrack, const G4Step& aStep, int Zt, int At, int Z1, int A1, int Z2, int A2, double Ex1, double Ex2) {
//      JEB version of binary reaction mechanism for beam + (Zt,At) -> (A1,Z1,Ex1) + (A2,Z2,Ex2)
        G4double Mt,M1,M2;
        G4double energy=aTrack.GetDynamicParticle()->GetKineticEnergy()/MeV;
        G4DynamicParticle* target = new G4DynamicParticle;
        G4ParticleDefinition* targetdef;
        targetdef=G4IonTable::GetIonTable()->GetIon(Zt,At,0.);
        target->SetDefinition(targetdef);
        Mt=targetdef->GetPDGMass()/CLHEP::amu_c2;
        G4DynamicParticle* part1 = new G4DynamicParticle;        G4ParticleDefinition* part1def;
        if(Z1==0) {
                G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
                G4String particleName;
                part1def = particleTable->FindParticle(particleName="neutron");
        }
        if(Z1!=0) {
                part1def=G4IonTable::GetIonTable()->GetIon(Z1,A1,Ex1*MeV);
        }
        part1->SetDefinition(part1def);
        M1=part1def->GetPDGMass()/CLHEP::amu_c2;
        G4DynamicParticle* part2 = new G4DynamicParticle;
        G4ParticleDefinition* part2def;
        part2def=G4IonTable::GetIonTable()->GetIon(Z2,A2,Ex2*MeV);
        part2->SetDefinition(part2def);
        M2=part2def->GetPDGMass()/CLHEP::amu_c2;
////Do some kinematics
        G4double CM_theta=acos(-1.+2.*G4UniformRand());//0->pi relative to the particle direction
	if(aTrack.GetDynamicParticle()->GetDefinition()->GetParticleName()=="deuteron")CM_theta=0.;
        G4double CM_thetah=4.*atan(1.)-CM_theta;//pi-cm_theta_light
        G4double CM_psi=G4UniformRand()*4.*2.*atan(1.);//0->2pi
        G4double CM_psih=8.*atan(1.)-CM_psi;//2pi-cm_psi_light
///Rotate so theta/psi are relative to the original beam direction
        G4ThreeVector momentumDirection = aTrack.GetMomentumDirection();//using rotation method from binaryreactionphysics
        G4ThreeVector v = G4ThreeVector(0.,0.,1.).cross(momentumDirection);
        G4double rotAngle = acos(momentumDirection.z());

        G4ThreeVector dir=G4ThreeVector(sin(CM_theta)*sin(CM_psi),sin(CM_theta)*cos(CM_psi),cos(CM_theta));
        if(v.getR()>0) dir.rotate(v,rotAngle);//rotate the direction to be relative to the beam axis


//Get Q-value
//      G4cout<<CM_theta<<G4endl;
        G4double Q_value=0;
        Q_value=aTrack.GetDynamicParticle()->GetDefinition()->GetPDGMass()+targetdef->GetPDGMass()-(part1def->GetPDGMass()+part2def->GetPDGMass());
//        G4cout<<"Effective Q-value= "<<Q_value/MeV<<G4endl;
//      G4cout<<"E_CM = "<<Mt*energy/(M1+M2)<<G4endl;
        G4double E_cm = (Mt*energy/(M1+M2))+Q_value;//new CM energy MeV
        if(E_cm<0.) {
		G4cout<<"Below threshold"<<G4endl;
		return &aParticleChange;//sub-threshold
	}
        G4double p_1 = sqrt(2.*part1->GetMass()*E_cm*(1.*M2/(M1+M2)));//E_1 = m2/(m1+m2) * E_t
//        G4double p_2 = sqrt(2.*part2->GetMass()*E_cm*(1.*M1/(M1+M2)));//E_2 = m1/(m1+m2) * E_t
//      G4cout<<"Free momentum: "<<p_1<<"\t"<<p_2<<"\tE_CM\t"<<E_cm<<G4endl;
        G4ThreeVector p_new_1 = p_1*dir;// new momentum of scattered Be in COM
        G4ThreeVector p_new_2 = -p_new_1;
        G4ThreeVector p_n = aTrack.GetMomentum();
        p_new_1+=p_n*(1.*M1/(M1+M2));
        p_new_2+=p_n*(1.*M2/(M1+M2));
//      G4cout<<"Part 1:\t"<<p_new_1<<G4endl;
//      G4cout<<"Part 2:\t"<<p_new_2<<G4endl;
//      G4cout<<"Orig:\t"<<p_n<<G4endl;
        part1->SetMomentum(p_new_1);
        part2->SetMomentum(p_new_2);
        G4double total_mom_1=p_new_1.getR();
        G4double total_mom_2=p_new_2.getR();
        G4double energy1=((total_mom_1*total_mom_1)/(2.*part1->GetMass()));
        G4double energy2=((total_mom_2*total_mom_2)/(2.*part2->GetMass()));
        part1->SetKineticEnergy(energy1);
        part2->SetKineticEnergy(energy2);
        G4double lab_theta=p_new_1.theta();
        G4double lab_psi=p_new_1.phi();
        G4double lab_thetah=p_new_2.theta();
        G4double lab_psih=p_new_2.phi();
//      aTrack.GetDynamicParticle()->DumpInfo();
//      G4cout<<"Alpha E: "<<(total_mom_1*total_mom_1)/(2.*part1->GetMass())<<"\tMom. vector: "<<part1->GetMomentumDirection()<<G4endl;
//        G4cout<<"Dumping info - comment out at end of src/BinaryReactionProcess.cc to stop"<<G4endl;

//      part1->DumpInfo();
  //    part2->DumpInfo();
	if(Z1==0 && Z2==2) G4cout<<"NeutronE "<<energy1<<G4endl;
	if(Z2==1) G4cout<<"ScatteredE "<<energy2<<G4endl;
        G4Track* sec1 = new G4Track(part1,
                aTrack.GetGlobalTime(),
                aTrack.GetPosition());
        G4double projectileMass = aTrack.GetDefinition()->GetAtomicMass();
        G4Track* sec2 = new G4Track(part2,
                aTrack.GetGlobalTime(),
                aTrack.GetPosition());

        aParticleChange.AddSecondary(sec1);
        aParticleChange.AddSecondary(sec2);
        aParticleChange.ProposeEnergy(0.);
        aParticleChange.ProposeTrackStatus(fStopAndKill);
        return &aParticleChange;

}

