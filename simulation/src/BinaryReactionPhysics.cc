#include "BinaryReactionPhysics.hh"

#include "BinaryReactionProcess.hh"
#include "G4GenericIon.hh"
#include "G4Neutron.hh"
#include "G4Proton.hh"
#include "G4Deuteron.hh"
#include "G4Alpha.hh"
#include "globals.hh"
#include "G4PhysicsListHelper.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(BinaryReactionPhysics);

BinaryReactionPhysics::BinaryReactionPhysics(G4int)
:  G4VPhysicsConstructor("BinaryReactionPhysics") {
}

BinaryReactionPhysics::BinaryReactionPhysics(const G4String& name)
:  G4VPhysicsConstructor(name) {
    G4cout<<"Creating Binary Reaction Physics"<<G4endl;
}

BinaryReactionPhysics::~BinaryReactionPhysics()
{
}

void BinaryReactionPhysics::ConstructParticle()
{
  G4GenericIon::GenericIon();
  G4Neutron::Neutron();
  G4Proton::Proton();
  G4Alpha::Alpha();
}

void BinaryReactionPhysics::ConstructProcess()
{
  BinaryReactionProcess* reactionProcess = new BinaryReactionProcess();
  reactionProcess->ParseParams(fReactionParams);
  G4PhysicsListHelper::GetPhysicsListHelper()->
    RegisterProcess(reactionProcess, G4GenericIon::GenericIon());
  G4PhysicsListHelper::GetPhysicsListHelper()->
    RegisterProcess(reactionProcess, G4Neutron::Definition());
  G4PhysicsListHelper::GetPhysicsListHelper()->
    RegisterProcess(reactionProcess, G4Proton::Definition());
  G4PhysicsListHelper::GetPhysicsListHelper()->
    RegisterProcess(reactionProcess, G4Deuteron::Definition());
  G4PhysicsListHelper::GetPhysicsListHelper()->
    RegisterProcess(reactionProcess, G4Alpha::Definition());


}

