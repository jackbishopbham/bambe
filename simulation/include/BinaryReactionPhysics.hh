#ifndef BinaryReactionPhysics_h
#define BinaryReactionPhysics_h

#include "G4VPhysicsConstructor.hh"
#include <string>
#include <map>

class BinaryReactionPhysics : public G4VPhysicsConstructor {
public:
  BinaryReactionPhysics(G4int verbose =1);
  BinaryReactionPhysics(const G4String& name);
  virtual ~BinaryReactionPhysics();

  virtual void ConstructParticle();
  virtual void ConstructProcess();

  void SetReactionParams() {
    G4cout<<"One may set reaction parameters here in include/BinaryReactionPhysics.hh"<<G4endl;
  };

private:
  std::map<std::string,double> fReactionParams;
};

#endif
