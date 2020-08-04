#include <iostream>
#include "ThreeVector.h"
#include "dataobjects/Particle.h"

#include "ConstraintConfiguration.h"
#include "FitManager.h"

using namespace sct::ana;

// todo: read from json
// https://github.com/VitalyVorobyev/decay-gen/tree/Generator

ParticlePtr generateDecay() {
    ParticlePtr pip = std::make_shared<Particle>(sct::kine::ThreeVector<double>{-0.1, 0.0, -0.3}, 0.139, 211);
    ParticlePtr pim = std::make_shared<Particle>(sct::kine::ThreeVector<double>{0.1, 0.0, 0.1}, 0.139, -211);
    ParticlePtr mother = std::make_shared<Particle>(std::vector{ pim, pip }, 310); // ks
    return mother;
}

int main() {
    auto mother = generateDecay();
    std::cerr << mother->mass() << std::endl;
    std::cerr << mother->pdgMass() << std::endl;
    
    FitManager fitmanager(mother, {});
    fitmanager.fit();
    std::cerr << mother->momentum() << std::endl;
    return 0;
}