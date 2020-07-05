#include <iostream>
#include "ThreeVector.h"
#include "dataobjects/Particle.h"

using namespace sct::ana;

// todo: read from json
// https://github.com/VitalyVorobyev/decay-gen/tree/Generator

ParticlePtr generateDecay() {
    ParticlePtr pip = std::make_shared<Particle>(sct::kine::FourVector({ 1, 1, 1, 4 }), 211);
    ParticlePtr pim = std::make_shared<Particle>(sct::kine::FourVector({ 1, 1, 1, 4 }), -211);
    ParticlePtr mother = std::make_shared<Particle>( std::vector{ pim, pip }, 310); // ks
    return mother;
}

int main() {
    auto mother = generateDecay();
    std::cerr << mother->momentum() << std::endl;
    return 0;
}