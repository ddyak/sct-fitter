#include "gtest/gtest.h"

#include "ThreeVector.h"
#include "dataobjects/Particle.h"
#include "treefitter/FitParams.h"
#include "treefitter/DecayChain.h"
#include "treefitter/ConstraintConfiguration.h"

#include <iostream>

using namespace sct::ana;

ParticlePtr generateDecay() {
    ParticlePtr pip = std::make_shared<Particle>(sct::kine::FourVector({ 1, 1, 1, 4 }), 211);
    ParticlePtr pim = std::make_shared<Particle>(sct::kine::FourVector({ 1, 1, 1, 4 }), -211);
    ParticlePtr mother = std::make_shared<Particle>(std::vector{ pim, pip }, 310); // ks
    return mother;
}


TEST(DecayChain, Dummy) {
    auto mother = generateDecay();
    DecayChain decay_chain(mother, {});
    
    EXPECT_EQ(decay_chain.dim(), 12);
}
