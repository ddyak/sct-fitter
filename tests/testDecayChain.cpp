#include "gtest/gtest.h"

#include "ThreeVector.h"
#include "dataobjects/Particle.h"
#include "treefitter/include/FitParams.h"
#include "treefitter/include/DecayChain.h"
#include "treefitter/include/ConstraintConfiguration.h"

using namespace sct::ana;

ParticlePtr generateDecay() {
    ParticlePtr pip = std::make_shared<Particle>(sct::kine::FourVector({ 1, 2, 3, 4 }), 211);
    ParticlePtr pim = std::make_shared<Particle>(sct::kine::FourVector({ -1, -2, -3, 4 }), -211);
    ParticlePtr mother = std::make_shared<Particle>(std::vector{ pim, pip }, 310); // ks
    return mother;
}

TEST(DecayChain, Initialize) {
    auto mother = generateDecay();
    DecayChain decay_chain(mother, {});
    EXPECT_EQ(decay_chain.dim(), 12);

    FitParams fitparams(decay_chain.dim());
    decay_chain.initialize(fitparams);
    Eigen::VectorXd target_state_vector(12);
    target_state_vector << -1, -2, -3, 4, 1, 2, 3, 4, 0, 0, 0, 8;
    EXPECT_EQ(fitparams.getStateVector(), target_state_vector);
    EXPECT_EQ(fitparams.getCovariance().rows(), 12);
    EXPECT_EQ(fitparams.testCovariance(), true);

    decay_chain.filter(fitparams);
}
