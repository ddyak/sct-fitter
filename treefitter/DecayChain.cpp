#include "DecayChain.h"

#include "ParticleBase.h"

using namespace sct::ana;

DecayChain::DecayChain(ParticlePtr particle, const ConstraintConfiguration& config):
    m_config(config),
    m_dim(0)
{
    m_headOfChain = ParticleBase::createParticle(particle, nullptr, config);
    m_headOfChain->updateIndex(m_dim);
}

bool DecayChain::initialize(FitParams& par) {
    return m_headOfChain->initParticle(par);
}