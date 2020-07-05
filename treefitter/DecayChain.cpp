#include "DecayChain.h"

#include "ParticleBase.h"

using namespace sct::ana;

DecayChain::DecayChain(ParticlePtr particle, const ConstraintConfiguration& config):
    m_config(config),
    m_dim(0)
{
    m_headOfChain = ParticleBase::createParticle(particle, nullptr);
    m_headOfChain->updateIndex(m_dim);
}
