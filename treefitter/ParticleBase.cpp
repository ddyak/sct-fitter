#include "ParticleBase.h"


#include "ConstraintConfiguration.h"
#include "InternalParticle.h"
#include "DummyParticle.h"


using namespace sct::ana;


ParticleBase::ParticleBase(ParticlePtr particle, const ParticleBase* mother, const ConstraintConfiguration& config):
    m_particle(particle),
    m_mother(mother),
    m_config(config),
    m_index(0) 
{}

/* static */ ParticleBase* ParticleBase::createParticle(ParticlePtr particle, const ParticleBase* mother, const ConstraintConfiguration& config)
{
    return new DummyParticle(particle, mother, config);
}

void ParticleBase::updateIndex(int& offset) {
    for (auto* daughter : m_daughters) {
        daughter->updateIndex(offset);
    }
    m_index = offset;
    offset += dim();
}