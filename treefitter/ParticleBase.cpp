#include "ParticleBase.h"
#include "InternalParticle.h"

using namespace sct::ana;


/* static */ ParticleBase* ParticleBase::createParticle(ParticlePtr particle, const ParticleBase* mother)
{
    return new InternalParticle();
}

void ParticleBase::updateIndex(int& offset) {
    for (auto* daughter : m_daughters) {
        daughter->updateIndex(offset);
    }
    m_index = offset;
    offset += dim();
}