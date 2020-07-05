#pragma once

#include "dataobjects/Particle.h"

namespace sct::ana {

class ParticleBase {
public:
    /** create the according treeFitter particle obj for a basf2 particle type  */
    static ParticleBase* createParticle(ParticlePtr particle, const ParticleBase* mother);

    /**  get dimension of constraint */
    virtual int dim() const = 0;

    /** this sets the index for momentum, position, etc. in the statevector  */
    virtual void updateIndex(int& offset);

private:
    /** index */
    int m_index;

    /** motherparticle */
    const ParticleBase* m_mother;

    /** daughter container  */
    std::vector<ParticleBase*> m_daughters;
};

}