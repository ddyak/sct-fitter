#pragma once

#include "dataobjects/Particle.h"

namespace sct::ana {

class ConstraintConfiguration;

class ParticleBase {
public:
    /** default constructor  */
    ParticleBase(ParticlePtr particle, const ParticleBase* mother, const ConstraintConfiguration& config);

    /** create the according treeFitter particle obj for a basf2 particle type  */
    static ParticleBase* createParticle(ParticlePtr particle, const ParticleBase* mother, const ConstraintConfiguration& config);

    /**  get dimension of constraint */
    virtual int dim() const = 0;

    /** get index  */
    int index() const { return m_index; }

    /** this sets the index for momentum, position, etc. in the statevector  */
    virtual void updateIndex(int& offset);

    /** get vertex index (in statevector!) */
    virtual int posIndex() const { return -1; }

    /** get tau index */
    virtual int tauIndex() const { return -1; }

    /** get momentum index */
    virtual int momIndex() const { return -1; }

    /** add daughter  */
    ParticleBase* addDaughter(ParticlePtr particle, const ConstraintConfiguration& config) {
        auto newDaughter = ParticleBase::createParticle(particle, this, config);
        m_daughters.push_back(newDaughter);
        return m_daughters.back();
    }

protected:
    /** pointer to framework type  */
    ParticlePtr m_particle;

    /** motherparticle */
    const ParticleBase* m_mother;

    /** daughter container */
    std::vector<ParticleBase*> m_daughters;

    /** has all the constraint config */
    const ConstraintConfiguration& m_config;

private:
    /** index */
    int m_index;
};

}