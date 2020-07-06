#pragma once

#include "ParticleBase.h"

namespace sct::ana {

class FitParams;
class ConstraintConfiguration;

class DummyParticle : public ParticleBase {
public:
    DummyParticle(ParticlePtr particle, const ParticleBase* mother, const ConstraintConfiguration& config):
        ParticleBase(particle, mother, config)
    {
        for (auto daughter : particle->daughters()) {
            addDaughter(daughter, config);
        }

        // add flags from config
    }

    /**  get dimension of constraint */
    virtual int dim() const override { return 4; }

    /** get momentum index */
    virtual int momIndex() const { return index(); }

protected:
    /** init momentum of *this and daughters */
    bool initMomentum(FitParams& fitparams) const { return false; }
};

}