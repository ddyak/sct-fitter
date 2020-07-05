#pragma once

#include "ParticleBase.h"

namespace sct::ana {

class InternalParticle : public ParticleBase {
    /**  get dimension of constraint */
    virtual int dim() const { return 4; }
};

}