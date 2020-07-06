#pragma once

#include "ParticleBase.h"

namespace sct::ana {

class InternalParticle : public ParticleBase {
public:
    /** space reserved in fit params, if has mother then it has tau */
    virtual int dim() const { return 4; }
    
};

}