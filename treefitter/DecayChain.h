#pragma once
#include "dataobjects/Particle.h"

namespace sct::ana {

class ConstraintConfiguration;
class FitParams;
class ParticleBase;

/** this class does a lot of stuff:
 Build decaytree structure allowing to index particles and handle the filtering of constraints across the tree
*/
class DecayChain {
public:
    /**  constructor   */
    DecayChain(ParticlePtr particle, const ConstraintConfiguration& config);

    /** initalize the chain */
    bool initialize(FitParams& par) { return false; }

    /** filter down the chain */
    bool filter(FitParams& par) { return false; }

    /** get dimension   */
    int dim() const { return m_dim;}

private:
    /** the dimension of constraint */
    mutable int m_dim; 

    /** config container */
    const ConstraintConfiguration& m_config;
    
    /** head of decay tree*/
    ParticleBase* m_headOfChain;
};

}
