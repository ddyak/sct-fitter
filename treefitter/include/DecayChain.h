#pragma once

#include "dataobjects/Particle.h"
#include "ErrCode.h"
#include "Constraint.h"
#include <map>

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
    ErrCode initialize(FitParams& par);

    /** filter down the chain */
    ErrCode filter(FitParams& par);

    /** filter with respect to a previous iteration for better stability */
    ErrCode filterWithReference(FitParams& par, const FitParams& ref);

    /** convert Belle2::particle into particle base(fitter base particle) */
    const ParticleBase* locate(ParticlePtr bc) const;

    /** get dimension */
    int dim() const { return m_dim; }

    /** init contraintlist */
    void initConstraintList();

private:
    /** the dimension of constraint */
    mutable int m_dim; 
    
    /** head of decay tree*/
    ParticleBase* m_headOfChain;

    /** list of constraints */
    std::vector<Constraint> m_constraintlist;
    
    /** config container */
    const ConstraintConfiguration& m_config;

    /** the map from Belle2::Particles to TreeFitter::ParticleBase */
    std::map<ParticlePtr, const ParticleBase*> m_particleMap;
};

}
