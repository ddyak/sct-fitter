#include "dataobjects/Particle.h"


namespace sct::ana {

class ConstraintConfiguration;
class FitParams;
class DecayChain;

/** class for fit and work with configuration */
class FitManager {
public:
    FitManager(ParticlePtr particle, const ConstraintConfiguration& config);
    
    /** main fit function */
    bool fit() { return false; };

    double chiSquare() const { return m_chiSquare; }

private:
    /** head of tree */
    ParticlePtr m_particle;

    const ConstraintConfiguration& m_config;

    double m_chiSquare;

    FitParams* m_fitparams;
    DecayChain* m_decaychain;
};

}
