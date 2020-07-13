#include "dataobjects/Particle.h"


namespace sct::ana {

class ConstraintConfiguration;
class FitParams;
class DecayChain;
class ParticleBase;

/** class for fit and work with configuration */
class FitManager {
public:
    FitManager(ParticlePtr particle, const ConstraintConfiguration& config);
    
    /** main fit function */
    bool fit();

    /** update particles parameters with the fit results */
    bool updateCand(ParticlePtr particle, const bool isTreeHead) const;

    /** locate particle base for a belle2 particle and update the particle with the values from particle base */
    void updateCand(const ParticleBase& pb, ParticlePtr cand, const bool isTreeHead) const;

    /** update the Belle2::Particles with the fit results  */
    void updateTree(ParticlePtr particle, const bool isTreeHead) const;

    /** getter for chi2 of the newton iteration */
    double chiSquare() const { return m_chiSquare; }

private:
    /** head of tree */
    ParticlePtr m_particle;

    const ConstraintConfiguration& m_config;
    bool m_useReferencing = false;
    bool m_updateDaugthers = true;

    double m_chiSquare;

    FitParams* m_fitparams;
    DecayChain* m_decaychain;
};

}
