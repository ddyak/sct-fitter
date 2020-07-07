#pragma once

#include "ParticleBase.h"
#include "FitParams.h"

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

    virtual bool initParticle(FitParams& fitparams) override {
        for (auto daughter : m_daughters) {
            daughter->initParticle(fitparams);
        }

        int momindex = momIndex();
        fitparams.getStateVector().segment(momindex, 4) = Eigen::Matrix<double, 4, 1>::Zero(4);

        // temp
        if (m_daughters.size() == 0) {
            fitparams.getStateVector().segment(momindex, 4) = m_particle->fourMomentum();
            return true;
        }

        for (auto daughter : m_daughters) {
            int daumomindex = daughter->momIndex();
            int maxrow = 4;//daughter->hasEnergy() ? 4 : 3;

            double e2 = fitparams.getStateVector().segment(daumomindex, maxrow).squaredNorm();
            fitparams.getStateVector().segment(momindex, maxrow) += fitparams.getStateVector().segment(daumomindex, maxrow);

            // if (maxrow == 3) {
            //     double mass = daughter->pdgMass();
            //     fitparams.getStateVector()(momindex + 3) += std::sqrt(e2 + mass * mass);
            // }
        }
        
        return true;
    }

};

}