#pragma once

#include "ParticleBase.h"
#include "FitParams.h"
#include "Projection.h"

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

    // does the particle have a 3-momentum or a 4-momentum ?
    /** get momentum dimension */
    virtual bool hasEnergy() const { return true; }

    /** get false  */
    virtual bool hasPosition() const { return false; }

    /** project constraint */
    virtual bool projectConstraint(const Constraint::Type type, 
                                    const FitParams& fitparams,
                                    Projection& p) const override {
        bool status;
        switch (type) {
        case Constraint::mass:
            status |= projectMassConstraintParticle(fitparams, p);
            break;
        // case Constraint::geometric:
            // status |= projectGeoConstraint(fitparams, p);
            // break;
        case Constraint::kinematic:
            status |= projectKineConstraint(fitparams, p);
            break;
        // default:
            // status |= ParticleBase::projectConstraint(type, fitparams, p);
        }

        return true;
    };

private:
    /** project kinematical constraint */
    bool projectKineConstraint(const FitParams& fitparams, Projection& p) const {
        const int momindex = momIndex();

        // `this` always has an energy row
        p.getResiduals().segment(0, 4) = fitparams.getStateVector().segment(momindex, 4);

        for (int imom = 0; imom < 4; ++imom) {
            p.getH()(imom, momindex + imom) = 1;
        }

        for (const auto daughter : m_daughters) {
            const int daumomindex = daughter->momIndex();
            const Eigen::Matrix<double, 1, 3> p3_vec = fitparams.getStateVector().segment(daumomindex, 3);

            // three momentum is easy just substract the vectors
            p.getResiduals().segment(0, 3) -= p3_vec;

            // energy depends on the parametrisation!
            if (daughter->hasEnergy()) {
                p.getResiduals()(3) -= fitparams.getStateVector()(daumomindex + 3);
                p.getH()(3, daumomindex + 3) = -1; // d/dE -E
            } else {
                // m^2 + p^2 = E^2
                // so
                // E = sqrt(m^2 + p^2)
                const double mass = daughter->pdgMass();
                const double p2 = p3_vec.squaredNorm();
                const double energy = std::sqrt(mass * mass + p2);
                p.getResiduals()(3) -= energy;

                for (unsigned i = 0; i < 3; ++i) {
                // d/dpx_i sqrt(m^2 + p^2)
                p.getH()(3, daumomindex + i) = -1 * p3_vec(i) / energy;
                }
            }

            // this has to be in any case
            // d/dp_i p_i
            for (unsigned i = 0; i < 3; ++i) {
                p.getH()(i, daumomindex + i) = -1;
            }
        }
        return true;
    }


public:
    /** add to constraint list */
    virtual void addToConstraintList(std::vector<Constraint>& list, int depth) const {
        for (auto daughter : m_daughters) {
           daughter->addToConstraintList(list, depth - 1);
        }
        // if (tauIndex() >= 0 && m_lifetimeconstraint) {
        //     list.push_back(Constraint(this, Constraint::lifetime, depth, 1));
        // }
        if (momIndex() >= 0 && !m_daughters.empty()) {
            list.push_back(Constraint(this, Constraint::kinematic, depth, 4, 3));
        }
        // if (m_geo_constraint) {
        //     assert(m_config);
        //     const int dim = m_config->m_originDimension == 2 && std::abs(m_particle->getPDGCode()) == m_config->m_headOfTreePDG ? 2 : 3;
        //     list.push_back(Constraint(this, Constraint::geometric, depth, dim, 3));
        // }
        if (m_massconstraint) {
            list.push_back(Constraint(this, Constraint::mass, depth, 1, 3));
        }
    };

protected:
    /** init momentum of *this and daughters */
    bool initMomentum(FitParams& fitparams) const { 
        // TODO: move code here from initParticle  
        return false; 
    }

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

    bool initCovariance(FitParams& fitparams) const override {
        ParticleBase::initCovariance(fitparams);
        for (auto daughter : m_daughters) {
            daughter->initCovariance(fitparams);
        }
        return true;
    }

private:
    bool m_massconstraint = true;

};


}