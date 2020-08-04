#include "InternalParticle.h"

#include "ParticleBase.h"
#include "FitParams.h"
#include "Projection.h"

using namespace sct::ana;


InternalParticle::InternalParticle(ParticlePtr particle, const ParticleBase* mother, const ConstraintConfiguration& config):
    ParticleBase(particle, mother, config)
{
    for (auto daughter : particle->daughters()) {
        addDaughter(daughter, config);
    }
    // add flags from config
}

int InternalParticle::dim() const
{
    // { x, y, z, tau, px, py, pz, E }
    // the last 4 always exists for composite particles
    // tau index conly exist with a vertex and geo constraint
    //
    if (m_onlymomentum) { return 4; }
    if (m_shares_vertex_with_mother) { return 4; }
    if (!m_shares_vertex_with_mother && !m_geo_constraint) { return 7; }
    if (!m_shares_vertex_with_mother && m_geo_constraint) { return 8; }

    // this case should not appear
    if (m_shares_vertex_with_mother && m_geo_constraint) { return -1; }
    return -1;
}

int InternalParticle::posIndex() const
{
    if (m_onlymomentum) return -1;

    // for example B0 and D* can share the same vertex
    return m_shares_vertex_with_mother ? this->mother()->posIndex() : index();
}

int InternalParticle::tauIndex() const
{
    if (m_onlymomentum) return -1;
    /** only exists if particle is geo cosntraint and has a mother */
    return m_geo_constraint ? index() + 3 : -1;
}

int InternalParticle::momIndex() const
{
    /** indexing in { x, y, z, tau, px, py, pz, E }
     * but tau is not existing for all InternalParticles
     * */

    if (m_onlymomentum) { return this->index(); };

    if (m_geo_constraint && !m_shares_vertex_with_mother) { return this->index() + 4; }

    if (m_shares_vertex_with_mother) { return this->index(); }

    if (!m_shares_vertex_with_mother && !m_geo_constraint) { return index() + 3; }

    // this will crash the initialisation
    return -1;
}

ErrCode InternalParticle::initCovariance(FitParams& fitparams) const {
    ParticleBase::initCovariance(fitparams);
    for (auto daughter : m_daughters) {
        daughter->initCovariance(fitparams);
    }
    return ErrCode(ErrCode::Status::success);
}

ErrCode InternalParticle::initMomentum(FitParams& fitparams) const
{
    int momindex = momIndex();
    fitparams.getStateVector().segment(momindex, 4) = Eigen::Matrix<double, 4, 1>::Zero(4);

    for (auto daughter : m_daughters) {
        int daumomindex = daughter->momIndex();
        int maxrow = daughter->hasEnergy() ? 4 : 3;

        fitparams.getStateVector().segment(momindex, maxrow) += fitparams.getStateVector().segment(daumomindex, maxrow);

        if (maxrow == 3) {
            double e2 = fitparams.getStateVector().segment(daumomindex, maxrow).squaredNorm();
            double mass = daughter->pdgMass();
            fitparams.getStateVector()(momindex + 3) += std::sqrt(e2 + mass * mass);
        }
    }
    return ErrCode(ErrCode::Status::success);
}

ErrCode InternalParticle::projectConstraint(const Constraint::Type type,
                                            const FitParams& fitparams,
                                            Projection& p) const
{
    ErrCode status;
    switch (type) {
        case Constraint::mass:
        status |= ParticleBase::projectMassConstraintParticle(fitparams, p);
        break;
        case Constraint::geometric:
        if (!m_onlymomentum) 
            status |= ParticleBase::projectGeoConstraint(fitparams, p);
        break;
        case Constraint::kinematic:
        status |= projectKineConstraint(fitparams, p);
        break;
        default:
        throw std::runtime_error("Internal: project Constraint unknown");
    }

    return status;
}

ErrCode InternalParticle::projectKineConstraint(const FitParams& fitparams, Projection& p) const {
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
    return ErrCode(ErrCode::Status::success);
}

void InternalParticle::addToConstraintList(std::vector<Constraint>& list, int depth) const {
    for (auto daughter : m_daughters) {
        daughter->addToConstraintList(list, depth - 1);
    }
    // if (tauIndex() >= 0 && m_lifetimeconstraint) {
        // list.push_back(Constraint(this, Constraint::lifetime, depth, 1));
    // }
    if (momIndex() >= 0) {
        list.push_back(Constraint(this, Constraint::kinematic, depth, 4, 3));
    }
    if (m_geo_constraint && !m_onlymomentum) {
        // assert(m_config);
        const int dim = 3;// m_config->m_originDimension == 2 && std::abs(m_particle->getPDGCode()) == m_config->m_headOfTreePDG ? 2 : 3;
        list.push_back(Constraint(this, Constraint::geometric, depth, dim, 3));
    }
    if (m_massconstraint) {
        list.push_back(Constraint(this, Constraint::mass, depth, 1, 3));
    }
}

ErrCode InternalParticle::initParticle(FitParams& fitparams) {
    for (auto daughter : m_daughters) {
        daughter->initParticle(fitparams);
    }

    int momindex = momIndex();
    fitparams.getStateVector().segment(momindex, 4) = Eigen::Matrix<double, 4, 1>::Zero(4);

    initMomentum(fitparams);
    
    return ErrCode(ErrCode::Status::success);
}