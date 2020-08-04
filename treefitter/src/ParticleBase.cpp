#include "ParticleBase.h"

#include "ConstraintConfiguration.h"
#include "DummyParticle.h"
#include "InternalParticle.h"
#include "Projection.h"
#include "FitParams.h"

using namespace sct::ana;


ParticleBase::ParticleBase(ParticlePtr particle, const ParticleBase* mother, const ConstraintConfiguration& config):
    m_particle(particle),
    m_mother(mother),
    m_config(config),
    m_index(0),
    m_pdgMass(particle->pdgMass())
{}

/* static */ ParticleBase* ParticleBase::createParticle(ParticlePtr particle, const ParticleBase* mother, const ConstraintConfiguration& config)
{
    if (particle->daughters().empty())
        return new DummyParticle(particle, mother, config);
    else
        return new InternalParticle(particle, mother, config);
}

void ParticleBase::updateIndex(int& offset) {
    for (auto* daughter : m_daughters) {
        daughter->updateIndex(offset);
    }
    m_index = offset;
    offset += dim();
}

ParticleBase* ParticleBase::addDaughter(ParticlePtr particle, const ConstraintConfiguration& config) {
    auto newDaughter = ParticleBase::createParticle(particle, this, config);
    m_daughters.push_back(newDaughter);
    return m_daughters.back();
}

ErrCode ParticleBase::initCovariance(FitParams& fitparams) const
{
    //this is very sensitive and can heavily affect the efficiency of the fit
    const int posindex = posIndex();
    if (posindex >= 0) {
        for (int i = 0; i < 3; ++i) {
            fitparams.getCovariance()(posindex + i, posindex + i) = 1;
        }
    }

    const int momindex = momIndex();
    if (momindex >= 0) {
        const int maxrow = hasEnergy() ? 4 : 3;
        for (int i = 0; i < maxrow; ++i) {
            fitparams.getCovariance()(momindex + i, momindex + i) = 10.;
        }
    }

    const int tauindex = tauIndex();
    if (tauindex >= 0) {
        fitparams.getCovariance()(tauindex, tauindex) = 1.;
    }

    return ErrCode(ErrCode::Status::success);
}

ErrCode ParticleBase::projectMassConstraintParticle(const FitParams& fitparams, Projection& p) const {
    const double mass = pdgMass();
    const double mass2 = mass * mass;
    const int momindex = momIndex();
    const double px = fitparams.getStateVector()(momindex);
    const double py = fitparams.getStateVector()(momindex + 1);
    const double pz = fitparams.getStateVector()(momindex + 2);
    const double E  = fitparams.getStateVector()(momindex + 3);

    /** be aware that the signs here are important
     * E-|p|-m extracts a negative mass and messes with the momentum !
     * */
    p.getResiduals()(0) = mass2 - E * E + px * px + py * py + pz * pz;

    p.getH()(0, momindex)     = 2.0 * px;
    p.getH()(0, momindex + 1) = 2.0 * py;
    p.getH()(0, momindex + 2) = 2.0 * pz;
    p.getH()(0, momindex + 3) = -2.0 * E;

    // TODO 0 in most cases -> needs special treatment if width=0 to not crash chi2 calculation
    // const double width = TDatabasePDG::Instance()->GetParticle(particle()->getPDGCode())->Width();
    // transport  measurement uncertainty into residual system
    // f' = sigma_x^2 * (df/dx)^2
    // p.getV()(0) = width * width * 4 * mass2;

    return ErrCode(ErrCode::Status::success);
}

const ParticleBase* ParticleBase::mother() const {
    const ParticleBase* rc = m_mother;
    // while (rc && rc->type() == kFeedthroughParticle) {
        // rc = rc->mother();
    // }
    return rc;
}

const ParticleBase* ParticleBase::locate(ParticlePtr particle) const {
    const ParticleBase* rc = (m_particle == particle) ? this : nullptr;
    if (!rc) {
        for (auto* daughter : m_daughters) {
            rc = daughter->locate(particle);
            if (rc) {break;}
        }
    }
    return rc;
}

ErrCode ParticleBase::projectGeoConstraint(const FitParams& fitparams, Projection& p) const
{
    //assert(m_config);
    // only allow 2d for head of tree particles that are beam constrained
    const int dim = 3;// m_config->m_originDimension == 2 && std::abs(m_particle->getPDGCode()) == m_config->m_headOfTreePDG ? 2 : 3;
    const int posindexmother = mother()->posIndex();
    const int posindex = posIndex();
    const int tauindex = tauIndex();
    const int momindex = momIndex();

    const double tau = fitparams.getStateVector()(tauindex);
    Eigen::Matrix < double, 1, -1, 1, 1, 3 > x_vec = fitparams.getStateVector().segment(posindex, dim);
    Eigen::Matrix < double, 1, -1, 1, 1, 3 > x_m = fitparams.getStateVector().segment(posindexmother, dim);
    Eigen::Matrix < double, 1, -1, 1, 1, 3 > p_vec = fitparams.getStateVector().segment(momindex, dim);
    const double mom = p_vec.norm();
    const double mom3 = mom * mom * mom;

    if (3 == dim) {
        // we can already set these
        //diagonal momentum
        p.getH()(0, momindex)     = tau * (p_vec(1) * p_vec(1) + p_vec(2) * p_vec(2)) / mom3 ;
        p.getH()(1, momindex + 1) = tau * (p_vec(0) * p_vec(0) + p_vec(2) * p_vec(2)) / mom3 ;
        p.getH()(2, momindex + 2) = tau * (p_vec(0) * p_vec(0) + p_vec(1) * p_vec(1)) / mom3 ;

        //offdiagonal momentum
        p.getH()(0, momindex + 1) = - tau * p_vec(0) * p_vec(1) / mom3 ;
        p.getH()(0, momindex + 2) = - tau * p_vec(0) * p_vec(2) / mom3 ;

        p.getH()(1, momindex + 0) = - tau * p_vec(1) * p_vec(0) / mom3 ;
        p.getH()(1, momindex + 2) = - tau * p_vec(1) * p_vec(2) / mom3 ;

        p.getH()(2, momindex + 0) = - tau * p_vec(2) * p_vec(0) / mom3 ;
        p.getH()(2, momindex + 1) = - tau * p_vec(2) * p_vec(1) / mom3 ;

    } else if (2 == dim) {

        // NOTE THAT THESE ARE DIFFERENT IN 2d
        p.getH()(0, momindex)     = tau * (p_vec(1) * p_vec(1)) / mom3 ;
        p.getH()(1, momindex + 1) = tau * (p_vec(0) * p_vec(0)) / mom3 ;

        //offdiagonal momentum
        p.getH()(0, momindex + 1) = - tau * p_vec(0) * p_vec(1) / mom3 ;
        p.getH()(1, momindex + 0) = - tau * p_vec(1) * p_vec(0) / mom3 ;
    } else {
        std::cerr << ("Dimension of Geometric constraint is not 2 or 3. This will crash many things. You should feel bad.") << std::endl;
    }

    for (int row = 0; row < dim; ++row) {
        double posxmother = x_m(row);
        double posx       = x_vec(row);
        double momx       = p_vec(row);

        /** the direction of the momentum is very well known from the kinematic constraints
         *  that is why we do not extract the distance as a vector here
         * */
        p.getResiduals()(row) = posxmother + tau * momx / mom - posx ;
        p.getH()(row, posindexmother + row) = 1;
        p.getH()(row, posindex + row) = -1;
        p.getH()(row, tauindex) = momx / mom;
    }

    return ErrCode(ErrCode::Status::success);
}