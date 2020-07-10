#include "ParticleBase.h"


#include "ConstraintConfiguration.h"
#include "InternalParticle.h"
#include "DummyParticle.h"
#include "Projection.h"


using namespace sct::ana;


ParticleBase::ParticleBase(ParticlePtr particle, const ParticleBase* mother, const ConstraintConfiguration& config):
    m_particle(particle),
    m_mother(mother),
    m_config(config),
    m_index(0),
    m_pdgMass(particle->mass())
{}

/* static */ ParticleBase* ParticleBase::createParticle(ParticlePtr particle, const ParticleBase* mother, const ConstraintConfiguration& config)
{
    return new DummyParticle(particle, mother, config);
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

bool ParticleBase::initCovariance(FitParams& fitparams) const
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

    return true;
}

bool ParticleBase::projectMassConstraintParticle(const FitParams& fitparams, Projection& p) const {
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

    return true;
}