#include "DecayChain.h"

#include "ParticleBase.h"
#include "FitParams.h"

using namespace sct::ana;


DecayChain::DecayChain(ParticlePtr particle, const ConstraintConfiguration& config):
    m_dim(0),
    m_config(config)
{
    m_headOfChain = ParticleBase::createParticle(particle, nullptr, config);
    m_headOfChain->updateIndex(m_dim);
}

void DecayChain::initConstraintList() {
    m_constraintlist.clear();
    m_headOfChain->addToConstraintList(m_constraintlist, 0);
    // removeConstraintFromList(); ddyak: remove from config, useless right now
    std::sort(m_constraintlist.begin(), m_constraintlist.end());
}

ErrCode DecayChain::initialize(FitParams& par) {
    ErrCode status;
    par.resetStateVector();
    status |= m_headOfChain->initParticle(par);
    par.resetCovariance();
    status |= m_headOfChain->initCovariance(par);
    initConstraintList();
    return status;
}

ErrCode DecayChain::filter(FitParams& par) {
    ErrCode status;
    par.resetCovariance();
    status |= m_headOfChain->initCovariance(par);
    par.resetChiSquare();
    for (auto constraint : m_constraintlist) {
        status |= constraint.filter(par);
    }
    return status;
}

ErrCode DecayChain::filterWithReference(FitParams& par, const FitParams& ref) {
    ErrCode status;
    par.resetCovariance();
    status |= m_headOfChain->initCovariance(par);
    par.resetChiSquare();
    for (auto constraint : m_constraintlist) {
        status |= constraint.filterWithReference(par, ref);
    }
    return status;
}

const ParticleBase* DecayChain::locate(ParticlePtr particle) const {
    const ParticleBase* rc(nullptr);
    const auto mapRow = m_particleMap.find(particle);

    if (mapRow == m_particleMap.end()) {
        //  take head of chain and recursively find particle in it
        rc = m_headOfChain->locate(particle);

        if (rc && rc->particle()) {
            const_cast<DecayChain*>(this)->m_particleMap[rc->particle()] = rc;
        }
    } else {
        //only used for "head of tree"
        rc = mapRow->second;// (B2::Particle, Particlebase)
    }
    return rc;
}