#include "DecayChain.h"

#include "ParticleBase.h"
#include "FitParams.h"

using namespace sct::ana;


DecayChain::DecayChain(ParticlePtr particle, const ConstraintConfiguration& config):
    m_config(config),
    m_dim(0)
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
        std::cerr << constraint.name() << std::endl;
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