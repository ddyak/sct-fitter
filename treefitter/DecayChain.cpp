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

bool DecayChain::initialize(FitParams& par) {
    bool status = true;
    par.resetStateVector();
    status |= m_headOfChain->initParticle(par);
    par.resetCovariance();
    status |= m_headOfChain->initCovariance(par);
    initConstraintList();
    return status;
}

bool DecayChain::filter(FitParams& par) {
    bool status;
    par.resetCovariance();
    status |= m_headOfChain->initCovariance(par);
    par.resetChiSquare();
    for (auto constraint : m_constraintlist) {
      status |= constraint.filter(par);
    }
    return status;
}