#include "FitManager.h"

#include "ConstraintConfiguration.h"
#include "DecayChain.h"
#include "FitParams.h"

using namespace sct::ana;

FitManager::FitManager(ParticlePtr particle, const ConstraintConfiguration& config):
    m_particle(particle),
    m_config(config),
    m_useReferencing(false)
{
    m_decaychain = new DecayChain(particle, config);
    m_fitparams = new FitParams(m_decaychain->dim());
}

bool FitManager::fit() {
    // 1. initialize decay_chain
    // 2. fit while criterion isn't reached
    // 3. update tree if necessary
    m_decaychain->initialize(*m_fitparams);

    std::size_t nitermax = 5;
    for (auto m_niter = 0; m_niter < nitermax; ++m_niter) {
        if (0 == m_niter) {
          m_decaychain->filter(*m_fitparams);
        } else if (m_niter > 0 && m_useReferencing) {
          auto* tempState = new FitParams(*m_fitparams);
          m_decaychain->filterWithReference(*m_fitparams, *tempState);
          delete tempState;
        }
      }

    if (!(m_fitparams->testCovariance()))
        std::cerr << "bad covariance" << std::endl;

    return false;
}