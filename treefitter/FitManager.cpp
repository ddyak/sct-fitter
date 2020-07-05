#include "FitManager.h"

#include "ConstraintConfiguration.h"
#include "DecayChain.h"
#include "FitParams.h"

using namespace sct::ana;

FitManager::FitManager(ParticlePtr particle, const ConstraintConfiguration& config):
    m_particle(particle),
    m_config(config) 
{
    m_decaychain = new DecayChain(particle, config);
    m_fitparams = new FitParams(m_decaychain->dim());
}

bool FitManager::fit() {
    // 1. initialize decay_chain
    // 2. fit while criterion isn't reached
    // 3. update tree if necessary
    return false;
}