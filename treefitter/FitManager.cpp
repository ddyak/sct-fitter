#include "FitManager.h"
#include "ConstraintConfiguration.h"

using namespace sct::ana;

FitManager::FitManager(ParticlePtr particle, const ConstraintConfiguration& config):
    m_config(config) {
    // m_decaychain = new DecayChain(particle, config);
    // m_fitparams = new FitParams(m_decaychain->dim());
}