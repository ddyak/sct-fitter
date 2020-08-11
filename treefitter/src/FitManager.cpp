#include "FitManager.h"

#include "ConstraintConfiguration.h"
#include "DecayChain.h"
#include "FitParams.h"
#include "ParticleBase.h"

using namespace sct::ana;
using namespace sct;


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
    for (std::size_t m_niter = 0; m_niter < nitermax; ++m_niter) {
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

    updateTree(m_particle, true);
    m_chiSquare = m_fitparams->getChiSquare();

    return false;
}

void FitManager::updateTree(ParticlePtr cand, const bool isTreeHead) const {
    const bool updateableMother = updateCand(cand, isTreeHead);
    
    if (updateableMother) {
        const int ndaughters = cand->nDaughters();
        for (int i = 0; i < ndaughters; i++) {
            ParticlePtr daughter = cand->daughter(i);
            updateTree(daughter, false);
        }
    }
}

bool FitManager::updateCand(ParticlePtr cand, const bool isTreeHead) const {
    const ParticleBase* pb = m_decaychain->locate(cand);
    if (pb) {
        updateCand(*pb, cand, isTreeHead);
    } else {
        std::cerr << "Can't find candidate" << std::endl;
        //B2ERROR("Can't find candidate " << cand.getName() << "in tree " << m_particle->getName());
    }
    return pb != nullptr;
}

void FitManager::updateCand(const ParticleBase& pb,
                            ParticlePtr cand, const bool isTreeHead) const
{
    int posindex = pb.posIndex();
    if (posindex < 0 && pb.mother()) {
        posindex = pb.mother()->posIndex();
    }
    if (m_updateDaugthers || isTreeHead) {
        // if (posindex >= 0) {
        //     const TVector3 pos(m_fitparams->getStateVector()(posindex),
        //                         m_fitparams->getStateVector()(posindex + 1),
        //                         m_fitparams->getStateVector()(posindex + 2));
        //     cand.setVertex(pos);
        //     if (&pb == m_decaychain->cand()) { // if head
        //         const double fitparchi2 = m_fitparams->chiSquare();
        //         cand.setPValue(TMath::Prob(fitparchi2, m_ndf));//if m_ndf<1, this is 0.
        //         setExtraInfo(&cand, "chiSquared", fitparchi2);
        //         setExtraInfo(&cand, "modifiedPValue", TMath::Prob(fitparchi2, 3));
        //         setExtraInfo(&cand, "ndf", m_ndf);
        //     }
        // }

        const int momindex = pb.momIndex();
        if (pb.hasEnergy()) {
            cand->set4Momentum(m_fitparams->getStateVector().segment(momindex, 4));
        } else {
            cand->set3Momentum(m_fitparams->getStateVector().segment(momindex, 3));
        }

        // TMatrixFSym cov7b2(7);
        // getCovFromPB(&pb, cov7b2);
        // cand.setMomentumVertexErrorMatrix(cov7b2);
    }

    // if (pb.tauIndex() > 0) {
    //     std::tuple<double, double>tau  = getDecayLength(cand);
    //     std::tuple<double, double>life = getLifeTime(cand);
    //     setExtraInfo(&cand, std::string("decayLength"), std::get<0>(tau));
    //     setExtraInfo(&cand, std::string("decayLengthErr"), std::get<1>(tau));
    //     setExtraInfo(&cand, std::string("lifeTime"), std::get<0>(life));
    //     setExtraInfo(&cand, std::string("lifeTimeErr"), std::get<1>(life));
    // }
}

FourVector FitManager::getMomentum(ParticlePtr particle) const {
    FourVector momentum;
    auto pb = m_decaychain->locate(particle);
    if (pb->hasEnergy()) {
        auto momIdx = pb->momIndex();
        momentum = m_fitparams->getStateVector().segment(momIdx, 4);
    } else {
        double mass = pb->pdgMass();
        auto momIdx = pb->momIndex();
        ThreeVector mom3 = m_fitparams->getStateVector().segment(momIdx, 3);
        momentum = FourVector::fromPMass(mom3, mass);
    }
    return momentum;
}
