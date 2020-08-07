/**************************************************************************
 * Aurora (SCT software framework) 2020                                   *
 *                                                                        *
 * Author: The SCT Collaboration                                          *
 * Contributors: Vitaly Vorobyev                                         *
 *                                                                        *
 * Based on:                                                              *
 * BASF2 (Belle Analysis Framework 2)                                     *
 * Copyright(C) 2010 - Belle II Collaboration                             *
 *                                                                        *
 * Author: The Belle II Collaboration                                     *
 * Contributors: Anze Zupanc, Marko Staric                                *
 *                                                                        *
 * This software is provided "as is" without any warranty.                *
 **************************************************************************/

#pragma once

#include <vector>
#include <memory>
#include <optional>
#include <unordered_set>
#include <unordered_map>

// #include "RecoDataObjects/Track.h"
// #include "RecoDataObjects/CaloCluster.h"

#include "external/KineTypedefs.h"

namespace sct::ana {

class Particle;
using ParticlePtr = std::shared_ptr<Particle>;
using ParticleConstPtr = std::shared_ptr<const Particle>;

/** Class to store reconstructed particles */
class Particle : public std::enable_shared_from_this<Particle> {
 public:
    enum EParticleType {
        c_Undefined, c_Track, c_CaloCluster,
        c_KLMCluster, c_MCParticle, c_Composite
    };
    enum EFlavorType {c_Unflavored = 0, c_Flavored = 1};
    using ExtraInfoMap = std::unordered_map<std::string, ftype>;

 private:
    int m_pdgCode;
    FourVector m_p;  // [px, py, pz, E]
    ThreeVector m_r;
    ErrMatrix m_errMatrix;  // [px, py, pz, E, x, y, z]
    std::vector<ParticlePtr> m_daughters;
    EFlavorType m_flavorType;
    EParticleType m_particleType;
    std::unordered_set<size_t> m_mdstAssociated;
    ExtraInfoMap m_extraInfo;
    mutable std::optional<ftype> m_mass;

    // reco::TrackFitResultPtr m_trkfit;
    // reco::KLMClusterPtr m_klmclu;
    // reco::CaloClusterPtr m_caloclu;
    // reco::CaloCluster::EHypothesisBit m_calohyp;
    ftype m_pvalue;

    // Copy ctor only for final-state-particles
    Particle(const Particle&);

 public:
    /** All private members are set to 0. Particle type is set to c_Undefined. */
    Particle();
    /** All private members are set to 0. */
    Particle(int pdgCode);
    /** All other private members are set to their default values (0). */
    Particle(const ThreeVector& momentum, ftype mass, int pdgCode, EParticleType particleType=c_Undefined);
    /** All other private members are set to their default values (0). */
    Particle(FourVector momentum, int pdgCode, EParticleType particleType=c_Undefined);

    /** Constructor for cartesian particles. */
    Particle(int pdgCode, const ThreeVector& momentum, const Eigen::Matrix<double, 3, 3>& momCovariance);
    /** Constructor for composite particles. */
    Particle(const std::vector<ParticlePtr>& daughters, int pdgCode);
    
    /** Constructor from MC Particle of podio tree */
    // Particle(const sct::MCParticle &mcparticle);
    /** Constructor from MC Particle of podio tree  */
    // Particle(const sct::Particle &particle, int pdg, size_t partIndex);

    // reco::TrackFitResultPtr track() const {return m_trkfit;}
    // reco::KLMClusterPtr klmCluster() const {return m_klmclu;}
    // reco::CaloClusterPtr caloCluster() const {return m_caloclu;}
    // reco::CaloCluster::EHypothesisBit caloClusterEHypothesisBit() const {return m_calohyp;}

    // The only way to get a copy of the decay tree
    ParticlePtr deepCopy() const;

    void set3Momentum(const ThreeVector& p3);
    void set4Momentum(const ThreeVector& p3, ftype mass);
    void set4Momentum(FourVector p4);
    void setVertex(ThreeVector vertex);
    /** Sets 7x7 error matrix (order: px,py,pz,E,x,y,z) */
    void setMomentumVertexErrorMatrix(ErrMatrix errMatrix);
    /** Sets chi^2 probability of fit */
    void setPValue(ftype pValue) { m_pvalue = pValue;}
    ftype pValue() const {return m_pvalue;}
    void updateMass(int pdgCode);
    bool appendDaughter(ParticlePtr daughter);

    int pdgCode(void) const;
    int charge(void) const;
    EFlavorType flavorType() const;
    EParticleType particleType() const;
    ftype mass() const;
    /** Particle nominal mass */
    ftype pdgMass(void) const;
    ftype energy() const;

    /** Returns four-momentum vector [px, py, pz, E] */
    const FourVector& fourMomentum() const;
    ThreeVector momentum() const;
    ftype momentumMagnitude() const;
    /** Returns momentum magnitude (same as momentumMagnitude but with shorter name) */
    ftype p() const;
    ftype px() const;
    ftype py() const;
    ftype pz() const;

    /** Returns vertex position (POCA for charged, IP for neutral FS particles) */
    const ThreeVector& vertex() const;
    ftype x() const;
    ftype y() const;
    ftype z() const;
    int mdstIndex() const;

    /** 7x7 error matrix (order: px,py,pz,E,x,y,z) */
    const ErrMatrix& momentumVertexErrorMatrix() const;
    /** 4x4 momentum error matrix (order: px,py,pz,E) */
    MomentumErrMatrix momentumErrorMatrix() const;
    /** 3x3 momentum error matrix (order: px,py,pz,E) */
    Eigen::Matrix<double, 3, 3> onlyMomentumErrorMatrix() const;
    /** 3x3 position error matrix (order: x,y,z) */
    PositionErrMatrix vertexErrorMatrix() const;
    size_t nDaughters(void) const;
    ParticlePtr daughter(size_t i) const;
    const std::vector<ParticlePtr>& daughters() const;

    /** Returns a vector of pointers to Final State daughter particles
      * @return vector of pointers to final state daughter particles */
    std::vector<ParticlePtr> finalStateDaughters() const;

    /** Returns true if final state ancessors of oParticle overlap
      * @param oParticle reference to particle
      * @return true if overlap, otherwise false */
    bool overlapsWith(const Particle &oParticle) const;

    /** Return name of this particle. */
    const std::string& name() const;

    /** Remove all stored extra info fields */
    void clearExtraInfo();

    /** Return given value if set.
      * throws std::runtime_error if variable is not set. */
    std::optional<ftype> extraInfo(const std::string &name) const;
    bool hasExtraInfo(const std::string &name) const;
    /** Sets the user-defined data of given name to the given value */
    void setExtraInfo(const std::string &name, ftype value);

 private:
    /** Fill final state particle daughters into a vector
      * Function is called recursively
      * @param fspDaughters vector of daughter particles */
    void fillFSPDaughters(std::vector<ParticlePtr> &fspDaughters) const;

    /** sets m_flavorType using m_pdgCode */
    void setFlavorType();

//    void setParticleType(const sct::Particle& particle);
};

}  // namespace sct::ana
