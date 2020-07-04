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

#include "Particle.h"
#include "external/EvtPdlDatabase.h"

#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <algorithm>  // remove, equal
#include <random>
#include <utility>

using namespace std;
using namespace sct::ana;

using sct::ThreeVector;
using sct::FourVector;
using sct::ErrMatrix;
using sct::MomentumErrMatrix;
using sct::PositionErrMatrix;
using sct::ftype;
using sct::comm::EvtPdlDatabase;

template<typename T>
bool setOverlap(const unordered_set<T>& lhs, const unordered_set<T>& rhs) {
    for (const auto& item : lhs) {
        if (rhs.find(item) != rhs.end()) return true;
    }
    return false;
}

Particle::Particle() :
    m_pdgCode(0), m_p(0,0,0,0), m_r(0,0,0),
    m_flavorType(c_Unflavored), m_particleType(c_Undefined) {
}

Particle::Particle(const ThreeVector& momentum, ftype mass, int pdgCode, EParticleType type) :
    m_pdgCode(pdgCode), m_p(FourVector::fromPMass(momentum, mass)),
    m_particleType(type), m_mass(mass) {
    setFlavorType();
}

Particle::Particle(FourVector momentum, int pdgCode, EParticleType type) :
    m_pdgCode(pdgCode), m_p(move(momentum)), m_particleType(type) {
    setFlavorType();
}

Particle::Particle(int pdgCode) : m_pdgCode(pdgCode),
    m_p(0,0,0,0), m_r(0,0,0), m_particleType(c_Composite) {
    setFlavorType();
}

Particle::Particle(const vector<ParticlePtr>& daughters, int pdgCode) :
    Particle(pdgCode) {
    for (auto daug : daughters) {
        appendDaughter(daug);
    }
}

int Particle::mdstIndex() const {
    if (m_mdstAssociated.size() != 1) return 0;
    return *m_mdstAssociated.cbegin();
}

const ThreeVector& Particle::vertex() const {
    return m_r;
}

void Particle::set3Momentum(const ThreeVector& p3) {
    m_p.head(3) = p3;
    m_mass = nullopt;
}

void Particle::set4Momentum(const ThreeVector& p3, ftype E) {
    m_p << p3, E;
    m_mass = nullopt;
}

void Particle::set4Momentum(FourVector p4) {
    m_p = move(p4);
    m_mass = nullopt;
}

/** Sets position (decay vertex) */
void Particle::setVertex(ThreeVector vertex) {
    m_r = move(vertex);
}

/** Sets 7x7 error matrix (order: px,py,pz,E,x,y,z) */
void Particle::setMomentumVertexErrorMatrix(ErrMatrix errMatrix) {
    m_errMatrix = move(errMatrix);
}

const ErrMatrix& Particle::momentumVertexErrorMatrix() const {
    return m_errMatrix;
}

MomentumErrMatrix Particle::momentumErrorMatrix() const {
    return m_errMatrix.topLeftCorner<c_DimMomentum, c_DimMomentum>();
}

PositionErrMatrix Particle::vertexErrorMatrix() const {
    return m_errMatrix.bottomRightCorner<c_DimPosition, c_DimPosition>();
}

void checkPdgId(int pdgCode) {
    if (!EvtPdlDatabase::Instance().isKnownId(pdgCode))
        throw runtime_error("PDG=" + to_string(pdgCode) + " is not in the EvtPdlDatabase");
}

void Particle::updateMass(int pdgCode) {
    checkPdgId(pdgCode);
    m_mass = EvtPdlDatabase::Instance().mass(pdgCode);
    m_p = FourVector::fromPMass(m_p.p(), m_mass.value());
}

ftype Particle::pdgMass(void) const {
    checkPdgId(m_pdgCode);
    return EvtPdlDatabase::Instance().mass(m_pdgCode);
}

int Particle::charge(void) const {
    checkPdgId(m_pdgCode);
    return EvtPdlDatabase::Instance().charge(m_pdgCode);
}

ParticlePtr Particle::daughter(size_t i) const {
    if (i < m_daughters.size()) return m_daughters[i];
    return nullptr;
}

vector<ParticlePtr> Particle::finalStateDaughters() const {
    vector<ParticlePtr> fspDaughters;
    fillFSPDaughters(fspDaughters);
    return fspDaughters;
}

bool Particle::appendDaughter(ParticlePtr daughter) {
    for (size_t idx : daughter->m_mdstAssociated) {
        if (m_mdstAssociated.find(idx) != m_mdstAssociated.end()) {
            return false;
        }
        m_mdstAssociated.insert(idx);
    }
    m_p += daughter->m_p;
    m_mass = nullopt;
    m_daughters.push_back(daughter);

    m_particleType = c_Composite;
    return true;
}

ftype Particle::energy() const {
    return m_p.E();
}

const FourVector& Particle::fourMomentum() const {
    return m_p;
}

ThreeVector Particle::momentum() const {
    return m_p.p();
}

ftype Particle::momentumMagnitude() const   {
    return m_p.norm();
}

ftype Particle::p() const {
    return m_p.norm();
}

ftype Particle::px() const {
    return m_p.px();
}

ftype Particle::py() const {
    return m_p.py();
}

ftype Particle::pz() const {
    return m_p.pz();
}

ftype Particle::x() const {
    return m_r.x();
}

ftype Particle::y() const {
    return m_r.y();
}

ftype Particle::z() const {
    return m_r.z();
}

bool Particle::overlapsWith(const Particle &oParticle) const {
    if (m_mdstAssociated.size() < oParticle.m_mdstAssociated.size()) {
        return setOverlap(m_mdstAssociated, oParticle.m_mdstAssociated);
    }
    return setOverlap(oParticle.m_mdstAssociated, m_mdstAssociated);
}

void Particle::fillFSPDaughters(vector<ParticlePtr>& fspDaughters) const {
    if (!nDaughters()) {  // this is FSP
        fspDaughters.push_back(const_cast<Particle*>(this)->shared_from_this());
        return;
    }
    // this is not FSP (go one level down)
    for (const auto d : daughters()) {
        d->fillFSPDaughters(fspDaughters);
    }
}

void Particle::setFlavorType() {
    m_flavorType = c_Flavored;
    if (m_pdgCode < 0) return;
    if (m_pdgCode == 22 || m_pdgCode == 310 || m_pdgCode == 130) {
        // gamma, K_S or K_L
        m_flavorType = c_Unflavored;
        return;
    }
    int nnn = m_pdgCode / 10;
    int q3 = nnn % 10; nnn /= 10;
    int q2 = nnn % 10; nnn /= 10;
    int q1 = nnn % 10;
    if (q1 == 0 && q2 == q3) m_flavorType = c_Unflavored; // unflavored meson
}

size_t Particle::nDaughters(void) const {
    return m_daughters.size();
}

ftype Particle::mass() const {
    if (m_mass == nullopt) {
        m_mass = m_p.mass();
    }
    return m_mass.value();
}

int Particle::pdgCode(void) const {
    return m_pdgCode;
}

Particle::EFlavorType Particle::flavorType() const {
    return m_flavorType;
}

Particle::EParticleType Particle::particleType() const {
    return m_particleType;
}

const std::vector<ParticlePtr>& Particle::daughters() const {
    return m_daughters;
}

const string& Particle::name() const {
    return EvtPdlDatabase::Instance().name(m_pdgCode);
}

bool Particle::hasExtraInfo(const string& name) const {
    return m_extraInfo.find(name) != m_extraInfo.cend();
}

void Particle::clearExtraInfo() {
    m_extraInfo.clear();
}

optional<ftype> Particle::extraInfo(const string& name) const {
    if (auto it = m_extraInfo.find(name); it != m_extraInfo.cend()) {
        return make_optional(it->second);
    }
    return nullopt;
}

void Particle::setExtraInfo(const string& name, ftype value) {
    m_extraInfo[name] = value;
}

Particle::Particle(const Particle& p) :
    enable_shared_from_this<Particle>(),
    m_pdgCode(p.m_pdgCode),
    m_p(p.m_p),
    m_r(p.m_r),
    m_errMatrix(p.m_errMatrix),
    m_flavorType(p.m_flavorType),
    m_particleType(p.m_particleType),
    m_mdstAssociated(p.m_mdstAssociated)
    // m_trkfit(p.m_trkfit),
    // m_klmclu(p.m_klmclu),
    // m_caloclu(p.m_caloclu),
    // m_calohyp(p.m_calohyp) 
    {}

ParticlePtr Particle::deepCopy() const {
    if (!nDaughters()) {  // FSP
        return ParticlePtr(new Particle(*this));
    }
    ParticlePtr p = make_shared<Particle>();
    for (size_t idx = 0; idx < nDaughters(); ++idx) {
        p->appendDaughter(daughter(idx)->deepCopy());
    }
    return p;
}
