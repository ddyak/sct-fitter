#pragma once

#include "dataobjects/Particle.h"
#include "Constraint.h"
#include "ErrCode.h"

namespace sct::ana {

class ConstraintConfiguration;
class FitParams;

class ParticleBase {
public:
    /** default constructor  */
    ParticleBase(ParticlePtr particle, const ParticleBase* mother, const ConstraintConfiguration& config);

    /** create the according treeFitter particle obj for a basf2 particle type  */
    static ParticleBase* createParticle(ParticlePtr particle, const ParticleBase* mother, const ConstraintConfiguration& config);

    /**  get dimension of constraint */
    virtual int dim() const = 0;

    /** get index  */
    int index() const { return m_index; }

    /** this sets the index for momentum, position, etc. in the statevector  */
    virtual void updateIndex(int& offset);

    /** get vertex index (in statevector!) */
    virtual int posIndex() const { return -1; }

    /** get tau index */
    virtual int tauIndex() const { return -1; }

    /** get momentum index */
    virtual int momIndex() const { return -1; }

    // does the particle have a 3-momentum or a 4-momentum ?
    /** get momentum dimension */
    virtual bool hasEnergy() const = 0;

    /** get false  */
    virtual bool hasPosition() const = 0;

    /** get pdg mass  */
    double pdgMass() const { return m_pdgMass ; }

    /** add daughter  */
    ParticleBase* addDaughter(ParticlePtr particle, const ConstraintConfiguration& config);

    /** initialize particle params */
    virtual ErrCode initParticle(FitParams& fitparams) = 0;

    /** init covariance matrix */
    virtual ErrCode initCovariance(FitParams&) const;

    /** add to constraint list */
    virtual void addToConstraintList(std::vector<Constraint>& alist, int depth) const = 0;

    /** project constraint */
    virtual ErrCode projectConstraint(Constraint::Type, const FitParams&, Projection&) const = 0;

    /** project mass constraint using the particles parameters */
    ErrCode projectMassConstraintParticle(const FitParams&, Projection&) const;

protected:
    /** pointer to framework type  */
    ParticlePtr m_particle;

    /** motherparticle */
    const ParticleBase* m_mother;

    /** daughter container */
    std::vector<ParticleBase*> m_daughters;

    /** has all the constraint config */
    const ConstraintConfiguration& m_config;

private:
    /** index */
    int m_index;

    /** pdg mass  */
    const double m_pdgMass;
};

}