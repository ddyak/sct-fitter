#pragma once

#include "ParticleBase.h"

namespace sct::ana {

class FitParams;
class ConstraintConfiguration;
class Projection;

class DummyParticle : public ParticleBase {
public:
    DummyParticle(ParticlePtr particle, const ParticleBase* mother, const ConstraintConfiguration& config);

    /**  get dimension of constraint */
    virtual int dim() const override { return 4; }

    /** get momentum index */
    virtual int momIndex() const { return index(); }

    // does the particle have a 3-momentum or a 4-momentum ?
    /** get momentum dimension */
    virtual bool hasEnergy() const { return true; }

    /** get false  */
    virtual bool hasPosition() const { return false; }

    /** project constraint */
    virtual ErrCode projectConstraint(const Constraint::Type type, 
                                    const FitParams& fitparams,
                                    Projection& p) const override;
private:
    /** project kinematical constraint */
    ErrCode projectKineConstraint(const FitParams& fitparams, Projection& p) const;

    /** project measurment constraint */
    ErrCode projectMeasurmentConstraint(const FitParams& fitparams, Projection& p) const;

public:
    /** add to constraint list */
    virtual void addToConstraintList(std::vector<Constraint>& list, int depth) const;

protected:
    /** init momentum of *this and daughters */
    bool initMomentum(FitParams& fitparams) const;

    virtual ErrCode initParticle(FitParams& fitparams) override;

    ErrCode initCovariance(FitParams& fitparams) const override;

private:
    bool m_massconstraint = true;

    /** only lower triangle filled! */
    Eigen::Matrix<double, 4, 4> m_covariance;

    /** only lower triangle filled! */
    Eigen::Matrix<double, 1, 4> m_momentum;
};


}