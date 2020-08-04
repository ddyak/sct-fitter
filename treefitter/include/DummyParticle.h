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
    virtual int dim() const override { return 3; }

    /** get momentum index */
    virtual int momIndex() const { return index(); }

    // does the particle have a 3-momentum or a 4-momentum ?
    /** get momentum dimension */
    virtual bool hasEnergy() const { return false; }

    /** get false  */
    virtual bool hasPosition() const { return false; }

    /** add to constraint list */
    virtual void addToConstraintList(std::vector<Constraint>& list, int depth) const;

    /** project constraint */
    virtual ErrCode projectConstraint(const Constraint::Type type, 
                                    const FitParams& fitparams,
                                    Projection& p) const override;
private:
    /** project measurment constraint */
    ErrCode projectMeasurmentConstraint(const FitParams& fitparams, Projection& p) const;

protected:
    virtual ErrCode initParticle(FitParams& fitparams) override;

    ErrCode initCovariance(FitParams& fitparams) const override;

private:
    bool m_massconstraint = true;

    /** only lower triangle filled! */
    Eigen::Matrix<double, 3, 3> m_covariance;

    /** only lower triangle filled! */
    Eigen::Matrix<double, 1, 3> m_momentum;
};


}