#include "DummyParticle.h"

#include "ParticleBase.h"
#include "FitParams.h"
#include "Projection.h"

using namespace sct::ana;


DummyParticle::DummyParticle(ParticlePtr particle, const ParticleBase* mother, const ConstraintConfiguration& config):
    ParticleBase(particle, mother, config)
{
    // todo: add flags from config
    // todo: read covariance from POD field
}

ErrCode DummyParticle::projectConstraint(const Constraint::Type type, 
                                const FitParams& fitparams,
                                Projection& p) const {
    ErrCode status;
    switch (type) {
    case Constraint::measurment:
        status |= projectMeasurmentConstraint(fitparams, p);
        break;
    default:
        throw std::runtime_error("Dummy: project Constraint unknown");
    }

    return status;
}

ErrCode DummyParticle::projectMeasurmentConstraint(const FitParams& fitparams, Projection& p) const {
    const int momindex = momIndex();

    p.getResiduals().segment(0, 3) = fitparams.getStateVector().segment(momindex, 3) - m_particle->momentum().transpose();

    for (int imom = 0; imom < 3; ++imom) {
        p.getH()(imom, momindex + imom) = 1;
    }

    p.getV().triangularView<Eigen::Lower>() = particle()->onlyMomentumErrorMatrix().triangularView<Eigen::Lower>();

    return ErrCode(ErrCode::Status::success);
}

void DummyParticle::addToConstraintList(std::vector<Constraint>& list, int depth) const {
    list.push_back(Constraint(this, Constraint::measurment, depth, 3, 3));
}

ErrCode DummyParticle::initParticle(FitParams& fitparams) {
    int momindex = momIndex();
    fitparams.getStateVector().segment(momindex, 3) = Eigen::Matrix<double, 3, 1>::Zero(3);

    fitparams.getStateVector().segment(momindex, 3) = m_particle->momentum();
    return ErrCode(ErrCode::Status::success);
}

ErrCode DummyParticle::initCovariance(FitParams& fitparams) const {
    ParticleBase::initCovariance(fitparams);
    return ErrCode(ErrCode::Status::success);
}