/**************************************************************************
 * BASF2 (Belle Analysis Framework 2)                                     *
 * Copyright(C) 2018 - Belle II Collaboration                             *
 *                                                                        *
 * Author: The Belle II Collaboration                                     *
 * Contributor: Wouter Hulsbergen, Francesco Tenchini, Jo-Frederik Krohn  *
 *                                                                        *
 * This software is provided "as is" without any warranty.                *
 **************************************************************************/
#pragma once

#include "ParticleBase.h"

namespace sct::ana {

  /** another unneccessary layer of abstraction */
  class InternalParticle : public ParticleBase {

  public:
    /** constructor */
    InternalParticle(ParticlePtr particle, const ParticleBase* mother, const ConstraintConfiguration& config);

    /** space reserved in fit params, if has mother then it has tau */
    virtual int dim() const override;

    /** position index in fit params*/
    virtual int posIndex() const override;

    /** tau index in fit params only if it has a mother */
    virtual int tauIndex() const override;

    /** momentum index in fit params depending on whether it has a mother  */
    virtual int momIndex() const override;

    /** has energy in fitparams  */
    virtual bool hasEnergy() const override { return true ; }

    /** has position index  */
    virtual bool hasPosition() const override { return !m_onlymomentum; }

    /** add to constraint list  */
    virtual void addToConstraintList(std::vector<Constraint>& list, int depth) const override;

    /** find out which constraint it is and project */
    virtual ErrCode projectConstraint(const Constraint::Type type, const FitParams& fitparams, Projection& p) const override;

private:
    /** project kinematical constraint */
    ErrCode projectKineConstraint(const FitParams&, Projection&) const;

    /** enforce conservation of momentum sum*/
  //  virtual void forceP4Sum(FitParams&) const override;

    /** set mass constraint flag */
    void setMassConstraint(bool b) { m_massconstraint = b ; }

    /** rotate in positive phi domain  */
    //double phidomain(const double);

  protected:

    /** init momentum of *this and daughters */
    ErrCode initMomentum(FitParams& fitparams) const;

    /** init particle in case it has a mother */
    //virtual ErrCode initParticleWithMother(FitParams& fitparams) override;

    /** init particle in case it has no mother */
    //virtual ErrCode initMotherlessParticle(FitParams& fitparams) override;

    virtual ErrCode initParticle(FitParams& fitparams) override;

    /** init covariance */
    virtual ErrCode initCovariance(FitParams&) const override;

  private:

    /** compare transverse track momentum*/
    // bool static compTrkTransverseMomentum(const RecoTrack* lhs, const RecoTrack* rhs);

    /** use dummy model without coordinate*/
    bool m_onlymomentum = true;

    /** has mass cosntraint */
    bool m_massconstraint = true;

    /** shares vertex with mother, that means decay vertex = productionvertex */
    bool m_shares_vertex_with_mother = false;

    /** use a geo metric constraint */
    bool m_geo_constraint = true;

    /** has lifetime constraint  */
    bool m_lifetimeconstraint = false;

    /** is conversion  */
    bool m_isconversion = false;

    /** automatically figure out if mother and particle vertex should be the same
     * and also add geometric constraints */
    bool m_automatic_vertex_constraining = false;
  } ;

}
