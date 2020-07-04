/**************************************************************************
 * Aurora Framework (SCT software)                                        *
 * Copyright(C) 2019 - SCT Collaboration                                  *
 *                                                                        *
 * Author: The SCT Collaboration                                          *
 * Contributors: Daniil Yakovlev, Vitaly Vorobyev                         *
 *                                                                        *
 * This software is provided "as is" without any warranty.                *
 **************************************************************************/

#pragma once

#include "ThreeVector.h"

#include <Eigen/Core>
#include <Eigen/Dense>

#include <iostream>
#include <cmath>

namespace sct::kine {

template <typename ftype = double>  // p, E;
class FourVector : public Eigen::Matrix<ftype, 1, 4> {
 public:
    static FourVector<ftype> fromPMass(const ThreeVector<ftype>& p, ftype mass) {
        return {p.x(), p.y(), p.z(), std::sqrt(mass*mass + p.squaredNorm())};
    }

    FourVector<ftype>(void) : Eigen::Matrix<ftype, 1, 4>() {}

    // This constructor allows you to construct MyVectorType from Eigen expressions
    template<typename OtherDerived>
    FourVector<ftype>(const Eigen::MatrixBase<OtherDerived>& other)
        : Eigen::Matrix<ftype, 1, 4>(other) {}

    // This method allows you to assign Eigen expressions to ThreeVector (without this it doesn't word)
    template<typename OtherDerived>
    FourVector<ftype>& operator=(const Eigen::MatrixBase <OtherDerived>& other) {
        this->Eigen::Matrix<ftype, 1, 4>::operator=(other);
        return *this;
    }

    using Eigen::Matrix<ftype, 1, 4>::x;
    using Eigen::Matrix<ftype, 1, 4>::y;
    using Eigen::Matrix<ftype, 1, 4>::z;
    using Eigen::Matrix<ftype, 1, 4>::norm;
    using Eigen::Matrix<ftype, 1, 4>::squaredNorm;
    using Eigen::Matrix<ftype, 1, 4>::head;

    FourVector<ftype>(ftype x, ftype y, ftype z, ftype e) :
        Eigen::Matrix<ftype, 1, 4>{x, y, z, e} {}

    FourVector<ftype>(const ThreeVector<ftype>& p, ftype e) :
        Eigen::Matrix<ftype, 1, 4> {p.x(), p.y(), p.z(), e} {}   

    inline ftype& px() {
        return x();
    }

    inline ftype px() const {
        return x();
    }

    inline ftype& py() {
        return y();
    }

    inline ftype py() const {
        return y();
    }

    inline ftype& pz() {
        return z();
    }

    inline ftype pz() const {
        return z();
    }

    inline ftype& E() {
        return (*this)[3];
    }

    inline ftype E() const {
        return (*this)[3];
    }

    inline ftype& t() {
        return (*this)[3];
    }

    inline ftype t() const {
        return (*this)[3];
    }

    ThreeVector<ftype> p() const {
        return head(3);
    }

    inline ThreeVector<ftype> r() const {
        return p();
    }

    void setP(ThreeVector<ftype> p) {
        head(3) = p.head(3);
    }

    ThreeVector<ftype> BoostVec() const {
        return -p() / E();
    }

    Eigen::Matrix<ftype, 4, 4> BoostMatrix() const {
        //  The Matrix for a Lorentz boosts is:
        //  ~~~
        //   | 1+gamma'*bx*bx  gamma'*bx*by   gamma'*bx*bz  gamma*bx |
        //   |  gamma'*by*bx  1+gamma'*by*by  gamma'*by*bz  gamma*by |
        //   |  gamma'*bz*bx   gamma'*bz*by  1+gamma'*bz*bz gamma*bz |
        //   |    gamma*bx       gamma*by       gamma*bz     gamma   |
        //  ~~~
        //  with the boost vector b=(bx,by,bz) and gamma=1/Sqrt(1-beta*beta)
        //  and gamma'=(gamma-1)/beta*beta.
        ThreeVector<ftype> bvec = BoostVec();
        const double betaSq = bvec.squaredNorm();
        const double gamma = 1. / std::sqrt(1. - betaSq);
        const double gampr = (gamma - 1.) / betaSq;
        Eigen::Matrix<ftype, 4, 4> bmtx;
        bmtx.topLeftCorner(3, 3) = Eigen::Matrix<ftype, 3, 3>::Identity() + gampr * bvec.transpose() * bvec;
        bmtx.topRightCorner(1, 3) = gamma * bvec;
        bmtx.bottomLeftCorner(3, 1) = gamma * bvec.transpose();
        bmtx(3, 3) = gamma;
        return bmtx;
    }

    double squaredMass() const {
        return E()*E() - p().norm2();
    }

    inline double mass() const {
        double msq = squaredMass();
        if (msq < 0) {
            std::cerr << "negative mass! " << msq << std::endl;
            return 0.;
        }
        return sqrt(msq);
    }

    FourVector<ftype> Boost(const Eigen::Matrix<ftype, 4, 4> bmtx) {
        return *this = bmtx * *this;
    }

    FourVector<ftype> Boosted(const Eigen::Matrix<ftype, 4, 4> bmtx) const {
        return bmtx * *this;
    }

    FourVector<ftype> Boosted(ThreeVector<ftype> bv) const {
        const ftype bbeta = bv.norm();
        if (fabs(bbeta) < 1e-15) {
            return *this;
        } else if (abs(bbeta) > 1.) {
            std::cerr << "boost with |beta| > 1: " << bbeta << std::endl;
            return *this;
        }
        const ftype bgamma = 1./sqrt(1. - bbeta*bbeta); 
        ThreeVector<ftype> bdirect = -bv / bbeta; 
        const ftype dotnr = p() * bdirect.transpose();
        return {p() + ((bgamma - 1.)*dotnr - bgamma*E()*bbeta)*bdirect, bgamma * (E() - bbeta * dotnr)};
    }

    FourVector<ftype>& Boost(ThreeVector<ftype> bv) {
        *this = Boosted(std::move(bv));
        return *this;
    }

    FourVector<ftype>& RotateX(double angle) {
        head(3) *= Eigen::AngleAxis(angle, ThreeVector<ftype>(1, 0, 0)).toRotationMatrix();
        return *this;
    }

    FourVector<ftype>& RotateY(double angle) {
        head(3) *= Eigen::AngleAxis(angle, ThreeVector<ftype>(0, 1, 0)).toRotationMatrix();
        return *this;
    }

    FourVector<ftype>& RotateZ(double angle) {
        head(3) *= Eigen::AngleAxis(angle, ThreeVector<ftype>(0, 0, 1)).toRotationMatrix();
        return *this;
    }
};

}  // namespace sct::kine
