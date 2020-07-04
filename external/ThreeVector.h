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

#include <Eigen/Core>
#include <Eigen/Dense>

#include <iostream>
#include <cmath>

namespace sct::kine {

template <typename ftype = double>
class ThreeVector : public Eigen::Matrix<ftype, 1, 3> {
 public:
    ThreeVector<ftype>(void):Eigen::Matrix<ftype, 1, 3>() {}
    
    // This constructor allows you to construct MyVectorType from Eigen expressions
    template<typename OtherDerived>
    ThreeVector<ftype>(const Eigen::MatrixBase<OtherDerived>& other)
        : Eigen::Matrix<ftype, 1, 3>(other) {}

    // This method allows you to assign Eigen expressions to ThreeVector (without this it doesn't word)
    template<typename OtherDerived>
    ThreeVector<ftype>& operator=(const Eigen::MatrixBase <OtherDerived>& other) {
        this->Eigen::Matrix<ftype, 1, 3>::operator=(other);
        return *this;
    }

    Eigen::Matrix<ftype, 1, 3>& asEigen() const {
        return *this;
    }

    using Eigen::Matrix<ftype, 1, 3>::x;
    using Eigen::Matrix<ftype, 1, 3>::y;
    using Eigen::Matrix<ftype, 1, 3>::z;
    using Eigen::Matrix<ftype, 1, 3>::norm;
    using Eigen::Matrix<ftype, 1, 3>::squaredNorm;
    using Eigen::Matrix<ftype, 1, 3>::head;

    ThreeVector<ftype>(Eigen::Matrix<ftype, 1, 3> a) : Eigen::Matrix<ftype, 1, 3>(a) {}
    ThreeVector<ftype>(ftype x, ftype y, ftype z) : Eigen::Matrix<ftype, 1, 3>{x,y,z} {}

    static ThreeVector<ftype> fromRPhiCosth(double r, double phi, double costh) {
        const double sphi = std::sin(phi);
        const double cphi = std::cos(phi);
        const double sinth = std::sqrt(1 - costh*costh);
        return {r*sinth*cphi, r*sinth*sphi, r*costh};
    }

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

    inline ftype costh() const {
        return z() / norm();
    }
    /** sine of polar angle for a three-vector */
    inline ftype sinth() const {
        return sqrt(1. - z()*z() / squaredNorm());
    }
    /** azimutal angle in the range [-pi, pi] */
    inline ftype phi() const {
        const ftype sth = sinth();
        const ftype cphi = x() / sth;
        return y() / sth > 0 ? acos(cphi) : -acos(cphi);
    }
    /** cosine of the azimutal angle */
    inline ftype cosphi() const {
        return x() / sinth();
    }
    /** sine of the azimutal angle*/
    inline ftype sinphi() const {
        return y() / sinth();
    }
    /** length of the transverse component */
    inline ftype pt() const {
        return head(2).norm();
    }
    
    inline ftype norm2() const {
        return squaredNorm();
    }

    ThreeVector<ftype>& RotateX(double angle) {
        return *this = Eigen::AngleAxis(angle, Eigen::Vector3f::UnitX()) * *this;
    }

    ThreeVector<ftype>& RotateY(double angle) {
        return *this = Eigen::AngleAxis(angle, Eigen::Vector3f::UnitY()) * *this;
    }

    ThreeVector<ftype>& RotateZ(double angle) {
        return *this = Eigen::AngleAxis(angle, Eigen::Vector3f::UnitZ()) * *this;
    }
};

template <typename ftype> 
inline ftype angle(const ThreeVector<ftype>& v1, const ThreeVector<ftype> v2) {
    ftype ptot2 = norm2(v1) * norm2(v2);
    if (fabs(ptot2) < 1e-15) {
        std::cerr << "SctKine::angle for zero vector" << std::endl;
        return 0.0;
    }
    return std::acos(v1.dot(v2) / sqrt(ptot2));
}

}  // namespace sct::kine
