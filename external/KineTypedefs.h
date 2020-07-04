/**************************************************************************
 * Aurora Framework (SCT software)                                        *
 * Copyright(C) 2019 - SCT Collaboration                                  *
 *                                                                        *
 * Author: The SCT Collaboration                                          *
 * Contributors: Vitaly Vorobyev                                          *
 *                                                                        *
 * This software is provided "as is" without any warranty.                *
 **************************************************************************/

#pragma once

#include <Eigen/Core>
#include <Eigen/Dense>

#include "ThreeVector.h"
#include "FourVector.h"

namespace sct {
/** error matrix dimensions and size of 1D representation */
enum {
    c_DimPosition = 3,
    c_DimMomentum = 4,
    c_DimHelJacob = 5,
    c_DimHelUncer = 5,  // Helix uncertainty matrix
    c_DimHelUncerCart = 6,
    c_DimMatrix = 7
};

/** Float type used in matrix calculations */
using ftype = double;

/// Vector types ///
/** vector type */
template<size_t S> using Vector1D = Eigen::Matrix<ftype, 1, S>;  // Eigen::Array3d;
/** Three-vector type */
// using ThreeVector = Vector1D<c_DimPosition>;  // Eigen::Array3d;
using ThreeVector = sct::kine::ThreeVector<>;
/** Lorentz-vector */
// using FourVector = Vector1D<c_DimMomentum>;  // Eigen::Array4d;
using FourVector = sct::kine::FourVector<>;

/// Matrix types ///
/** Matrix D1xD2 */
template<size_t D1, size_t D2> using Matrix = Eigen::Matrix<ftype, D1, D2>;
/** Square Matrix SxS */
template<size_t S> using SquareMatrix = Eigen::Matrix<ftype, S, S>;
/** 3x3 matrix */
using ThreeMatrix = SquareMatrix<c_DimPosition>;
/** 4x4 matrix */
using FourMatrix = SquareMatrix<c_DimMomentum>;

/** 7x7 covariance matrix for particle position (x, y, z) and momentum (E, px, py, pz) */
using ErrMatrix = SquareMatrix<c_DimMatrix>;
/** 4x4 covariance matrix for particle momentum (E, px, py, pz) */
using MomentumErrMatrix = SquareMatrix<c_DimMomentum>;
/** 3x3 covariance matrix for particle position (x, y, z) */
using PositionErrMatrix = SquareMatrix<c_DimPosition>;
/** 3x3 rotation matrix in 3D space */
using RotationMatrix = SquareMatrix<c_DimPosition>;
/**  5x5 jacobian matrix for the transport of the helix parameters
 *  Used in: Helix.h */
using HelixJacobian = SquareMatrix<c_DimHelJacob>;
/** 5x5 Covariance matrix for the five helix parameters.
 *  Indices correspond to the order d0, phi0, omega, z0, tanLambda
 *  Used in: UncertainHelix.h */
using HelixUncert = SquareMatrix<c_DimHelUncer>;
/** 6x6 Covariance matrix for position and momentum of the track at the perigee
 *  Used in: UncertainHelix.h */
using HelixUncertCart = SquareMatrix<c_DimHelUncerCart>;

}  // namespace sct
