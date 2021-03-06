/**************************************************************************
 *
 * BASF2 (Belle Analysis Framework 2)                                     *
 * Copyright(C) 2018 - Belle II Collaboration                             *
 *                                                                        *
 * Author: The Belle II Collaboration                                     *
 * Contributor: Wouter Hulsbergen, Francesco Tenchini, Jo-Frederik Krohn  *
 *                                                                        *
 * This software is provided "as is" without any warranty.                *
 **************************************************************************/

#include "Projection.h"

using namespace sct::ana;


Projection::Projection(int dimP, int dimC) :
  m_H(dimC, dimP),
  m_V(dimC, dimC),
  m_residual(dimC),
  m_dimCov(dimC),
  m_dimProj(dimP),
  m_offset(0) {}

void Projection::resetProjection()
{
  m_residual = Eigen::Matrix < double, -1, 1, 0, 5, 1 >::Zero(m_dimCov);
  m_H = Eigen::Matrix < double, -1, -1, 0, 5, MAX_MATRIX_SIZE >::Zero(m_dimCov, m_dimProj);
  m_V = Eigen::Matrix < double, -1, -1, 0, 5, 5 >::Zero(m_dimCov, m_dimCov);
  m_offset = 0;
}
