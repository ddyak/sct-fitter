#include "FitParams.h"

using namespace sct::ana;


FitParams::FitParams(const int dim):
    m_dim(dim),
    m_chiSquare(1e10),
    m_nConstraints(0),
    // m_dimensionReduction(0),
    // m_nConstraintsVec(dim, 0),
    m_globalState(dim),
    m_globalCovariance(dim, dim)
{
    resetStateVector();
    resetCovariance();
}

void FitParams::resetStateVector()
{
    m_globalState = VectorXd::Zero(m_dim);
}

void FitParams::resetCovariance()
{
    m_globalCovariance = MatrixXd::Zero(m_dim, m_dim);
    // std::fill(m_nConstraintsVec.begin(), m_nConstraintsVec.end(), 0);
    m_chiSquare = 0;
    m_nConstraints = 0;
}

bool FitParams::testCovariance() const
{
    for (int row = 0; row < m_dim; ++row) {
        if (m_globalCovariance(row, row) <= 0) 
            return false;
    }
    return true;
}