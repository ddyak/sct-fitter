#pragma once

#include <Eigen/Core>

using Eigen::MatrixXd;
using Eigen::VectorXd;

namespace sct::ana {

/** Class to store and manage fitparams (statevector) */
class FitParams {
public:
    /** constructor */
    explicit FitParams(const int dim);

    /** copy constructor */
    FitParams(const FitParams& rhs)
      : m_dim(rhs.m_dim),
        // m_chiSquare(rhs.m_chiSquare),
        // m_nConstraints(toCopy.m_nConstraints),
        // m_dimensionReduction(toCopy.m_dimensionReduction),
        m_globalState(VectorXd(rhs.m_globalState)),
        m_globalCovariance(MatrixXd(rhs.m_globalCovariance))
    {}

    /** Assignment operator. */
    FitParams& operator=(const FitParams& rhs) {
      m_dim = rhs.m_dim;
    //   m_chiSquare = rhs.m_chiSquare;
    //   m_nConstraints = other.m_nConstraints;
    //   m_dimensionReduction = other.m_dimensionReduction;
      m_globalState = rhs.m_globalState;
      m_globalCovariance = rhs.m_globalCovariance;
      return *this;
    }

    /** getter for the states covariance */
    MatrixXd& getCovariance() { return m_globalCovariance; }

    /** const getter for the states covariance */
    const MatrixXd& getCovariance() const { return m_globalCovariance; }

    /** getter for the fit parameters/statevector */
    VectorXd& getStateVector() { return m_globalState; }

    /** const getter for the fit parameters/statevector */
    const VectorXd& getStateVector() const { return m_globalState; }

    /** reset the staevector */
    void resetStateVector();

    /** reset the staevector */
    void resetCovariance();

    /** test if the covariance makes sense */
    bool testCovariance() const;

    /** get the states dimension */
    int getDimensionOfState() const { return m_dim; }

    /** increment global chi2 */
    void addChiSquare(double chisq, int nconstraints) {
        m_chiSquare += chisq;
        m_nConstraints += nconstraints;
    }

    /** reset chi2 */
    void resetChiSquare() {
        m_chiSquare = 0;
        m_nConstraints = 0;
    }
    
private:
    /** dimension of statevector */
    int m_dim;

    /** chi2 */
    double m_chiSquare;

    /** number of conatraints */
    int m_nConstraints;
    
    /** vector holding all parameters of this fit */
    Eigen::VectorXd m_globalState;

    /** covariance of the global state */
    Eigen::MatrixXd m_globalCovariance;
};

}