#include <Eigen/Core>

using Eigen::MatrixXd;
using Eigen::VectorXd;

namespace sct::ana {

class FitParams {
public:
    /** constructor */
    explicit FitParams(const int dim);

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

private:
    /** dimension of statevector */
    int m_dim;
    
    /** vector holding all parameters of this fit */
    Eigen::VectorXd m_globalState;

    /** covariance of the global state */
    Eigen::MatrixXd m_globalCovariance;
};

}