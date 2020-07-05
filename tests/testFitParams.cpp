#include "gtest/gtest.h"

#include "treefitter/FitParams.h"

using namespace sct::ana;

TEST(FitParams, Getters) {
    std::size_t size = 20;
    FitParams fit_params(size);
    EXPECT_EQ(fit_params.getStateVector().rows(), size);
    EXPECT_EQ(fit_params.getCovariance().cols(), size);
}

TEST(FitParams, Resets) {
    std::size_t size = 20;
    FitParams fit_params(size);

    fit_params.resetCovariance();
    fit_params.resetStateVector();

    EXPECT_EQ(fit_params.getCovariance(), Eigen::MatrixXd::Zero(size, size));
    EXPECT_EQ(fit_params.getStateVector(), Eigen::VectorXd::Zero(size));
}