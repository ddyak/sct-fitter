#include "gtest/gtest.h"

#include "external/ThreeVector.h"


TEST(ThreeVector, Example) {
	sct::kine::ThreeVector<double> vec(1,2,3);
	EXPECT_EQ(vec.sum(), 6);
}