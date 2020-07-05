#include "gtest/gtest.h"

#include "external/EvtPdlDatabase.h"

using namespace sct::comm;
using namespace std::string_literals;

TEST(EvtPdlDatabase, electron) {
    int pdgCode = 11; // e-
    EXPECT_EQ(EvtPdlDatabase::Instance().isKnownId(pdgCode), true);
    EXPECT_EQ(EvtPdlDatabase::Instance().isKnownName("e-"s), true);
    EXPECT_EQ(EvtPdlDatabase::Instance().name(pdgCode), "e-");
    EXPECT_NEAR(EvtPdlDatabase::Instance().mass(pdgCode), 0.000510999, 1e-8);
    EXPECT_NEAR(EvtPdlDatabase::Instance().lifetime(pdgCode), 0, 1e-8);
}
