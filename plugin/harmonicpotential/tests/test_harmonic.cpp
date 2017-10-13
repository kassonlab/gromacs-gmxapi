//
// Created by Eric Irrgang on 10/13/17.
//

#include <gtest/gtest.h>
#include "harmonicpotential.h"
#include "gromacs/pulling/vectortype.h"

namespace {

using gmx::detail::vec3;

TEST(HarmonicPotentialPlugin, Build)
{
    ASSERT_TRUE(true);
    ASSERT_FALSE(false);

    plugin::Harmonic puller;
}

TEST(HarmonicPotentialPlugin, ForceCalc)
{
    const vec3<real> e1{real(1), real(0), real(0)};
    const vec3<real> e2{real(0), real(1), real(0)};

    plugin::Harmonic puller;

    const vec3<real> expectedForce{real(1), real(1), real(1)};

    auto f = puller.calculateForce(e1, e2);
    ASSERT_EQ(f, expectedForce);
}

} // end anonymous namespace

int main(int argc, char* argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
