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
    constexpr vec3<real> zerovec = gmx::detail::make_vec3<real>(0, 0, 0);
    // define some unit vectors
    const vec3<real> e1{real(1), real(0), real(0)};
    const vec3<real> e2{real(0), real(1), real(0)};
    const vec3<real> e3{real(0), real(0), real(1)};

    plugin::Harmonic puller;

    // When input vectors are equal, output vector is meaningless and magnitude is set to zero.
    ASSERT_EQ(real(0.0), norm(puller.calculateForce(e1, e1)));

    // Default equilibrium distance is 1.0, so force should be zero when norm(r12) == 1.0.
    ASSERT_EQ(zerovec, puller.calculateForce(zerovec, e1));
    ASSERT_EQ(zerovec, puller.calculateForce(e1, zerovec));
    ASSERT_EQ(zerovec, puller.calculateForce(e1, 2*e1));

    // -kx should give vector (1, 0, 0) when vector r1 == r2 - (2, 0, 0)
    ASSERT_EQ(real(1), puller.calculateForce(-2*e1, zerovec).x);
    ASSERT_EQ(e1, puller.calculateForce(-2*e1, zerovec));

    // -kx should give vector (-2, 0, 0) when vector r1 == r2 + (2, 0, 0)
    ASSERT_EQ(-2*e1, puller.calculateForce(2*e1, -e1));
}

} // end anonymous namespace

int main(int argc, char* argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
