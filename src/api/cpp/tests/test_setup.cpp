//
// Created by Eric Irrgang on 10/18/17.
//

#include "gmxpre.h"

#include <gtest/gtest.h>

#include "moduletest.h"

namespace
{

//! Test fixture for grompp
class GromppTest :
        public gmx::test::MdrunTestFixture
{
    public:
        //! Execute the trajectory writing test
        void runTest()
        {
            runner_.useTopGroAndNdxFromDatabase("spc-and-methanol");
            EXPECT_EQ(0, runner_.callGrompp());
        }
};

/* TODO Now that Verlet is the default, change the implementation
   of useEmptyMdpFile() to do that. */
std::string theMdpFile{
R"""(cutoff-scheme = Verlet
pull = yes
pull-ngroups = 2
pull-ncoords = 1
pull-group1-name = SOL
pull-group2-name = Methanol
pull-coord1-groups = 1 2
pull-coord1-type = umbrella
pull-coord1-geometry = distance
)"""};

/* This test ensures that an empty .mdp file (ie. all default values) works. */
TEST_F(GromppTest, EmptyMdpFileWorks)
{
    runner_.useStringAsMdpFile(theMdpFile);
    runTest();
}

} // end anonymous namespace
