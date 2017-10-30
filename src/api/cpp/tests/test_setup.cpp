//
// Created by Eric Irrgang on 10/18/17.
//

#include "gmxpre.h"

#include <gtest/gtest.h>
#include <fstream>

#include "moduletest.h"

namespace
{

static const std::string atomicConfiguration{R"delimiter(spc-and-methanol
    6
    1MeOH   Me1    1   1.965   1.457   1.206 -0.8621 -0.4370 -0.5161
    1MeOH    O2    2   1.978   1.412   1.078 -0.0083 -0.4996 -0.6302
    1MeOH    H3    3   1.903   1.451   1.025 -0.2944 -1.5927 -1.1407
    2SOL     OW    4   1.559   1.516   0.709  0.7434  0.8949  1.0702
    2SOL    HW1    5   1.502   1.495   0.789  0.5852 -0.0293  0.7108
    2SOL    HW2    6   1.501   1.532   0.630  0.9009  1.8670  1.1441
   30.0000   30.0000   30.0000
)delimiter"};

//! Test fixture for grompp
class GromppTest :
        public gmx::test::MdrunTestFixture
{
    public:
        void runGrompp()
        {
            runner_.useTopGroAndNdxFromDatabase("spc-and-methanol");
            {
                std::fstream groFile(runner_.groFileName_, groFile.trunc | groFile.out);
                groFile << atomicConfiguration;
            }
            EXPECT_EQ(0, runner_.callGrompp());
        }

        void runMD()
        {
            EXPECT_EQ(0, runner_.callMdrun());
        }
};

/* Initial target distance is `pull-coord1-init`. Rate is in nm / ps.
 *
 * Apply stochastic dynamics (Langevin?) at very low temperature. Pull from 0.6 to 10.6nm
 * distance over 10ps. This is a small system, but use a big box.
 */
std::string theMdpFile{
R"""(cutoff-scheme = Verlet
integrator = sd
nsteps = 10000
tc-grps = System
tau-t = 0.1
ref-t = 10
pull = yes
pull-ngroups = 2
pull-ncoords = 1
pull-group1-name = SOL
pull-group2-name = Methanol
pull-coord1-groups = 1 2
pull-coord1-type = umbrella
pull-coord1-geometry = distance
pull-coord1-k = 1000
pull-coord1-init = 0.6
pull-coord1-rate = 1
)"""};

TEST_F(GromppTest, MdpFileWorks)
{
    runner_.useStringAsMdpFile(theMdpFile);
    runGrompp();
    runMD();
}

} // end anonymous namespace
