#include "gromacs/compat/make_unique.h"
#include "runner-impl.h"
#include "gmxapi/gmxapi.h"
#include <gtest/gtest.h>

#include "gromacs/fileio/tpxio.h"
#include "gromacs/mdlib/integrator.h"
#include "gromacs/mdlib/mdrun.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/keyvaluetree.h"
#include "programs/mdrun/md.h"


TEST(ApiRunner, BasicMD)
{
    const std::string filename = "topol.tpr";
    // Configure an MD simulation from tpr file
    auto              runner = gmx::compat::make_unique<gmxapi::RunnerImpl>(filename);
    // Get the initial atom positions
    auto              x = runner->getX();
    // Check for success
    ASSERT_NE(x, nullptr);
    // Check for expected result. Note this is padded. Not sure what to do about that.
    ASSERT_EQ(x->size(), 7);
    // Stash one of the atom positions
    auto atom1 = (*x)[1];
    // Run from initial configuration to 1000 steps
    runner->run(1000);
    // Continue run for an additional 1000 steps
    runner->run(1000);
    // Get new atom positions
    x = runner->getX();
    // Check for success
    ASSERT_NE(x, nullptr);
    // Make sure atom moved
    ASSERT_NE((*x)[1][0], atom1[0]);
    // Clean up
    runner->close();
    runner.reset();
    ASSERT_NE(x, nullptr);
}
