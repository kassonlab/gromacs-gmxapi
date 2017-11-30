//
// Created by Eric Irrgang on 11/28/17.
//

#include "workflow.h"
#include "workflow-impl.h"

#include <memory>

#include "gmxapi/context.h"
#include "gmxapi/runner.h"
#include "gmxapi/md/mdmodule.h"
#include "gmxapi/md/runnerstate.h"
#include "gmxapi/md.h"
#include "gmxapi/status.h"
#include "gmxapi/system.h"

#include "gromacs/compat/make_unique.h"
#include "gromacs/mdtypes/tpxstate.h"
#include "gromacs/restraint/restraintpotential.h"
#include "gromacs/utility/arrayref.h"

#include "testingconfiguration.h"
#include <gtest/gtest.h>

namespace
{

const auto filename = gmxapi::testing::sample_tprfilename;

TEST(ApiWorkflowImpl, Build)
{
    // Create a work spec, then the implementation graph, then the container
    auto node = gmx::compat::make_unique<gmxapi::MDNodeSpecification>(filename);
    ASSERT_NE(node, nullptr);
}

TEST(ApiWorkflow, Creation)
{
    // Create from create() method(s)
}

} // end anonymous namespace
