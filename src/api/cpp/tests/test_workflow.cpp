//
// Created by Eric Irrgang on 11/28/17.
//

#include "workflow.h"
#include "workflow-impl.h"

#include <memory>

#include "gmxapi/context.h"
#include "gmxapi/md/mdmodule.h"
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

// Create a work spec, then the implementation graph, then the container
TEST(ApiWorkflowImpl, Build)
{
    // Create work spec
    auto node = gmx::compat::make_unique<gmxapi::MDNodeSpecification>(filename);
    ASSERT_NE(node, nullptr);

    // Create key
    std::string key{"MD"};
    key.append(filename);

    // Create graph (workflow implementation object)
    gmxapi::Workflow::Impl impl;
    impl[key] = std::move(node);
    ASSERT_EQ(impl.count(key), 1);
    ASSERT_EQ(impl.size(), 1);

    // Create workflow container
    gmxapi::Workflow work{std::move(impl)};
}

TEST(ApiWorkflow, Creation)
{
    // Create from create() method(s)
    auto work = gmxapi::Workflow::create(filename);
    ASSERT_NE(work, nullptr);
}

TEST(ApiWorkflow, Accessors)
{
    auto work = gmxapi::Workflow::create(filename);
//    work->addNode()
//    work->getNode()
}

} // end anonymous namespace
