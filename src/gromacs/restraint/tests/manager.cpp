// Copyright 2017
// Author: M. Eric Irrgang
#include <gtest/gtest.h>
#include "gromacs/restraint/manager.h"

namespace
{

TEST(RestraintManager, singleton)
{
    ASSERT_TRUE(true);
    auto managerInstance = gmx::restraint::Manager::instance();
    ASSERT_TRUE(managerInstance.get() != nullptr);
}

}
