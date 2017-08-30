// Mimic a simple invocation of `gmx mdrun` from a TPR file.

#include <memory>

#include "gmxapi/gmxapi.h"
#include "gmxapi/system.h"
#include <gtest/gtest.h>

using std::shared_ptr;

namespace
{


// Test 1: create objects in the traditional scopes

TEST(ApiRunner, Basic)
{
//    ASSERT_TRUE(gmxapi::test_impl_from_file("membrane.tpr"));
    {
        auto  system = gmxapi::from_tpr_file("topol.tpr");
        ASSERT_TRUE(system != nullptr);
    }
    ASSERT_TRUE(true);
}

//TEST(ApiRunner, Run)
//{
//    auto system = gmxapi::from_tpr_file("topol.tpr");
//    ASSERT_EQ(system->run().success(), true);
//}

// Test 2: create objects with handles outside of the execution flow.


} // end anonymous namespace
