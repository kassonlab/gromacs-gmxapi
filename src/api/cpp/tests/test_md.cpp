#include "gmxapi/md.h"
#include <gtest/gtest.h>

#include "gmxapi/status.h"
#include "md-impl.h"
#include "gromacs/utility/keyvaluetree.h"

TEST(ApiModuleMD, Construction)
{
//    {   // Check default construction and destruction
//        gmxapi::MDProxy proxy{};
//    }
    // Helper function should return a non-null unique_ptr
//    ASSERT_TRUE(module.get() != nullptr);
//    ASSERT_EQ(module->info(), "MDStatePlaceholder initialized with filename: \"topol.tpr\"\n");
}

TEST(ApiModuleMD, Build)
{
//    auto mdBuilder = gmxapi::MDProxy().builder();
//    ASSERT_TRUE(mdBuilder != nullptr);
//    auto md = mdBuilder->build();
//    ASSERT_TRUE(md != nullptr);
//    ASSERT_STREQ("Generic MDEngine object", md->info().c_str());
}

TEST(ApiModuleMD, Binding)
{
//    // Register MD with a runner
//    auto module = gmxapi::mdFromTpr("topol.tpr");
//    auto runner = MyRunner();
//    module->bind(&runner);
//    runner.run();
}
