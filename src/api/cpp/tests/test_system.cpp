//#include "atoms.h"

#include "testingconfiguration.h"

#include "gmxapi/gmxapi.h"
#include "gmxapi/system.h"
#include "gmxapi/md.h"
#include <gtest/gtest.h>

namespace
{

const auto filename = gmxapi::testing::sample_tprfilename;

TEST(ApiSystem, Construction)
{
    {   // Construction
        auto system = gmxapi::System();
    }   // Destruction

    auto system = gmxapi::fromTprFile(filename);
    ASSERT_TRUE(system != nullptr);
}

TEST(ApiSystem, Accessors)
{
    auto system = gmxapi::fromTprFile(filename);
//    ASSERT_TRUE(system->md() != nullptr);
//    ASSERT_NO_THROW(system->md()->info());
//    ASSERT_STREQ("Generic MDEngine object", system->md()->info().c_str());
//
//    ASSERT_TRUE(system->runner() != nullptr);

//    ASSERT_EQ(system->atoms()->x()->size(), 7);
//    ASSERT_TRUE(system->atoms() != nullptr);
//    ASSERT_TRUE(system->atoms()->x() != nullptr);
}

} // end anonymous namespace
