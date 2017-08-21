#include "atoms.h"
#include "gmxapi/gmxapi.h"
#include "gmxapi/system.h"
#include <gtest/gtest.h>

TEST(ApiSystem, Construction)
{
    {   // Construction
        auto system = gmxapi::System();
    }   // Destruction

    auto system = gmxapi::from_tpr_file("topol.tpr");
    ASSERT_TRUE(system != nullptr);
//    ASSERT_EQ(system->atoms()->x()->size(), 7);
}

TEST(ApiSystem, Accessors)
{
    auto system = gmxapi::from_tpr_file("topol.tpr");
//    ASSERT_TRUE(system->atoms() != nullptr);
//    ASSERT_TRUE(system->atoms()->x() != nullptr);
}
