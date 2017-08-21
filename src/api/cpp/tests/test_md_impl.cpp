// Test the implementation classes not visible in the public API.

#include "gmxapi/md.h"
#include "../md-impl.h"
#include <gtest/gtest.h>

#include "gromacs/utility/keyvaluetree.h"

TEST(ApiMDInput, Construction)
{
    {
        auto mdInput = new gmxapi::MDInput();
        delete mdInput;
    }
    {
        auto mdInput = gmxapi::MDInput::from_tpr_file("topol.tpr");
        ASSERT_TRUE(mdInput.get() != nullptr);
    }
}

TEST(ApiMDInput, Accessors)
{
    auto mdInput = gmxapi::MDInput::from_tpr_file("topol.tpr");
    ASSERT_EQ(mdInput->nAtoms(), 6);
}
