#include "gmxapi/context.h"
#include "gmxapi/gmxapi.h"
#include <gtest/gtest.h>

TEST(ApiContext, Construction)
{
    {
        auto context = new gmxapi::Context();
        context->initialize();
        ASSERT_TRUE(context->isInitialized());
        context->deinitialize();
        ASSERT_FALSE(context->isInitialized());
        delete context;
    }

    // {
    // auto inputRecord = gmxapi::MDInput::from_tpr_file("topol.tpr");
    // auto context = gmxapi::Context::fromInputRec(inputRecord);
    // // Right now, the context is initialized when creating from TPR.
    // ASSERT_TRUE(context.isInitialized());
    // }
}
