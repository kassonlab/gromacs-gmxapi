#include "gmxapi/runner.h"
#include "gmxapi/md/runnerstate.h"
#include "gmxapi/md.h"

#include "gromacs/compat/make_unique.h"

#include "gmxapi/system.h"
#include "gmxapi/context.h"
#include <gtest/gtest.h>

namespace
{

class DummyMD : public gmxapi::MDEngine
{
};

TEST(ApiRunner, Build)
{
    auto md = std::make_shared<gmxapi::MDEngine>();
    auto runnerBuilder = gmxapi::UninitializedMDRunnerState::Builder();
    runnerBuilder.mdEngine(md);
    runnerBuilder.topology(gmx::compat::make_unique<gmx_mtop_t>());
    runnerBuilder.state(gmx::compat::make_unique<t_state>());
    runnerBuilder.inputRecord(gmx::compat::make_unique<t_inputrec>());
    auto runner = runnerBuilder.build();
    auto session = runner->initialize(gmxapi::defaultContext());
    auto status = session->run();
    // Just make sure we made it this far...
    ASSERT_TRUE(!status.success());
}

TEST(ApiRunner, BasicMD)
{
//    const std::string filename = "topol.tpr";
    const std::string filename = "/Users/eric/build/gromacs-dev/src/api/cpp/tests/topol.tpr";
    auto system = gmxapi::fromTprFile(filename);

    {
        std::shared_ptr<gmxapi::Context> context = gmxapi::defaultContext();
        ASSERT_TRUE(context != nullptr);
        ASSERT_TRUE(system != nullptr);
        ASSERT_TRUE(system->runner() != nullptr);
        auto runner = system->runner();
        auto session = runner->initialize(context);
        ASSERT_TRUE(session != nullptr);
        ASSERT_NO_THROW(session->run());
//        ASSERT_NO_THROW(session->run(1000));
    }

}

}