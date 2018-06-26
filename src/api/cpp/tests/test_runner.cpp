#include "testingconfiguration.h"

#include <memory>

#include "gmxapi/context.h"
#include "gmxapi/md/mdmodule.h"
#include "gmxapi/md.h"
#include "gmxapi/session.h"
#include "gmxapi/status.h"
#include "gmxapi/system.h"

#include "gromacs/compat/make_unique.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/iforceprovider.h"
#include "gromacs/mdtypes/imdmodule.h"
#include "gromacs/mdtypes/imdpoptionprovider.h"
#include "gromacs/mdtypes/imdoutputprovider.h"
#include "gromacs/mdtypes/tpxstate.h"
#include "gromacs/restraint/restraintpotential.h"
#include "gromacs/utility/arrayref.h"

#include <gtest/gtest.h>

namespace
{

const auto filename = gmxapi::testing::sample_tprfilename;

class DummyMDModule final : public gmx::IMDModule
{
    private:
        class OptionProvider : public gmx::IMdpOptionProvider
        {};
        std::shared_ptr<OptionProvider> optionprovider{std::make_shared<OptionProvider>()};

        class OutputProvider : public gmx::IMDOutputProvider
        {};
        std::shared_ptr<OutputProvider> outputprovider{std::make_shared<OutputProvider>()};

        class ForceProvider : public gmx::IForceProvider
        {
            public:
                void calculateForces(const t_commrec *cr,
                                     const t_mdatoms *mdatoms,
                                     const matrix box,
                                     double t,
                                     const rvec *x,
                                     gmx::ArrayRef<gmx::RVec> force)
                override
                {
                    force_called++;
                };

                unsigned int force_called{0};
        };
        std::shared_ptr<ForceProvider> forceprovider{std::make_shared<ForceProvider>()};

        gmx::IForceProvider* getForceProvider()
        {
            return forceprovider.get();
        };
    public:
        gmx::IMdpOptionProvider *mdpOptionProvider() override
        {
            return optionprovider.get();
        }

        gmx::IMDOutputProvider *outputProvider() override
        {
            return outputprovider.get();
        }

        void initForceProviders(ForceProviders *forceProviders) override
        {
            forceProviders->addForceProvider(getForceProvider());
        }

        unsigned int force_called() { return forceprovider->force_called; };
};

TEST(ApiRunner, BasicMD)
{

    auto system = gmxapi::fromTprFile(filename);

    {
        std::shared_ptr<gmxapi::Context> context = gmxapi::defaultContext();
        ASSERT_TRUE(context != nullptr);
        ASSERT_TRUE(system != nullptr);
        gmxapi::MDArgs args;
        args.emplace_back("-nsteps");
        args.emplace_back("10");
        context->setMDArgs(args);
        auto session = system->launch(context);
        ASSERT_TRUE(session != nullptr);
        gmxapi::Status status;
        ASSERT_NO_THROW(status = session->run());
//        ASSERT_NO_THROW(session->run(1000));
        ASSERT_TRUE(status.success());
        status = session->close();
        ASSERT_TRUE(status.success());
    }
}

TEST(ApiRunner, ContinuedMD)
{
    // Run a simulation, then extend the target number of steps and continue the simulation
    auto system = gmxapi::fromTprFile(filename);

    {
        std::shared_ptr<gmxapi::Context> context = gmxapi::defaultContext();

        {
            ASSERT_TRUE(context != nullptr);
            ASSERT_TRUE(system != nullptr);
            gmxapi::MDArgs args;
            args.emplace_back("-nsteps");
            args.emplace_back("20");
            context->setMDArgs(args);
            auto session = system->launch(context);
            ASSERT_TRUE(session != nullptr);
            gmxapi::Status status;
            ASSERT_NO_THROW(status = session->run());
            ASSERT_TRUE(status.success());
            ASSERT_NO_THROW(status = session->close());
            ASSERT_TRUE(status.success());
        }

        // Reuse the context. Add MD parameters. Run a new session extending the previous trajectory.
        {
            gmxapi::MDArgs args;
            args.emplace_back("-nsteps");
            args.emplace_back("20");
            context->setMDArgs(args);
            auto session = system->launch(context);
            ASSERT_TRUE(session != nullptr);
            gmxapi::Status status;
            ASSERT_NO_THROW(status = session->run());
            ASSERT_TRUE(status.success());
            ASSERT_NO_THROW(status = session->close());
            ASSERT_TRUE(status.success());
        }
    }
}

} // end anonymous namespace
