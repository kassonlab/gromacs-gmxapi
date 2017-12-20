#include "testingconfiguration.h"

#include <memory>

#include "gmxapi/context.h"
#include "gmxapi/md/mdmodule.h"
#include "gmxapi/md.h"
#include "gmxapi/session.h"
#include "gmxapi/status.h"
#include "gmxapi/system.h"

#include "gromacs/math/vectypes.h"
#include "gromacs/restraint/restraintpotential.h"
#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/arrayref.h"

#include <gtest/gtest.h>

namespace
{

const auto filename = gmxapi::testing::sample_tprfilename;

class NullRestraint : public gmx::IRestraintPotential
{
    public:
        gmx::PotentialPointData evaluate(gmx::Vector r1,
                                         gmx::Vector r2,
                                         double t) override
        {
            return {};
        }

        std::array<unsigned long, 2> sites() const override
        {
            return {{0,1}};
        }
};

class SimpleApiModule : public gmxapi::MDModule
{
    public:
        const char *name() override
        {
            return "SimpleApiModule";
        }

        std::shared_ptr<gmx::IRestraintPotential> getRestraint() override
        {
            auto restraint = std::make_shared<NullRestraint>();
            return restraint;
        }
};

// This testing is currently performed elsewhere.
//TEST(ApiRestraint, MdAndNullPlugin)
//{
//
//    {
//        // \todo move slow validation tests to a separate testing suite.
//        // Automate evaluation of the results of these test.
//        std::string waterfile = "water.tpr";
//        auto system = gmxapi::fromTprFile(waterfile);
//
//        auto module = std::make_shared<SimpleApiModule>();
//        assert(module != nullptr);
//        auto status = system->setRestraint(module);
//        ASSERT_TRUE(status.success());
//
//        std::shared_ptr<gmxapi::Context> context = gmxapi::defaultContext();
//
//        auto session = system->launch(context);
//        ASSERT_TRUE(session->isOpen());
//
//        ASSERT_NO_THROW(status = session->run());
////        ASSERT_TRUE(module->force_called() > 0);
////        ASSERT_NO_THROW(session->run(1000));
//        ASSERT_TRUE(status.success());
//        ASSERT_TRUE(session->isOpen());
//        status = session->close();
//        ASSERT_TRUE(status.success());
//    }
//
//}


} // end anonymous namespace
