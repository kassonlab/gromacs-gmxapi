// Copyright 2017
// Author: M. Eric Irrgang
#include <gtest/gtest.h>
#include "gromacs/restraint/manager.h"

namespace
{

class DummyRestraint: public gmx::IRestraintPotential
{
    public:
        ~DummyRestraint() override = default;

        gmx::PotentialPointData evaluate(gmx::Vector r1,
                                         gmx::Vector r2,
                                         double t) override
        {
            (void) r1;
            (void) r2;
            (void) t;
            return {};
        }

        void update(gmx::Vector v,
                    gmx::Vector v0,
                    double t) override
        { (void)v; (void)v0; (void)t; }

        std::vector<unsigned long> sites() const override
        {
            return std::vector<unsigned long>();
        }

        void bindSession(gmxapi::SessionResources *session) override
        {
            (void)session;
        }
};

TEST(RestraintManager, singleton)
{
    auto managerInstance = gmx::restraint::Manager::instance();
    ASSERT_TRUE(managerInstance.get() != nullptr);
}

TEST(RestraintManager, restraintList)
{
    auto managerInstance = gmx::restraint::Manager::instance();
    managerInstance->addToSpec(std::make_shared<DummyRestraint>(), "a");
    managerInstance->addToSpec(std::make_shared<DummyRestraint>(), "b");
    ASSERT_EQ(managerInstance->countRestraints(), 2);
    managerInstance->clear();
    ASSERT_EQ(managerInstance->countRestraints(), 0);
    managerInstance->addToSpec(std::make_shared<DummyRestraint>(), "c");
    managerInstance->addToSpec(std::make_shared<DummyRestraint>(), "d");
    ASSERT_EQ(managerInstance->countRestraints(), 2);
}

}
