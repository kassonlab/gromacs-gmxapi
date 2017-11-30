#include "gmxapi/system.h"

#include <cstdio>

#include <array>

#include "gmxapi/gmxapi.h"
#include "gmxapi/md.h"
#include "gmxapi/context.h"
#include "gmxapi/runner.h"
#include "gmxapi/status.h"
#include "gmxapi/md/runnerstate.h"

#include "gromacs/compat/make_unique.h"
#include "md-impl.h"
#include "session.h"
#include "system-impl.h"
#include "workflow.h"
#include "gromacs/utility.h"
#include "gromacs/mdtypes/tpxstate.h"
#include "gromacs/utility/keyvaluetree.h"
#include "programs/mdrun/runner.h"

namespace gmxapi
{

System::Impl::~Impl() = default;


// Constructor and destructor needs to be defined after Impl is defined so that we can
// use unique_ptr
System::System() :
    impl_ {gmx::compat::make_unique<System::Impl>()}
{
    assert(impl_ != nullptr);
}

std::unique_ptr<Session> System::launch(std::shared_ptr<Context> context)
{
    (void)context;
    auto session = gmx::compat::make_unique<Session>();
    return session;
}

std::unique_ptr<Session> System::launch()
{
    return launch(std::make_shared<Context>());
}

Status System::status()
{
    assert(impl_ != nullptr);
    return impl_->status();
}

System::System(std::unique_ptr<System::Impl>&& implementation) :
    impl_{std::forward<std::unique_ptr<System::Impl>>(implementation)}
{
    assert(impl_ != nullptr);
}

System::~System() = default;

System::System(System &&) noexcept = default;

System &System::operator=(System &&) noexcept = default;

std::unique_ptr<gmxapi::System> fromTprFile(std::string filename)
{
    // Confirm the file is readable and parseable and note unique identifying information
    // for when the work spec is used in a different environment.

    // Create a new Workflow instance.
    auto workflow = Workflow::create(filename);
    auto systemImpl = gmx::compat::make_unique<System::Impl>(std::move(workflow));
    if (systemImpl == nullptr)
    {

    }
    auto system = gmx::compat::make_unique<System>(std::move(systemImpl));

    // The TPR file has enough information for us to
    //  1. choose an MD engine
    //  2. Get structure information
    //  3. Get topology information
    //  4. Get a lot of simulation and runtime parameters, but not all.
    // It does not have enough information on its own to determine much about the
    // necessary computation environment. That comes from environment
    // introspection and user runtime options.

    // for what it's worth, we can choose/configure a builder based
    // on the sort of system we are building.

    // For now have very limited execution environment abstraction
    // or flexibility.
//    builder->defaultContext();

//    // Build MDEngine member
//    auto mdBuilder = MDProxy().builder();
//    assert(mdBuilder != nullptr);
//    std::shared_ptr<MDEngine> md = mdBuilder->build();
//    assert(md != nullptr);
//    assert(!md->info().empty());
//
//    auto runnerBuilder = UninitializedMDRunnerState::Builder();
//    runnerBuilder.mdEngine(md);
//    auto tpxState = gmx::TpxState::initializeFromFile(filename);
//    assert(!tpxState->isDirty());
//    assert(tpxState->isInitialized());
//    runnerBuilder.tpxState(std::move(tpxState));
//    auto runner = runnerBuilder.build();
//    assert(runner != nullptr);

    // One way to get around the required CLI argument is to make sure topol.tpr is in the working directory...

//
//    auto builder = gmx::compat::make_unique<System::Builder>();
//    builder->mdEngine(md);
//    builder->runner(std::move(runner));
//
//    auto system = builder->build();

    return system;
}


System::Impl::Impl() :
    context_{std::make_shared<Context>()},
    workflow_{nullptr},
    status_{gmx::compat::make_unique<Status>()}
{
    assert(context_ != nullptr);
    assert(status_ != nullptr);
}

Status System::Impl::status() const
{
    return *status_;
}

System::Impl::Impl(std::unique_ptr<gmxapi::Workflow> &&workflow) :
    context_{nullptr},
    workflow_{std::move(workflow)},
    status_{gmx::compat::make_unique<Status>(true)}
{
    //What should context_ be initialized to? How about the default local context.
    //assert(context_ != nullptr);
    assert(workflow_ != nullptr);
    assert(status_ != nullptr);
}

} // end namespace gmxapi
