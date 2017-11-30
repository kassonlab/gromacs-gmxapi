//
// Created by Eric Irrgang on 7/31/17.
//

#include "gmxapi/runner.h"

#include <string>
#include <cassert>
#include "gromacs/utility.h"
#include "programs/mdrun/runner.h"
#include "gromacs/compat/make_unique.h"
#include "gromacs/mdtypes/tpxstate.h"
#include "gromacs/mdtypes/inputrec.h"

#include "gmxapi/exceptions.h"
#include "gmxapi/status.h"
#include "gmxapi/md.h"
#include "gmxapi/md/runnerstate.h"
#include "gmxapi/md/mdmodule.h"

namespace gmxapi
{

// Delegate construction to RunnerProxy::RunnerProxy(std::shared_ptr<MDProxy> md)
RunnerProxy::RunnerProxy() :
    RunnerProxy{std::shared_ptr<MDProxy>(nullptr)}
{};

RunnerProxy::RunnerProxy(std::shared_ptr<MDProxy> md) :
    module_{std::move(md)},
    instanceState_{std::make_shared<EmptyMDRunnerState>()}
{
}

void RunnerProxy::setState(std::shared_ptr<IMDRunner> state)
{
    instanceState_ = std::move(state);
}

void RunnerProxy::registerMDBuilder(std::unique_ptr<MDBuilder> builder)
{
    (void)builder;
}

Status RunnerProxy::run()
{
    assert(instanceState_ != nullptr);
    return instanceState_->run();
}

std::shared_ptr<IMDRunner> RunnerProxy::initialize(std::shared_ptr<Context> context)
{
    std::shared_ptr<IMDRunner> initializedRunner = instanceState_->initialize(context);

    instanceState_.swap(initializedRunner);
    return instanceState_;
}

void RunnerProxy::setRestraint(std::shared_ptr<gmxapi::MDModule> restraint)
{
    assert(instanceState_ != nullptr);
    instanceState_->setRestraint(std::move(restraint));
}

void EmptyMDRunnerState::registerMDBuilder(std::unique_ptr<MDBuilder> builder)
{
    (void)builder;
    // nothing to bind to
    throw ProtocolError("EmptyMDRunnerState has nothing to bind to");
}

Status EmptyMDRunnerState::run()
{
    // Nothing to run...
    throw ProtocolError("EmptyMDRunnerState has nothing to run");
    return Status();
}

std::shared_ptr<IMDRunner> EmptyMDRunnerState::initialize(std::shared_ptr<Context> context)
{
    (void)context;
    // Runner proxy should have been configured by a builder with something like UninitializedMDRunnerState
    throw ProtocolError("EmptyMDRunnerState cannot be initialized.");
    return nullptr;
}

class UninitializedMDRunnerState::Impl
{
    public:
        std::shared_ptr<MDEngine> mdProxy_;
        std::shared_ptr<gmx::TpxState> tpxState_;
};

Status UninitializedMDRunnerState::run()
{
    // We could be more helpful about suggesting the user initialize the runner first...
    throw ProtocolError("UninitializedMDRunnerState cannot be run.");
    return Status();
}


void UninitializedMDRunnerState::registerMDBuilder(std::unique_ptr<MDBuilder> builder)
{
    (void)builder;
}

std::shared_ptr<IMDRunner> UninitializedMDRunnerState::initialize(std::shared_ptr<Context> context)
{
    RunningMDRunnerState::Builder builder{};
    builder.tpxState(impl_->tpxState_);

    std::shared_ptr<IMDRunner> initializedRunner = builder.build();
    return initializedRunner;
}

UninitializedMDRunnerState::UninitializedMDRunnerState() :
    impl_{gmx::compat::make_unique<UninitializedMDRunnerState::Impl>()}
{
}

UninitializedMDRunnerState::~UninitializedMDRunnerState() = default;

UninitializedMDRunnerState::Builder::Builder() :
    runner_{nullptr}
{
    // Cannot use make_unique() because constructor is private
    try
    {
        runner_.reset(new UninitializedMDRunnerState());
    }
    catch(const std::exception& e)
    {
        // How should we report memory errors?
        throw e;
    }
}

std::unique_ptr<UninitializedMDRunnerState> UninitializedMDRunnerState::Builder::build()
{
    std::unique_ptr<UninitializedMDRunnerState> runnerState;
    if (runner_->impl_->mdProxy_ && runner_->impl_->tpxState_)
    {
        runnerState.swap(runner_);
    }
    else
    {
        // Todo: codify build protocol and/or provide more helpful error-checking
        throw(ProtocolError("Builder has insufficient input for a valid product."));
    }
    return runnerState;
}

UninitializedMDRunnerState::Builder &UninitializedMDRunnerState::Builder::mdEngine(std::shared_ptr<MDEngine> md)
{
    assert(md != nullptr);
    assert(runner_->impl_ != nullptr);
    runner_->impl_->mdProxy_ = std::move(md);
    return *this;
}

UninitializedMDRunnerState::Builder &
UninitializedMDRunnerState::Builder::tpxState(std::shared_ptr<gmx::TpxState> input)
{
    assert(input != nullptr);
    assert(runner_->impl_ != nullptr);
    runner_->impl_->tpxState_ = std::move(input);
    return *this;
}

UninitializedMDRunnerState::Builder::~Builder() = default;

class RunningMDRunnerState::Impl
{
    public:
        std::shared_ptr<MDEngine> mdProxy_;
        std::shared_ptr<gmx::Mdrunner> runner_;
        // Make sure we don't have incompatible types getting implicitly converted behind our backs.
        decltype(t_inputrec::nsteps) nSteps_;

        Impl();
        Status run();
};

RunningMDRunnerState::Impl::Impl() :
    mdProxy_{nullptr},
    runner_{nullptr},
    nSteps_{0}
{
}

Status RunningMDRunnerState::Impl::run()
{
    assert(runner_ != nullptr);
    Status status{};
    // Todo: check the number of steps to run
    if (runner_->mdrunner() == 0)
    {
        status = true;
    }
    return status;
}

RunningMDRunnerState::~RunningMDRunnerState() = default;

RunningMDRunnerState::RunningMDRunnerState() :
    impl_{gmx::compat::make_unique<RunningMDRunnerState::Impl>()}
{
}

Status RunningMDRunnerState::run()
{
    if (impl_ == nullptr)
    {
        throw gmxapi::ProtocolError("Runner not initialized.");
    }
    return impl_->run();
}


std::shared_ptr<IMDRunner> RunningMDRunnerState::initialize(std::shared_ptr<Context> context)
{
    (void)context;
    // Should we reinitialize somehow?
    throw gmxapi::NotImplementedError("Initializing a running Mdrunner is not defined.");
    return std::shared_ptr<IMDRunner>();
}

void RunningMDRunnerState::registerMDBuilder(std::unique_ptr<MDBuilder> builder)
{
    (void)builder;
    // implement the runner--mdengine binding protocol
}

// Implement the protocol
// protocol:
/* void gmxapi::IMDRunner::register(std::shared_ptr<gmxapi::MDModule> module)
 * {
 *      auto runner = impl_->runner_;
 *      auto restraint = module->getRestraint();
 *      runner->addPullPotential(restraint);
 * };
 */
void RunningMDRunnerState::setRestraint(std::shared_ptr<gmxapi::MDModule> module)
{
    assert(impl_ != nullptr);
    assert(impl_->runner_ != nullptr);

    // Follow the protocol to register a gmxapi::MDModule with a gmxapi::IMDRunner
    // by passing an gmx::IRestraintPotential to a gmx::MdRunner
    // \todo we should be registering a Spec or Factory instead of creating the IRestraint now.
    // Note that thread-MPI threads have not yet been spawned at this phase as of 20 November 2017
    auto runner = impl_->runner_;
    auto restraint = module->getRestraint();
    if(restraint != nullptr)
    {
        runner->addPullPotential(restraint,
                                 module->name());
    }
}

//void RunningMDRunnerState::addModule(std::shared_ptr<gmx::IMDModule> module)
//{
//    assert(impl_ != nullptr);
//    assert(impl_->runner_ != nullptr);
//    impl_->runner_->addModule(module);
//}


RunningMDRunnerState::Builder::Builder() :
    runner_{nullptr}
{
    // make_unique() not available for private constructors
    std::unique_ptr<RunningMDRunnerState> newRunner{new RunningMDRunnerState()};
    runner_ = std::move(newRunner);
}

std::unique_ptr<RunningMDRunnerState> RunningMDRunnerState::Builder::build()
{
    std::unique_ptr<RunningMDRunnerState::Impl> runnerImpl{nullptr};

    if (tpxState_ != nullptr)
    {
        auto newMdrunner = gmx::compat::make_unique<gmx::Mdrunner>();
        newMdrunner->setTpx(tpxState_);
        // Right now we need to borrow the CLII code...
        newMdrunner->initFromAPI();

        runnerImpl = gmx::compat::make_unique<RunningMDRunnerState::Impl>();
        runnerImpl->runner_ = std::move(newMdrunner);
    }
    runner_->impl_ = std::move(runnerImpl);

    // Not implemented
    //runner_->mdEngine(md_);

    std::unique_ptr<RunningMDRunnerState> activeRunner{nullptr};
    if ((runner_->impl_ != nullptr) /* && other validity checks... */)
    {
        activeRunner = std::move(runner_);
    }

    return activeRunner;
}

RunningMDRunnerState::Builder &RunningMDRunnerState::Builder::tpxState(std::shared_ptr<gmx::TpxState> input)
{
    assert(input != nullptr);
    tpxState_ = std::move(input);
    return *this;
}

RunningMDRunnerState::Builder::~Builder() = default;

void IMDRunner::setRestraint(std::shared_ptr<gmxapi::MDModule> restraint)
{
        (void)restraint;
        throw ::gmxapi::ProtocolError("setRestraint not implemented for this class.");
}

} // end namespace gmxapi
