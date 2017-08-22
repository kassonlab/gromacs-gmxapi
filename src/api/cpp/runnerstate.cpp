//
// Created by Eric Irrgang on 7/31/17.
//

#include "gmxapi/runner.h"
#include "api/cpp/gmxapi/md/runnerstate.h"
#include "gmxapi/exceptions.h"

#include <string>
#include <cassert>
#include <programs/mdrun/runner.h>

#include "gmxapi/md.h"
#include "gromacs/compat/make_unique.h"

namespace gmxapi
{

/*!
 * \brief State implementation for first draft hack.
 *
 * Uses an outdated cut-and-pasted mdrunner with minor encapsulation and a filename parameter.
 */
//RunnerImplState::RunnerImplState(std::string filename) :
//                impl_{std::make_shared<gmxapi::RunnerImpl>(filename)},
//                mdproxy_{gmxapi::mdFromTpr(filename)}
//        {};
//
//RunnerImplState::~RunnerImplState()
//        {
//        }
//
//RunnerImplState::RunnerImplState(std::shared_ptr<RunnerProxy> owner, std::shared_ptr<MDProxy> mdproxy) :
//                impl_{nullptr},
//                mdproxy_{std::move(mdproxy)},
//                owner_{owner}
//        {};
//
//std::shared_ptr<IMDRunnerBuilder> RunnerImplState::builder()
//        {
//            if (mdproxy_ == nullptr)
//            {
//                // usage error?
//                throw gmxapi::Exception();
//            }
//            // When build() is called, the high-level runner handle needs to be updated.
//            std::shared_ptr<RunnerProxy> owner;
//            if (owner_.expired())
//            {
//                owner = std::make_shared<RunnerProxy>();
//                owner_ = owner;
//            }
//            else {
//                owner = owner_.lock();
//            }
//            auto builder = std::make_shared<RunnerImplBuilder>(owner, mdproxy_);
//            return builder;
//        }
//
//gmxapi::Status RunnerImplState::run()
//        {
//            return gmxapi::Status(impl_->run() == 0);
//        }
//
//gmxapi::Status RunnerImplState::run(long int nsteps)
//{
//    return gmxapi::Status(impl_->run(nsteps) == 0);
//}
//
//void RunnerImplState::registerMDBuilder(std::unique_ptr<MDBuilder> builder)
//{
//    // This logic is a mess, but should clean up with the redesign underway...
//    throw gmxapi::Exception();
//}
//
//RunnerImplBuilder::RunnerImplBuilder(std::shared_ptr<RunnerProxy> owner, std::shared_ptr<MDProxy> md) :
//    md_{std::move(md)},
//    owner_{std::move(owner)}
//{}
//
//std::shared_ptr<IMDRunner> RunnerImplBuilder::build()
//{// This is where we would normally provide a registerMDBuilder to md via bind.
//    class RegistrationHelper : public IMDRunner
//    {
//        public:
//            std::string filename;
//            virtual void registerMDBuilder(std::unique_ptr<MDBuilder> builder) override
//            {
//
//                filename = builder->inputAsTprFilename();
//            };
//
//            virtual Status run() override
//            {
//                return Status();
//            }
//    };
//    RegistrationHelper registerer;
//    md_->bind(&registerer);
//    auto filename = registerer.filename;
//    auto newState = std::make_shared<RunnerImplState>(filename);
//    owner_->setState(newState);
//    return newState;}

// Delegate construction to RunnerProxy::RunnerProxy(std::shared_ptr<MDProxy> md)
RunnerProxy::RunnerProxy() :
    RunnerProxy{std::make_shared<MDProxy>()}
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

void EmptyMDRunnerState::registerMDBuilder(std::unique_ptr<MDBuilder> builder)
{
    // nothing to bind to
    throw(Exception());
}

Status EmptyMDRunnerState::run()
{
    // Nothing to run...
    throw(Exception());
    return Status();
}

std::shared_ptr<IMDRunner> EmptyMDRunnerState::initialize(std::shared_ptr<Context> context)
{
    // Runner proxy should have been configured by a builder with something like UninitializedMDRunnerState
    throw(Exception());
    return nullptr;
}

class UninitializedMDRunnerState::Impl
{
    public:
        std::shared_ptr<MDEngine> mdProxy_;
        std::shared_ptr<t_inputrec> inputRecord_;
        std::shared_ptr<t_state> state_;
        std::shared_ptr<gmx_mtop_t> topology_;
};

Status UninitializedMDRunnerState::run()
{
    // We could be more helpful about suggesting the user initialize the runner first...
    throw(Exception());
    return Status();
}


void UninitializedMDRunnerState::registerMDBuilder(std::unique_ptr<MDBuilder> builder)
{

}

std::shared_ptr<IMDRunner> UninitializedMDRunnerState::initialize(std::shared_ptr<Context> context)
{
    RunningMDRunnerState::Builder builder{};


    std::shared_ptr<IMDRunner> initializedRunner = builder.build();
    return initializedRunner;
}

UninitializedMDRunnerState::UninitializedMDRunnerState() :
    impl_{gmx::compat::make_unique<UninitializedMDRunnerState::Impl>()}
{
}

UninitializedMDRunnerState::~UninitializedMDRunnerState() = default;

UninitializedMDRunnerState::Builder::Builder() :
    runner_{}
{
    runner_.reset(new UninitializedMDRunnerState());
}

std::unique_ptr<UninitializedMDRunnerState> UninitializedMDRunnerState::Builder::build()
{
    std::unique_ptr<UninitializedMDRunnerState> runnerState;
    if (runner_->impl_->mdProxy_ && runner_->impl_->inputRecord_ && runner_->impl_->state_ && runner_->impl_->topology_)
    {
        runnerState.swap(runner_);
    }
    else
    {
        throw(Exception());
    }
    return runnerState;
}

UninitializedMDRunnerState::Builder &UninitializedMDRunnerState::Builder::mdEngine(std::shared_ptr<MDEngine> md)
{
    assert(md != nullptr);
    runner_->impl_->mdProxy_ = std::move(md);
    return *this;
}

UninitializedMDRunnerState::Builder &
UninitializedMDRunnerState::Builder::inputRecord(std::shared_ptr<t_inputrec> inputRecord)
{
    assert(inputRecord != nullptr);
    runner_->impl_->inputRecord_ = std::move(inputRecord);
    return *this;
}

UninitializedMDRunnerState::Builder &UninitializedMDRunnerState::Builder::state(std::shared_ptr<t_state> state)
{
    assert(state != nullptr);
    runner_->impl_->state_ = std::move(state);
    return *this;
}

UninitializedMDRunnerState::Builder &UninitializedMDRunnerState::Builder::topology(std::shared_ptr<gmx_mtop_t> topology)
{
    assert(topology != nullptr);
    runner_->impl_->topology_ = std::move(topology);
    return *this;
}


UninitializedMDRunnerState::Builder::~Builder() = default;

class RunningMDRunnerState::Impl
{
    public:
        std::shared_ptr<MDEngine> mdProxy_;
        std::shared_ptr<t_inputrec> inputRecord_;
        std::shared_ptr<t_state> state_;
        std::shared_ptr<gmx_mtop_t> topology_;
        std::shared_ptr<gmx::Mdrunner> runner_;
        decltype(t_inputrec().nsteps) nSteps;
        Impl();
        Status run();
};

RunningMDRunnerState::Impl::Impl()
{
    runner_ = std::make_shared<gmx::Mdrunner>();
    runner_->inputRecord(inputRecord_);
    runner_->stateInput(state_);
    runner_->molecularTopologyInput(topology_);
    // Not implemented
    //runner_->mdEngine(md_);
}

Status RunningMDRunnerState::Impl::run()
{
    Status status{};
    // Todo: check the number of steps to run
    if (runner_->mdrunner())
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
    return impl_->run();
}


std::shared_ptr<IMDRunner> RunningMDRunnerState::initialize(std::shared_ptr<Context> context)
{
    // Should we reinitialize somehow?
    return std::shared_ptr<IMDRunner>();
}

void RunningMDRunnerState::registerMDBuilder(std::unique_ptr<MDBuilder> builder)
{
    // implement the runner--mdengine binding protocol
}


RunningMDRunnerState::Builder::Builder() :
    runner_{new RunningMDRunnerState()}
{
}

std::unique_ptr<RunningMDRunnerState> RunningMDRunnerState::Builder::build()
{
    std::unique_ptr<RunningMDRunnerState> activeRunner = std::move(runner_);
    return activeRunner;
}

RunningMDRunnerState::Builder::~Builder() = default;

} // end namespace gmxapi