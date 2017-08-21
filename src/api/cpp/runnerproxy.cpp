//
// Created by Eric Irrgang on 7/31/17.
//

#include "gmxapi/runner.h"
#include "runnerproxy.h"
#include "runner-impl.h"

#include <string>

#include "gmxapi/md.h"

namespace gmxapi
{

/*!
 * \brief State implementation for first draft hack.
 *
 * Uses an outdated cut-and-pasted mdrunner with minor encapsulation and a filename parameter.
 */
RunnerImplState::RunnerImplState(std::string filename) :
                impl_{std::make_shared<gmxapi::RunnerImpl>(filename)},
                mdproxy_{gmxapi::mdFromTpr(filename)}
        {};

RunnerImplState::~RunnerImplState()
        {
        }

RunnerImplState::RunnerImplState(std::shared_ptr<SingleNodeRunnerProxy> owner, std::shared_ptr<MDProxy> mdproxy) :
                impl_{nullptr},
                mdproxy_{std::move(mdproxy)},
                owner_{owner}
        {};

std::shared_ptr<IRunnerBuilder> RunnerImplState::builder()
        {
            if (mdproxy_ == nullptr)
            {
                // usage error?
                throw gmxapi::Exception();
            }
            // When build() is called, the high-level runner handle needs to be updated.
            std::shared_ptr<SingleNodeRunnerProxy> owner;
            if (owner_.expired())
            {
                owner = std::make_shared<SingleNodeRunnerProxy>();
                owner_ = owner;
            }
            else {
                owner = owner_.lock();
            }
            auto builder = std::make_shared<RunnerImplBuilder>(owner, mdproxy_);
            return builder;
        }

gmxapi::Status RunnerImplState::run()
        {
            return gmxapi::Status(impl_->run() == 0);
        }

gmxapi::Status RunnerImplState::run(long int nsteps)
{
    return gmxapi::Status(impl_->run(nsteps) == 0);
}

void RunnerImplState::registerMDBuilder(std::unique_ptr<MDBuilder> builder)
{
    // This logic is a mess, but should clean up with the redesign underway...
    throw gmxapi::Exception();
}

RunnerImplBuilder::RunnerImplBuilder(std::shared_ptr<SingleNodeRunnerProxy> owner, std::shared_ptr<MDProxy> md) :
    md_{std::move(md)},
    owner_{std::move(owner)}
{}

std::shared_ptr<IMDRunner> RunnerImplBuilder::build()
{// This is where we would normally provide a registerMDBuilder to md via bind.
    class RegistrationHelper : public IMDRunner
    {
        public:
            std::string filename;
            virtual void registerMDBuilder(std::unique_ptr<MDBuilder> builder) override
            {

                filename = builder->inputAsTprFilename();
            };

            virtual Status run() override
            {
                return Status();
            }
    };
    RegistrationHelper registerer;
    md_->bind(&registerer);
    auto filename = registerer.filename;
    auto newState = std::make_shared<RunnerImplState>(filename);
    owner_->setState(newState);
    return newState;}

SingleNodeRunnerProxy::SingleNodeRunnerProxy() :
    SingleNodeRunnerProxy{std::make_shared<MDProxy>()}
{};

SingleNodeRunnerProxy::SingleNodeRunnerProxy(std::shared_ptr<MDProxy> md) :
    module_{std::move(md)},
    state_{nullptr}
{
}

std::shared_ptr<IRunnerBuilder> SingleNodeRunnerProxy::builder()
{
    // First draft only uses runner implementation initialized from a filename...
    state_ = std::make_shared<RunnerImplState>(this->shared_from_this(), module_);
    return state_->builder();
}

void SingleNodeRunnerProxy::setState(std::shared_ptr<State> state)
{
    state_ = std::move(state);
}

SingleNodeRunnerProxy::~SingleNodeRunnerProxy() = default;


} // end namespace gmxapi