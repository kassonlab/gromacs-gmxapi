//
// Created by Eric Irrgang on 5/18/18.
//

#include "gmxapi/md/mdsignals.h"

#include <atomic>
#include <algorithm>
#include "gmxapi/exceptions.h"

#include "gromacs/compat/make_unique.h"
#include "gromacs/mdlib/simulationsignal.h"
#include "programs/mdrun/runner.h"

#include "gmxapi/session.h"
#include "sessionresources-impl.h"
#include "mdsignals-impl.h"

namespace gmxapi {


Signal::Signal(Signal &&signal) noexcept = default;
Signal &Signal::operator=(Signal &&signal) noexcept = default;

Signal::Signal(std::unique_ptr<gmxapi::Signal::SignalImpl>&& impl) :
    impl_{std::move(impl)}
{
}

Signal::~Signal() = default;

void Signal::operator()()
{
    impl_->call();
}

StopSignal::StopSignal(gmx::Mdrunner *runner) : runner_{runner}
{}

void StopSignal::call()
{
    auto signals = runner_->signals();
    // sig > 0 stops at next NS step. sig < 0 stops at next step.
    signals->at(eglsSTOPCOND).sig = -1;
}

void SignalManager::addSignaller(std::string name)
{
    called_[name].store(false);
}

class SignalManager::LogicalAND : public Signal::SignalImpl
{
    public:
        LogicalAND(SignalManager* manager, std::string name) :
            name_{std::move(name)},
            manager_{manager}
        {}

        void call() override
        {
            auto& callCounter = manager_->called_.at(name_);
            callCounter.store(true);
            using pairType = typename decltype(manager_->called_)::value_type;
            if (std::all_of(manager_->called_.cbegin(),
                            manager_->called_.cend(),
                            [](const pairType& p){ return p.second.load(); }))
            {
                StopSignal(manager_->runner_).call();
            }
        }

    private:
        const std::string name_;
        SignalManager* manager_;
};

Signal SignalManager::getSignal(std::string name,
                                md::signals signal)
{
    if (called_.find(name) == called_.end())
    {
        std::string message = name + " is not registered for this signal.";
        throw gmxapi::ProtocolError(std::move(message));
    }

    if(signal != md::signals::STOP)
    {
        throw gmxapi::NotImplementedError("This signaller only handles stop signals.");
    }

    auto signalImpl = gmx::compat::make_unique<LogicalAND>(this, name);
    Signal functor{std::move(signalImpl)};
    return functor;
}

Signal getMdrunnerSignal(SessionResources *resources,
                         md::signals signal)
{
//// while there is only one choice...
//    if (signal == md::signals::STOP)
//    {
    if(signal != md::signals::STOP)
    {
        throw gmxapi::NotImplementedError("This signaller only handles stop signals.");
    };

    assert(resources);

    auto signaller = resources->getMdrunnerSignal(signal);

    return signaller;
}

} // end namespace gmxapi
