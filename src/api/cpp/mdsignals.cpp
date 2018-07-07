//
// Created by Eric Irrgang on 5/18/18.
//

#include "gmxapi/md/mdsignals.h"

#include <atomic>

#include "gromacs/compat/make_unique.h"
#include "gromacs/mdlib/simulationsignal.h"
#include "programs/mdrun/runner.h"

#include "gmxapi/session.h"

#include "session-impl.h"

namespace gmxapi {

class Signal::SignalImpl
{
    public:
        virtual void call() = 0;


};

Signal::Signal(std::unique_ptr<gmxapi::Signal::SignalImpl>&& impl) :
    impl_{std::move(impl)}
{
}

Signal::~Signal() = default;

void Signal::operator()()
{
    impl_->call();
}

Signal::Signal(Signal &&signal) = default;

Signal &Signal::operator=(Signal &&signal) = default;

class StopSignal : public Signal::SignalImpl
{
    public:
        explicit StopSignal(gmx::Mdrunner* runner) : runner_{runner} {};

        StopSignal(gmx::Mdrunner* runner, unsigned int numParticipants) : StopSignal(runner)
        {
            StopSignal::numParticipants_.store(numParticipants);
            StopSignal::numCalls_.store(0);
        }

        void call() override
        {
            unsigned int n{++StopSignal::numCalls_};
            if (n >= StopSignal::numParticipants_.load())
            {   
                auto signals = runner_->signals();
                // sig > 0 stops at next NS step. sig < 0 stops at next step.
                signals->at(eglsSTOPCOND).sig = -1;
            }
        }

    private:
        gmx::Mdrunner* runner_;

        // Number of participants in this signal
        static std::atomic<unsigned int> numParticipants_;

        // Number of times the signal has been called.
        static std::atomic<unsigned int> numCalls_;
};

std::atomic<unsigned int> StopSignal::numParticipants_{0};
std::atomic<unsigned int> StopSignal::numCalls_{0};

Signal getMdrunnerSignal(Session* session, md::signals signal)
{
//// while there is only one choice...
//    if (signal == md::signals::STOP)
//    {
    assert(signal == md::signals::STOP);
    assert(session);

    auto impl = session->getRaw();
    assert(impl);

    auto runner = impl->getRunner();
    assert(runner);

 //   std::unique_ptr<Signal::SignalImpl> signalImpl = gmx::compat::make_unique<StopSignal>(runner);
    std::unique_ptr<Signal::SignalImpl> signalImpl = gmx::compat::make_unique<StopSignal>(runner, impl->numRestraints);

    Signal functor{std::move(signalImpl)};

    return functor;
//    }
//    else
//    {
//    }
}

} // end namespace gmxapi
