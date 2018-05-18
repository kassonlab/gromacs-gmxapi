//
// Created by Eric Irrgang on 5/18/18.
//



#include <gromacs/compat/make_unique.h>
#include "gmxapi/md/mdsignals.h"

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

        void call() override
        {
            auto signals = runner_->signals();
            signals->at(eglsSTOPCOND).sig = true;
        }

    private:
        gmx::Mdrunner* runner_;
};


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

        std::unique_ptr<Signal::SignalImpl> signalImpl = gmx::compat::make_unique<StopSignal>(runner);

        Signal functor{std::move(signalImpl)};

        return functor;
//    }
//    else
//    {
//    }
}

} // end namespace gmxapi
