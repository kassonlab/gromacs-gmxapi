//
// Created by Eric Irrgang on 5/18/18.
//

/* WARNING
 * This whole file is not intended to make it into a public release and is not part of the gmxapi API. It is for
 * prototyping only. Please don't let it slip into a release without serious design considerations.
 */

#ifndef GMXAPI_MDSIGNALS_H
#define GMXAPI_MDSIGNALS_H

#include <memory>

namespace gmxapi {

namespace md
{

enum class signals {
        STOP
};

} // end namespace md


class Session;  // defined in gmxapi/session.h

/*!
 * \brief Proxy for signalling function objects.
 *
 * Objects of this type are simple callables that issue a specific signal.
 */
class Signal
{
    public:
        class SignalImpl;
        explicit Signal(std::unique_ptr<SignalImpl>&& signal);
        Signal(Signal&& signal);
        Signal& operator=(Signal&& signal);
        ~Signal();

        void operator()();

    private:
        std::unique_ptr<SignalImpl> impl_;
};

/*!
 * \brief Get a function object that issues a signal to the currently active MD runner.
 *
 * \param session pointer to the active Session.
 * \return Callable function object handle
 */
Signal getMdrunnerSignal(Session* session, md::signals signal);

} // end namespace md

#endif //GMXAPI_MDSIGNALS_H
