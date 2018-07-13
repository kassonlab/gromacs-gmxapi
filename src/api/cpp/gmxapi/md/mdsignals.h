//
// Created by Eric Irrgang on 5/18/18.
//

#ifndef GMXAPI_MDSIGNALS_H
#define GMXAPI_MDSIGNALS_H

/*! \internal
 * \file
 * \brief Temporary infrastructure for signalling MD simulations.
 *
 * These interfaces are not considered to be part of the gmxapi spec, but will exist in 0.0.6 and
 * possibly 0.0.7 until more abstract data flow is available to MD plugin developers, at which point
 * any remaining functionality here will be moved to private implementation details.
 *
 * \ingroup gmxapi_md
 */

#include <memory>

namespace gmxapi {

/*!
 * \brief Internal details of gmxapi MD functionality
 */
namespace md
{

/*!
 * \brief Symbolic signal slots for MD signalling.
 *
 * \see getMdrunnerSignal()
 */
enum class signals {
        STOP
};

} // end namespace md


class SessionResources;  // reference gmxapi/session/resources.h

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
        Signal(Signal&&) noexcept ;
        Signal& operator=(Signal&&) noexcept ;
        ~Signal();

        void operator()();

    private:
        std::unique_ptr<SignalImpl> impl_;
};

/*!
 * \brief Get a function object that issues a signal to the currently active MD runner.
 *
 * \param resources pointer to the active Session resources.
 * \param signal type of signal the client would like to issue.
 * \return Callable function object handle
 *
 * \throws gmxapi::NotImplementedError for unknown values of signal.
 */
Signal getMdrunnerSignal(SessionResources *resources,
                         md::signals signal);

} // end namespace md

#endif //GMXAPI_MDSIGNALS_H
