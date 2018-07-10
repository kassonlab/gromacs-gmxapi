//
// Created by Eric Irrgang on 7/10/18.
//

#ifndef GMXAPI_MDSIGNALS_IMPL_H
#define GMXAPI_MDSIGNALS_IMPL_H

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

class StopSignal : public Signal::SignalImpl
{
    public:
        explicit StopSignal(gmx::Mdrunner* runner);

        void call() override;

    private:
        gmx::Mdrunner* runner_;
};

/*!
 * \brief Manage signal paths exposed through session resources to gmxapi operations.
 *
 * Manages signals for a single gmx::Mdrunner. Currently only supports a stop signal that
 * is required to be issued by all registered possible issuers before the signal is sent to
 * the associated runner. This is not what we want in the long run. This class should handle
 * signal inputs to operations that take signals as input (like Mdrunner) and should allow
 * multiple subscribers. For additional signal processing, such as boolean operations,
 * additional operations should be inserted in a chain.
 */
class SignalManager
{
    public:
        explicit SignalManager(gmx::Mdrunner* runner);
        ~SignalManager();

        /*!
         * \brief Add a name to the list of operations that will be using this signal.
         */
        void addSignaller(std::string name);

        /*!
         * \brief Allow a registered signaller to retrieve a functor.
         *
         * \param name Registered signal issuer.
         * \return Generic Signal object.
         *
         * \throws gmxapi::ProtocolError if named signaller was not previously registered.
         */
        Signal getSignal(std::string name,
                         md::signals signal);

        /*!
         * \brief A member class that can access SignalManager's private members.
         */
        class LogicalAND;

    private:
        gmx::Mdrunner* runner_;
        /*!
         * \brief Track whether the signal has been issued by each registrant.
         */
        std::map<std::string, std::atomic_bool> called_;
};



} //end namespace gmxapi

#endif //GMXAPI_MDSIGNALS_IMPL_H
