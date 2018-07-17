//
// Created by Eric Irrgang on 7/10/18.
//

#ifndef GMXAPI_MDSIGNALS_IMPL_H
#define GMXAPI_MDSIGNALS_IMPL_H

/*! \file
 * \brief Implementation details for gmxapi::Signal and related gmxapi::md infrastructure.
 *
 * \ingroup gmxapi_md
 */

#include "gmxapi/md/mdsignals.h"

#include <atomic>

#include "gromacs/compat/make_unique.h"
#include "gromacs/mdlib/simulationsignal.h"
#include "programs/mdrun/runner.h"

#include "gmxapi/session.h"

#include "session-impl.h"

namespace gmxapi {


/*!
 * \brief The Signal Implementation interface.
 *
 * A SignalImpl concrete class must implement a `call()` method that issues the signal.
 */
class Signal::SignalImpl
{
    public:
        virtual void call() = 0;
};

/*!
 * \brief Signal implementation for MD simulation stop signals.
 *
 * Provides a call() operator that sets the stop condition for the MD simulation.
 *
 * Client code is not expected to create objects of this type directly, but to retrieve
 * one, wrapped in a gmxapi::Signal for immediate use, from a SignalManager with SignalManager::getSignal()
 *
 * It is the responsibility of the client code to make sure that the address of the Mdrunner
 * remains valid for the lifetime of this object.
 */
class StopSignal : public Signal::SignalImpl
{
    public:
        /*!
         * \brief Create short-lived signal implementation.
         *
         * \param runner non-owning handle
         *
         * The object is constructed with a handle to the runner associated with the SignalManager and
         * owned by the owner of the SignalManager.
         */
        explicit StopSignal(gmx::Mdrunner* runner);

        /*!
         * \brief Set a stop condition for the attached runner.
         */
        void call() override;

    private:
        /// non-owning handle to a runner owned by the owner of the SignalManager.
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
 *
 * SignalManager objects are created during Session launch and are owned exclusively by session
 * implementation objects. If Session::isOpen() is true, the SignalManager should still be valid,
 * but the intended use case is for SignalManager handles to be retrieved immediately before use
 * by implementation code within the library with SessionImpl::getSignalManager().
 *
 * A SignalManager should be created for each consumer (each gmx::Mdrunner) in a Session.
 * This occurs in the SessionImpl::create() function.
 *
 * \ingroup gmxapi_md
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
         * \param signal type of signal the client would like to issue.
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
        /// Non-owning handle to the associated runner.
        gmx::Mdrunner* runner_;
        /*!
         * \brief Track whether the signal has been issued by each registrant.
         */
        std::map<std::string, std::atomic_bool> called_;
};



} //end namespace gmxapi

#endif //GMXAPI_MDSIGNALS_IMPL_H
