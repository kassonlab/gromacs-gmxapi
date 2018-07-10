//
// Created by Eric Irrgang on 12/1/17.
//

#ifndef GMXAPI_SESSION_IMPL_H
#define GMXAPI_SESSION_IMPL_H
/*! \file
 * \brief Declare implementation interface for Session API class(es).
 *
 * \ingroup gmxapi
 */

#include <map>
#include "programs/mdrun/runner.h"

#include "gmxapi/session/resources.h"

namespace gmxapi
{

// Forward declaration
class MpiContextManager; // Locally defined in session.cpp
class ContextImpl; // locally defined in context.cpp
class SignalManager; // defined in mdsignals-impl.h

/*!
 * \brief Implementation class for executing sessions.
 *
 * In 0.0.3, there is only one context and only one session type, but this will likely change soon.
 * \ingroup gmxapi
 */
class SessionImpl
{
    public:
        /// Use create() factory to get an object.
        SessionImpl() = delete;
        ~SessionImpl();

        /*!
         * \brief Check if the session is (still) running.
         *
         * When a session is launched, it should be returned in an "open" state by the launcher function.
         * \return True if running, false if already closed.
         */
        bool isOpen() const noexcept;

        /*!
         * \brief Get the current / most recent status.
         *
         * \return Copy of current status object reflecting most recent operation.
         */
        Status status() const noexcept;

        /*!
         * \brief Explicitly close the session.
         *
         * Sessions should be explicitly `close()`ed to allow for exceptions to be caught by the client
         * and because closing a session involves a more significant state change in the program than
         * implied by a typical destructor. If close() can be shown to be exception safe, this protocol may be removed.
         *
         * On closing a session, the status object is transfered to the caller.
         * \return
         */
        std::unique_ptr<Status> close();

        /*!
         * \brief Run the configured workflow to completion or error.
         *
         * \return copy of the resulting status.
         *
         * \internal
         * By the time we get to the run() we shouldn't have any unanticipated exceptions.
         */
        Status run() noexcept;

        static std::unique_ptr<SessionImpl> create(std::shared_ptr<ContextImpl> context,
                                                   std::unique_ptr<gmx::Mdrunner> runner);

        Status setRestraint(std::shared_ptr<gmxapi::MDModule> module);

        /*! \internal
         * \brief API implementation function to retrieve the current runner.
         *
         * \return non-owning pointer to the current runner or nullptr if none.
         */
        gmx::Mdrunner* getRunner();

        /*!
         * \brief Get a handle to the resources for the named session operation.
         *
         * \param name unique name of element in workflow
         * \return temporary access to the resources.
         *
         * If called on a non-const Session, creates the resource if it does not yet exist. If called on a const Session,
         * returns nullptr if the resource does not exist.
         */
        gmxapi::SessionResources* getResources(const std::string& name) const noexcept;

        gmxapi::SessionResources* createResources(std::shared_ptr<gmxapi::MDModule> module) noexcept;

        SignalManager* getSignalManager();

    private:
        /*!
         * \brief Private constructor for use by create()
         *
         * \param context specific context to keep alive during session.
         * \param runner ownership of live Mdrunner object.
         */
        SessionImpl(std::shared_ptr<ContextImpl> context,
                    std::unique_ptr<gmx::Mdrunner> runner);

        /*!
         * \brief Manage session resources for named workflow elements.
         */
        std::map<std::string, std::unique_ptr<SessionResources>> resources_;

        /*!
         * \brief Current / most recent Status for the session.
         *
         * \internal
         * An open session has a valid status object. A closed session has status_ == nullptr.
         */
        std::unique_ptr<Status> status_;

        /*!
         * \brief Extend the life of the owning context. The session will get handles for logging, UI status messages, and other facilities through this interface.
         */
        std::shared_ptr<Context> context_;

        /*!
         * \brief RAII management of gmx::init() and gmx::finalize()
         */
        std::unique_ptr<MpiContextManager> mpiContextManager_;

        std::unique_ptr<gmx::Mdrunner> runner_;

        std::unique_ptr<SignalManager> signal_;
};

} //end namespace gmxapi

#endif //GMXAPI_SESSION_IMPL_H
