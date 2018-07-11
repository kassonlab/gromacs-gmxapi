//
// Created by Eric Irrgang on 11/29/17.
//

#ifndef GROMACS_SESSION_H
#define GROMACS_SESSION_H

/*! \file
 * \brief Declarations for a workflow execution session and helper functions.
 *
 * \ingroup gmxapi
 */

#include <memory>

namespace gmxapi
{

// forward declaration
class MDModule;
class Status;
class Context;
class Workflow;

/*!
 * \brief Private implementation class for a session.
 *
 * Actual implementation class may depend on the execution context, but this should be irrelevant to
 * a client. The implementation details are not exposed in the high-level API, but may be in the extension
 * API. See developer documentation for details.
 */
class SessionImpl;

/*!
 * \brief Workflow execution session.
 *
 * When a workflow is launched in an execution context, the result is a Session object
 * that serves as a handle to interact with the running workflow. The handle allows dynamic changes to
 * the workflow, or control or data crossing the API boundary.
 *
 * Separating run() from construction allows the client to examine the running execution
 * environment or to retrieve the communicator before beginning long-running computation.
 *
 * The session should be explicitly `close()`ed before destruction to allow exception
 * handling during shutdown. The destructor will close() the session if it is not
 * closed already, but errors may be hidden.
 *
 * The session owns a Status object that can be shared and retained by client code
 * for further inspection.
 *
 * \ingroup gmxapi
 */
class Session
{
    public:
        /// A session must be created by launching a workflow in an execution context.
        Session() = delete;

        /*! \brief A session cannot be copied, only moved.
         *
         * For shared ownership of a session, use a shared pointer.
         */
        Session(const Session&) = delete;
        Session& operator=(const Session&) = delete;

        /*!
         * \brief Pass ownership of a Session.
         */
        Session(Session&&) noexcept = default;
        Session& operator=(Session&&) noexcept = default;

        /*!
         * \brief Construct by taking ownership of an implementation object.
         *
         * \param impl Concrete object to take ownership of.
         */
        explicit Session(std::unique_ptr<SessionImpl>&& impl) noexcept;

        /*!
         * \brief Destroy Session.
         *
         * If the session is still active (has not been closed) then it will be closed
         * with exceptions suppressed. If possible, problems encountered will be
         * noted in the Status object, which the client may have retained shared
         * ownership of.
         */
        ~Session();

        /*!
         * \brief Close a running session.
         *
         * close() should be called before destroying the Session object so that the
         * client can catch any exceptions thrown during shut down that may be
         * unavoidable in the parallel computing environment.
         *
         * \return status of close() operation.
         */
        Status close();

        /*!
         * \brief Run the current workflow to completion.
         *
         * The client should examine the Status object for errors and resulting workflow state.
         * \return the Session's Status after the run.
         */
        Status run() noexcept;

        bool isOpen() const noexcept;

        /*! \internal
         * \brief Get a non-owning handle to the implementation object.
         *
         * Get a raw pointer to the implementation object. The pointer is valid only during the lifetime of the Session,
         * so retain a shared pointer to this Session object or only hold the pointer for the duration of a code block
         * guaranteed to exist entirely within the lifetime of a Session object.
         *
         * \return opaque pointer used by gmxapi implementation and extension code.
         */
        SessionImpl* getRaw() const noexcept;

    private:
        /// \brief opaque pointer to implementation
        std::unique_ptr<SessionImpl> impl_;
};

/*!
 * \brief Set a uniquely identifiable restraint instance on the MD simulator.
 *
 * \param session
 * \param module
 * \return
 */
Status setSessionRestraint(Session *session,
                           std::shared_ptr<gmxapi::MDModule> module);

/*!
 * \brief Launch a workflow in the provided execution context.
 *
 * \param context Execution environment
 * \param work Directed acyclic graph defining workflow and data flow.
 * \return non-unique ownership of the session or nullptr in failure.
 *
 * The provided context maintains a weak reference to the executing session, while
 * the session extends the life of the context.
 *
 * If session fails to launch (nullptr returned) examine the Status of the context
 * for details.
 */
std::shared_ptr<Session> launchSession(Context* context, const Workflow& work) noexcept;


} // end namespace gmxapi

#endif //GROMACS_SESSION_H