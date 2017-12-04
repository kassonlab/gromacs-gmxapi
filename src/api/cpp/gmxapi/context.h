#ifndef GMXAPI_CONTEXT_H
#define GMXAPI_CONTEXT_H
/*! \file
 * \brief Declares classes representing execution context.
 *
 * \ingroup gmxapi
 */


#include <memory>

namespace gmxapi
{

class Status;
class Workflow;
class Session;

/*!
 * \brief Context implementation abstract base class.
 *
 * Context Implementations depend on the execution environment, hardware resources, and
 * possibly other factors, and so are not constructed directly but by helper functions.
 * Their details are not exposed at the high level API.
 */
class ContextImpl;

/// Execution context.
/*!
 * The execution context represents computing resources and zero, one, or more
 * workflows to execute. All API objects exist in some context, which determines
 * how the API objects interact underneath the hood.
 *
 * A proxy can be configured with information needed to initialize a runtime
 * environment capable of executing a work load, independently of defining the
 * work.
 * The actual execution
 * environment is not necessarily instantiated / configured until the work is
 * performed.
 * Thus, construction of a Context object does not necessarily imply
 * initialization of compute resources, but any active compute resources are
 * appropriately deinitialized when the object is destroyed. However, to allow
 * opportunities for exception handling, resources should be deinitialized when
 * and as work is completed by the Runner.
 *
 * Ultimately, the class in this header file should just define an interface.
 * Implementations for different execution environments will be provided by the
 * library or extension code and documented separately,
 * and it should be straight-forward to
 * write extension code for other execution environments.
 *
 * In the first draft, it is difficult to separate object definition from
 * initialization.
 * \ingroup gmxapi
 */
class Context
{
    public:
        /*!
         * \brief Get a handle to a new default context object.
         */
        Context();
        ~Context();

        // Nearly trivial copy
        Context(const Context&) = default;
        Context& operator=(const Context&) = default;

        // Allow move
        Context(Context&&) = default;
        Context& operator=(Context&&) = default;

        explicit Context(std::shared_ptr<ContextImpl> &&impl);

        /*!
         * \brief Launch a workflow in the current context, if possible.
         *
         * \param work Configured workflow to instantiate.
         * \return Ownership of a new session or nullptr if not possible.
         *
         * Context maintains a weak reference to the running session and a Status object
         * that can be examined if launch fails due to an invalid work specification or
         * incompatible resources.
         */
        std::shared_ptr<Session> launch(const Workflow& work);

    private:
        /*!
         * \brief Private implementation may be shared by several interfaces.
         */
        std::shared_ptr<ContextImpl> impl_;
};

/*!
 * \brief Construct a context appropriate for the current environment.
 *
 * \return ownership of a new context handle.
 */
std::unique_ptr<Context> defaultContext();

}      // end namespace gmxapi

#endif // header guard
