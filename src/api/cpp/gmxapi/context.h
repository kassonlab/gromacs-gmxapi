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
class IRunner;
class MDInput;

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

        // Disallow copy
        Context(const Context&) = delete;
        Context& operator=(const Context&) = delete;

        // Allow move
        Context(Context&&) = default;
        Context& operator=(Context&&) = default;

    private:
        class Impl;
        std::unique_ptr<Impl> impl_;
};

std::unique_ptr<Context> defaultContext();

}      // end namespace gmxapi

#endif // header guard
