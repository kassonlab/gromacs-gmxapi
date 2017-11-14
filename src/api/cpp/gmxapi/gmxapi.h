/*! \file
 * \brief Header for public Gromacs C++ API
 *
 * API clients include this header. It is intended to provide a minimal set of
 * declarations to allow independently implemented API clients to make compatible
 * references to gmxapi objects. Each client will still need to include additional
 * headers and to link against gmxapi, but clients don't need to be completely ABI
 * or API compatible with each other. Clients should use the gmxapi versioning utilities
 * to check for compatibility before accessing members of an object passed by
 * another client.
 *
 * \ingroup gmxapi
 */
/*! \mainpage
 * The `gmxapi` library allows GROMACS extension and client code to interact with
 * GROMACS internals without compile-time dependence on the GROMACS source code.
 * It is sufficient to use the headers described here and to link against `libgmxapi` (or whatever it ends up being called in your installation). CMake helpers make it
 * easy to compile your own code against an installed copy of GROMACS.
 *
 * Using the `Gromacs::gmxapi` CMake target, the headers described here can be
 * included from `gmxapi/...` and symbols in the `::gmxapi` C++ namespace can be
 * resolved in the library by link time.
 *
 * For API stability and unambiguous versioning, headers in the top-level `gmxapi/`
 * directory should not refer to the ::gmx namespace used by the core GROMACS library,
 * except for gromacsfwd.h, which consolidates forward declarations in the ::gmx namespace
 * required by the extension APIs in subdirectories like `gmxapi/md`.
 * \todo gromacsfwd.h probably shouldn't be in the top level either...
 * To reduce dependencies, headers should not depend on headers at a deeper level
 * than themselves, where versioning and compatibility guarantees are weaker
 * (possibly dependent on GROMACS versions).
 *
 * This API is developed primarily to support a Python interface to GROMACS, developed
 * as a completely separate software repository.
 *
 * Refer to the <a href="modules.html">modules</a> section for a hierarchical overview of the API documentation.
 */
/*!
 * \defgroup gmxapi gmxapi
 *
 * \brief Provide external access to an installed GROMACS binary through the gmxapi library.
 *
 * API client code primarily operates on Proxy objects, which can be copied,
 * serialized, and configured independently of execution context, data location,
 * or parallelism.
 *
 * When data needs to be accessed or moved, smart Handle objects are used.
 * DataHandle objects explicitly describe ownership and access characteristics,
 * and so are not copyable, but can have ownership transferred. DataHandle objects
 * should be held briefly, since the implementation can depend on the current
 * execution context and can affect performance.
 *
 * Computation Modules and Runners have Proxy objects whose lifetime is independent
 * of execution context, and Implementation objects whose lifetimes are necessarily
 * bounded by the execution Context in which they run.
 *
 * Proxy and Implementation classes use Interfaces defined by the API, which can also be
 * subclassed in client code to extend the library functionality.
 *
 * The Implementation class supporting a given proxy object is likely to change between
 * successive handle requests, using a "state" behavioral design pattern.
 *
 * API classes participating in behavioral protocols beyond RAII construction and
 * destruction use implementation inheritance or separate manager classes to
 * ensure proper behavior by classes
 * implementing the relevant Interface. E.g. if the API specifies that an
 * Interface requires calls to object->initialize() to be followed by calls to
 * object->deinitialize(), the API library provides mechanisms to guarantee
 * that clients of the interface execute an appropriate state machine.
 *
 * The runner for a workflow cannot be fully instantiated until the execution context
 * is active, but the execution context cannot be configured without some knowledge
 * of the work to be performed. Thus, the generic runner API object participates in initializing
 * the execution context, and a particular context implementation can provide
 * specialized factory extensions to instantiate and initialize the concrete runner
 * implementation, which exists at the library API level and is not directly exposed to
 * the external API.
 *
 * For example:
 *
 * A client may create a MDProxy to pass to a RunnerProxy to pass to a Context, optionally
 * keeping references to the proxies or not.
 * The client asks the Context for a Session.
 * Context asks RunnerProxy for a RunnerBuilder, which asks MDProxy for an MD Engine builder.
 * Each builder (and product) gets a handle to the Context.
 * MDEngineBuilder produces and returns a handle to a new MDEngine.
 * RunnerBuilder produces a Runner with a handle to the MDEngine and returns runner handle to session.
 * Session makes API calls to runner handle until it no longer needs it and releases it.
 * Releasing the handle can provide notification to the runner.
 * When no references remain (as via handles), the runner is destroyed (same for MDEngine).
 * In the case of complex shutdown, the Context holds the last references to the runners and manages shutdown.
 *
 *
 */

#ifndef GMXAPI_H
#define GMXAPI_H

#include <memory>

/*! \brief Contains the external C++ Gromacs API.
 *
 * High-level interfaces for client code is provided in the gmxapi namespace.
 *
 * \ingroup gmxapi
 */
namespace gmxapi
{

// Forward declarations for other gmxapi classes.
class MDEngine;
class Context;
class Status;

// In order to create Python bindings, the class needs to be defined, but its members don't.
class MDHolder {
    public:
        explicit MDHolder(std::shared_ptr<MDEngine> md);

        std::shared_ptr<MDEngine> getMDEngine();
        const std::shared_ptr<MDEngine> getMDEngine() const;
    private:
        std::shared_ptr<MDEngine> md_;
};

}

#endif // header guard
