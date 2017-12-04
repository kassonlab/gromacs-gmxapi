#ifndef GMXAPI_MD_H
#define GMXAPI_MD_H
/*! \file
 * \brief Declare base classes and API for MD simulation engines.
 *
 * Helper functions, standard concrete classes, and implementation interfaces are in gmxapi/md/
 * \ingroup gmxapi_md
 */
/*! \defgroup gmxapi_md Molecular Dynamics
 * \brief API access to Molecular Mechanics and Molecular Dynamics calculation in GROMACS
 *
 * \ingroup gmxapi
 */
/*! \dir md
 * \brief Additional declarations for API implementation and extension.
 *
 * This directory contains headers that require gromacsfwd.h and that declare
 * objects that have stronger dependencies on GROMACS to fully define or extend.
 * \ingroup gmxapi_md
 */
#include "exceptions.h"
#include <memory>
#include <string>

namespace gmxapi
{

class MDModule;

/*! \addtogroup gmxapi_md

 # Extending MD with a custom restraint potential

 Below, we show code that extends the GROMACS library, code that interfaces with this
 API, and client code that calls GROMACS with the custom code.
 Refer to the sample plugin repository for a more complete example.

 \todo Move this to an example source code file that we can compile and test.

 ## Example

 In client code, extend the GROMACS library by implementing a new restraint
 potential (see library documentation).

 The gmxapi protocol to register a gmxapi::MDModule with a gmxapi::IMDRunner
 by passing an gmx::IRestraintPotential to a gmx::MdRunner is described in the
 gmxapi::IMDRunner docs. To exercise it, we need to call gmxapi::IMDRunner::setRestraint(),
 passing a std::shared_ptr<gmxapi::MDModule> argument.

        class NullRestraint : public gmx::IRestraintPotential
        {
            public:
                gmx::PotentialPointData evaluate(gmx::Vector r1,
                                                 gmx::Vector r2,
                                                 double t) override
                {
                    return {};
                }
        };

 Use gmxapi::Module to define an API object class that we can pass around.

        class SimpleApiModule : public gmxapi::MDModule
        {
            public:
                const char *name() override
                {
                    return "NullApiModule";
                }

                // Implement the MDModule protocol.
                std::shared_ptr<gmx::IRestraintPotential> getRestraint() override
                {
                    auto restraint = std::make_shared<NullRestraint>();
                    return restraint;
                }
        };

 C++ client code to run an MD simulation with the custom restraint.

        bool mysim() {
                auto system = gmxapi::fromTprFile(filename);
                std::shared_ptr<gmxapi::Context> context = gmxapi::defaultContext();
                auto runner = system->runner();

                auto session = runner->initialize(context);

                auto module = std::make_shared<SimpleApiModule>();
                session->setRestraint(module);

                gmxapi::Status status;
                status = session->run();
                return status.success();
        };

*/

/*!
 * \brief Container for Molecular Dynamics simulation setup.
 *
 * Client code provides the specification for MD work through an object of this type and registers
 * the object in the computing context when an execution session is launched. The contents of the
 * MDWorkSpec are used to pass appropriate parameters to the MD runner.
 *
 * \ingroup gmxapi_md
 */
class MDWorkSpec
{
    public:
        MDWorkSpec();
        ~MDWorkSpec();

        /*!
         * \brief Grant shared ownership of a modular MD computation object
         *
         * \param module instance that can produce a IRestraintPotential at runtime.
         */
        void addModule(std::shared_ptr<gmxapi::MDModule> module);

        /*!
         * \brief Get a handle to the stored list of modules
         *
         * Future versions of MDWorkSpec will not directly hold and grant access to module instances.
         * \return reference that is only valid for the life of this object.
         */
        std::vector<std::shared_ptr<gmxapi::MDModule>>& getModules();
    private:
        /// \cond internal
        /// \brief Private implementation class
        class Impl;
        /// \brief Opaque pointer to implementation object.
        std::unique_ptr<Impl> impl_{nullptr};
        /// \endcond
};

// Forward declaration. Defined in runner.h
class IMDRunner;
// Forware declaration. Defined below.
class MDBuilder;

/// \cond internal
/*!
 * \brief Base class for Molecular Dynamics engine implementations and states.
 *
 * An MD task can have a handle before, during, or after execution, and the local handle may refer
 * to a different implementation class depending on whether execution takes place locally or remotely,
 * and whether distributed data structures are cached locally, etc.
 *
 * Pure virtual base class is not instantiated by clients directly. Objects are created
 * by other API objects or helper functions. \see mdFromTpr()
 *
 * A State instance provides the instantaneous implementation for the owning object.
 *
 * The MDState implementation is responsible for implementing the bind() method, allowing a
 * runner proxy and MD proxy to be translated into an actual runner and MD Engine instance.
 *
 * \ingroup gmxapi_md
 */
class MDEngine
{
    public:
        MDEngine() = default;
        virtual ~MDEngine();
        MDEngine(const MDEngine&) = delete;
        MDEngine& operator=(const MDEngine&) = delete;
        MDEngine(MDEngine&&) = default;
        MDEngine& operator=(MDEngine&&) = default;
        /*!
         * \brief Get a builder for an MD Engine
         *
         * This method allows a caller to convert a proxy object or uninitialized MDEngine into
         * a runnable functor or for unprivileged code to advance the state engine.
         *
         * The default implementation produces an empty proxy, but it might
         * be important to provide a different behavior (or make pure virtual) to help catch usage errors when called on
         * states that do not have clear semantics for builder().
         *
         * \return ownership of a MD Engine builder implementing the gmxapi::MDBuilder interface.
         */
        virtual std::unique_ptr<MDBuilder> builder();

        /// Allow implementing classes to provide information in a generic way.
        virtual const std::string info() const;

        /*
         * \brief Bind to a runner.
         *
         * Implement the runner binding protocol. \See gmxapi::IMDRunner::registerMDBuilder().
         * An object implementing the IMDRunner interface may bind. MDProxy will register a
         * function pointer with which to
         * request a builder for an actual MDEngine.
         */
        virtual // Implemented in libgmxapi. Maybe it should be easier to get the behavioral template?...
        void bind(IMDRunner* runner);
};

/*! \brief Proxy object for an MD engine.
 *
 * Not instantiated by clients directly. Objects are created
 * by other API objects or helper functions. \see mdFromTpr()
 *
 * \note The copy semantics of copying proxy objects are not yet clear, particularly as to whether copied
 * proxies ought to continue to point to the same state member. Since proxy objects are not necessary
 * to manage the underlying resources, other API objects may keep alive a proxy's state member, potentially
 * issuing a new proxy for it at some point.
 *
 * \ingroup gmxapi_md
 */
class MDProxy : public MDEngine
{
    public:
        MDProxy();

        ~MDProxy() override = default;
        MDProxy(MDProxy &) = delete;
        MDProxy &operator=(MDProxy &) = delete;
        MDProxy(MDProxy &&proxy) noexcept;
        MDProxy &operator=(MDProxy &&proxy) noexcept;

        std::unique_ptr<MDBuilder> builder() override;

        /// Note: the caller can retain access to the state argument through whatever interfaces it implements...
        void setState(std::shared_ptr<MDEngine> state);

        /// Get some human-readable status information.
        const std::string info() const override;
    private:
        std::shared_ptr<MDEngine> instanceState_;
};

/*!
 * \brief Build an MD Engine functor at runtime.
 *
 * Using the interface defined by MDBuilder, an implementing class provides a runner with a builder
 * with which to construct an MD engine at run time. The providing object can provide an appropriate
 * concrete builder for the configured task. By the time an MDBuilder reference is returned to the
 * calling code, the MD engine may already be largely configured.
 *
 * For example, an MD task loaded from a TPR file can be provided by a builder that only requires
 * build() to be called to produce a ready-to-run simulation.
 *
 * This builder is part of the MDRunner protocol and the IMDRunner interface.
 *
 * \ingroup gmxapi_md
 */
class MDBuilder
{
    public:
        virtual ~MDBuilder() = default;
        virtual std::unique_ptr<MDEngine> build() = 0;

        // helper functions
        virtual std::string inputAsTprFilename()
        {
            // This is probably a temporary function or at the least would not be generally implemented
            // for a while. We just need it for 0.0.1 right now.
            throw gmxapi::Exception();
        };

        // piece-by-piece construction functions
};

/// \endcond

/*! \brief Get a proxy by reading a TPR file.
 *
 * \param filename TPR input file name.
 * \returns Ownership of a new simulation task proxy object.
 *
 * \ingroup gmxapi_md
 */
std::unique_ptr<MDProxy> mdFromTpr(const std::string filename);


}      // end namespace gmxapi

#endif // header guard
