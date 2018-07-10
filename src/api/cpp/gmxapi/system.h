#ifndef GMXAPI_SYSTEM_H
#define GMXAPI_SYSTEM_H
/*! \file
 * \brief Declare container for molecular systems
 *
 * \ingroup gmxapi
 */
#include <memory>

#include "gmxapi/gmxapi.h"

namespace gmxapi
{

// Forward declaration for a return type defined elsewhere.
class Session;

/// Container for molecular model and simulation parameters.
/*!
 * \brief Deprecated: A wrapper for gmx::Mdrunner
 *
 * It was not intended to end up as such, but gmxapi::System has ended up as basically a
 * wrapper for the gmx::Mdrunner simulator object. As it is, this class does not fit the
 * current gmxapi paradigm and will be removed, reworked, or renamed soon.
 *
 * # Protocol
 *
 * As of gmxapi 0.0.6, a simulation is configured and launched as follows.
 *
 * 1. Caller gets a System handle with gmxapi::fromTprFile().
 * 2. Caller optionally attaches additional MD Modules with getSpec()->addModule(std::shared_ptr<gmxapi::MDModule> module). See gmxapi::MDHolder
 * 3. Caller gets a runnable object by passing a Context to System::launch()
 *
 * During launch() configured gmxapi::MDModules are attached to the simulator, which is
 * then run by calling run() on the object returned by launch().
 *
 * \ingroup gmxapi
 */
class System final
{
    public:
        /*! \brief Private implementation class.
         *
         * System::Impl does not have a public interface and is only exposed in opaque pointers.
         */
        class Impl;

        /// A blank system object is possible, but not yet useful.
        System();
        /// No copy.
        /*! The semantics of copying a System are ambiguous, so disallow implicit
         * copy. Some sort of prototype or clone idiom is probably useful, but
         * needs to explicitly identify any expensive operations.
         */
        System(const System &)            = delete;
        /// No copy.
        System &operator=(const System &) = delete;

        /// Allow move.
        System(System &&) noexcept;
        /// Allow move.
        System &operator=(System &&) noexcept;

        /*!
         * \brief Create by taking ownership of an implementation object.
         *
         * \param implementation
         */
        explicit System(std::unique_ptr<Impl>&& implementation);

        /// \cond internal
        /// Destructor defined later to allow unique_ptr members of partially-defined types.
        ~System();
        /// \endcond

        Status setRestraint(std::shared_ptr<gmxapi::MDModule> module);

        /*!
         * \brief Borrow shared ownership of the System's container of associated modules.
         *
         * Used with gmxapi::MDHolder to add MD Modules to the simulation to be run.
         *
         * \return handle to be passed to gmxapi::MDHolder
         *
         * \todo This is used in gmxpy but not tested here.
         */
        std::shared_ptr<MDWorkSpec> getSpec();

        /*!
         * \brief Configure the computing environment for the specified workflow.
         *
         * \return Ownership of a ready-to-run workflow or nullptr if there were errors.
         *
         * If errors occur, they will be stored in the context object. If run without
         * an argument, launch() uses the current context of the System object. If a
         * context argument is given, the system and its configured workflow are
         * translated to the provided context and launched.
         *
         * \param context (optional) execution context in which to launch.
         *
         * \note The Session object does not "own" the Context, but must be able to extend the lifetime of the Context in
         * which it is running.
         *
         * \todo Policy: does System then track the (potentially remote) context or should
         * it be considered to have "forked", and the new session object retrieved from
         * the session handle if needed.
         *
         * \cond internal
         * # Protocol
         *
         * The current implementation of System::launch() performs the following actions.
         *
         * When launch() is called, a new gmxapi::Session is created by passing a gmxapi::Workflow to context->launch(). The Workflow basically just contains the TPR filename.
         * 1. A new Mdrunner is created from the information in the gmxapi::Workflow
         * 2. A new Session is created using the ContextImpl and the runner
         *
         * Then, for each module available through getSpec()->getModules(), the session and module
         * are passed to gmxapi::setSessionRestraint().
         * 1. A gmx::IRestraintPotential is retrieved from the module.
         * 2. A unique, named SessionResources is created for the module and attached to the SessionImpl.
         * 3. The SessionResources is passed to IRestraintPotential::bindSession(). Currently, the only thing
         *    the restraint could do at this point is to save a copy of the pointer and later pass it to
         *    gmxapi::getMdrunnerSignal().
         * 4. The restraint is passed to gmx::Mdrunner::addPullPotential(), which adds the restraint to
         *    the global gmx::restraint::Manager, which then needs to be `clear()`ed after the runner completes.
         *
         * Shared ownership of the Session is returned to the caller of launch().
         *
         * \todo gmx::restraint::Manager should not be a singleton.
         * It currently performs a similar role to the gmxapi::System::spec_ member and should be passed to
         * the runner, but we don't have a way to do that right now.
         *
         * \endcond
         *
         * \{
         */
        std::shared_ptr<Session> launch();
        std::shared_ptr<Session> launch(std::shared_ptr<Context> context);
        /// \}

        /*!
         * \brief Get the status of the last API call involving this system.
         *
         * \return copy of the most recent status.
         */
        Status status();

//        /// Get a handle to system atoms.
//        std::unique_ptr<Atoms> atoms();

    private:
        /*!
         * \brief Opaque pointer to implementation.
         */
        std::unique_ptr<Impl> impl_;
};


/// Defines an MD workflow from a TPR file.
/*! The TPR file has sufficient information to fully specify an MD run, though
 * various parameters are implicit until the work is launched. The TPR filename
 * provided must refer to identical TPR files at the API client and at the
 * master rank of the execution host.
 *
 * \param filename Filesystem path of TPR file.
 * \returns gmxapi::System object with the specified workflow.
 * \ingroup gmxapi
 */
std::unique_ptr<gmxapi::System> fromTprFile(std::string filename);

}      // end namespace gmxapi

#endif // include guard
