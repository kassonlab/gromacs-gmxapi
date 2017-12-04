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
//
//class Atoms;
//class MDProxy;
//class IMDRunner;
//class MDEngine;

// Forward declaration for a return type defined elsewhere.
class Session;

/// Container for molecular model and simulation parameters.
/*!
 * \cond
 * A system instance is sort of a container of builders, and a Context is sort of a factory. together they allow a simulation to be constructed and initialized with the appropriate implementations of runner, integrator, and data objects.
 *
 * Ultimately, the system needs to pass serialized data sufficient to reconstruct
 * itself as part of the workflow it contains when the work is launched.
 * \endcond
 * \ingroup gmxapi
 */
class System final
{
    public:
        /*! \brief Private implementation class
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
        // Note there is confusing overlap in the use of these two functions that should be normalized.
        std::shared_ptr<MDWorkSpec> getSpec();

        /*!
         * \brief Configure the computing environment for the specified workflow.
         *
         * \return Ownership of a ready-to-run workflow or nullptr if there were errors.
         *
         * If errors occur, they will be stored in the context object. If run without
         * and argument, launch() uses the current context of the System object. If a
         * context argument is given, the system and its configured workflow are
         * translated to the provided context and launched.
         *
         * \todo Policy: does System then track the (potentially remote) context or should
         * it be considered to have "forked", and the new session object retrieved from
         * the session handle if needed.
         */
        std::shared_ptr<Session> launch();
        std::shared_ptr<Session> launch(std::shared_ptr<Context> context);

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
