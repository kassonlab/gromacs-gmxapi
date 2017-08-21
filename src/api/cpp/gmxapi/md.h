#ifndef GMXAPI_MD_H
#define GMXAPI_MD_H

/*! \file
 * \brief Declare Proxy and API objects for MD simulation engines.
 * \ingroup gmxapi
 */
#include "exceptions.h"
#include <memory>
#include <string>

namespace gmxapi
{
class IMDRunner;
class MDState;

/*! \brief Proxy object for an MD engine.
 *
 * Not instantiated by clients directly. Objects are created
 * by other API objects or helper functions. \see mdFromTpr()
 *
 * \note The copy semantics of copying proxy objects are not yet clear, particularly as to whether copied
 * proxies ought to continue to point to the same state member. Since proxy objects are not necessary
 * to manage the underlying resources, other API objects may keep alive a proxy's state member, potentially
 * issuing a new proxy for it at some point.
 */
class MDProxy
{
    public:
        MDProxy();

        ~MDProxy();

        MDProxy(const MDProxy &proxy);

        MDProxy(MDProxy &&proxy) noexcept;

        MDProxy &operator=(const MDProxy &proxy);

        MDProxy &operator=(MDProxy &&proxy) noexcept;

        /*!
         * \brief Bind to a runner.
         *
         * An object implementing the IMDRunner interface may bind. MDProxy will register a function pointer with which to
         * request a builder for an actual MDEngine.
         */
        void bind(IMDRunner* runner);

        /// \endcond


        /// Note: the caller can retain access to the state argument through whatever interfaces it implements...
        void setState(std::shared_ptr<MDState> state);

        /// Get some human-readable status information.
        std::string info();
    private:
        std::shared_ptr<MDState> state_;
};

class MDEngine
{
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

/*! \brief Get a proxy by reading a TPR file.
 *
 * \param filename TPR input file name.
 * \returns Ownership of a new simulation task proxy object.
 */
std::unique_ptr<MDProxy> mdFromTpr(const std::string filename);

}      // end namespace gmxapi

#endif // header guard
