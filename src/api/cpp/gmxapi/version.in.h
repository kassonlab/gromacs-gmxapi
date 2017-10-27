#ifndef GMXAPI_VERSION_H
#define GMXAPI_VERSION_H
/*! \file
 * \brief Implement versioning API for C++ external Gromacs interface.
 *  Versioning follows semantic versioning scheme in which the major version
 *  specifies the API compatibility level, and minor version indicates additional
 *  features that may or may not be ABI compatible.
 *
 *  Defines a class Version in the gmxapi namespace providing static methods
 *  for use by API client code at compile time and run time.
 *
 * Goals: do better than compiler error, build time linker errors, and run time linker errors.
 *
 *  Todo: provide versioning for headers and library so clients can do abstract comparison of build versus runtime.
 *
 *  Todo: Provide better boilerplate (at least) for client self-test of API version compatibility at build time.
 *
 *  Todo: Add compatibility test/warning at module load when client was compiled against a different libgmxapi.
 *  \ingroup gmxapi
 */

#include <string>

namespace gmxapi
{

// Todo: it may be preferable for CMake to get the version from the header instead of the other way around.
// It would be nice to be able to pull the headers straight from the repository...
static constexpr unsigned int GMXAPI_MAJOR   = @GMXAPI_MAJOR@;
static constexpr unsigned int GMXAPI_MINOR   = @GMXAPI_MINOR@;
static constexpr unsigned int GMXAPI_PATCH   = @GMXAPI_PATCH@;
static const char GMXAPI_RELEASE[] = "@GMXAPI_RELEASE@";

/*!
 * \brief Provide API library version information for client code.
 *
 * Allow client code to query the currently loaded gmxapi library object to find the built version. Provide helpers
 * to compare against the features for which the client was written and the headers against which it was compiled.
 *
 * \ingroup gmxapi
 */
class Version
{
    public:
        /// Query gmxapi major version.
        /// \returns major version number
        static unsigned int major();
        /// Query gmxapi minor version.
        /// \returns minor version number
        static unsigned int minor();
        /// Query gmxapi patch level.
        /// \returns patch level number
        static unsigned int patch();
        /// Get formatted release string.
        /// Format is major.minor.patch
        /// \returns release as string
        static std::string release();
        /// Check features availability
        /// \returns `true` if the named feature is available.
        /// Features introduced after release 1.0.0 may be named in the documentation
        /// to improve readability of client code and simplify development. Prefer
        /// this mechanism when checking for features still under development or to
        /// distinguish between interface levels of a specific feature.
        /// \param featurename Feature name described in the feature's documentation.
        static bool has_feature(std::string featurename);
        /// Check for sufficiently high API version number.
        /// \returns `true` if gmxapi library version is the same or greater than the argument(s).
        /// \param major gmxapi major version number.
        /// \param minor gmxapi minor version number (optional).
        /// \param patch patch level of the api library (optional).
        static bool is_at_least(unsigned int major, unsigned int minor = 0, unsigned int patch = 0);
};

}      // namespace gmxapi

#endif // version.h include guard
