#ifndef GMXAPI_MODULES_H
#define GMXAPI_MODULES_H
/// \cond

/*! \file
 * \brief Public C++ API for Gromacs computational modules.
 *
 * \ingroup gmxapi
 */

#include <memory>

#include "gmxapi/gmxapi.h"

namespace gmxapi
{

/*! \brief Specifies the interface to runners of trivial graphs.
 *
 * The requirements of running only a single Gromacs tool once are substantially
 * different than those of a runner for a chain or data graph. This interface
 * allows classes to offer only a simple implementation.
 *
 * Implementations for this interface depend on execution context and, possibly,
 * on the module to be run.
 */
class MDEngine
{
    public:
        class Factory;
        Status run();
    private:
        class Impl;
        std::unique_ptr<Impl> impl_;
};

}; // namespace gmxapi

/// \endcond
#endif // header guard
