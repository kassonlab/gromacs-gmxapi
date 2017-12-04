#ifndef GMXAPI_MD_IMPL_H
#define GMXAPI_MD_IMPL_H
/*! \file
 * \brief Declarations for molecular dynamics API implementation details.
 *
 * \ingroup gmxapi
 */

#include "gmxapi/md.h"
#include "gmxapi/gmxapi.h"

#include <memory>

namespace gmxapi
{

class MDWorkSpec;

/*!
 * \brief Implementation class to hide guts of MDHolder
 *
 * Holds the gmxapi interface for an object that can help instantiate the gmx::MdRunner
 */
class MDHolder::Impl
{
    public:
        explicit Impl(std::shared_ptr<MDWorkSpec>&& spec);

        std::shared_ptr<MDWorkSpec> spec_{nullptr};
};

} // namespace gmxapi

#endif // header guard
