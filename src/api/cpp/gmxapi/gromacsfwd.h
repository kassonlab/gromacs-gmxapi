//
// Created by Eric Irrgang on 11/9/17.
//

#ifndef GMXAPI_GROMACSFWD_H
#define GMXAPI_GROMACSFWD_H

/*! \ingroup gmxapi
 * \file
 * \brief Provide forward declarations for symbols in the GROMACS public headers.
 *
 * Basic API clients only need to compile
 * and link against the gmxapi target, but some gmxapi classes use opaque pointers to
 * library classes that are forward-declared here.
 *
 * We don't want to include ::gmx headers if we don't have to, but we need to declare
 * some things in the ::gmx namespace somewhere. These are forward declarations for
 * opaque pointers in libgromacs for client code building against libgmxapi.
 * Client code that is
 * more entwined with libgromacs can include headers from there.
 *
 * For maximal compatibility with other libgmxapi clients (such as third-party
 * Python modules), client code should use the wrappers and protocols in the
 * gmxapi.h header. Note that there is a separate CMake target to build the full
 * developer documentation for gmxapi.
 *
 * Refer to GMXAPI developer docs for the protocols that map gmxapi interfaces to
 * GROMACS library interfaces.
 * Refer to the GROMACS developer
 * documentation for details on library interfaces forward-declared in this header.
 *
 * \todo It would be nice to include links to the documentation for these classes, too.
 */

namespace gmx
{

class IRestraintPotential;
class TpxState;

} // end namespace gmx

#endif //GMXAPI_GROMACSFWD_H
