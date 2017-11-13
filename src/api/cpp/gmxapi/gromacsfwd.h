//
// Created by Eric Irrgang on 11/9/17.
//

#ifndef GMXAPI_GROMACSFWD_H
#define GMXAPI_GROMACSFWD_H

/*! \file
 * \brief Provide forward declarations for symbols in the GROMACS public headers.
 *
 * We don't want to include ::gmx headers if we don't have to, but we need to declare
 * some things in the ::gmx namespace somewhere. These are forward declarations for
 * opaque pointers for client code building against libgmxapi. Client code that is
 * more entwined with libgromacs can include headers from there.
 *
 * It would be nice to include links to the documentation for these classes, too.
 * \ingroup gmxapi
 */

namespace gmx
{

class IRestraintPotential;
class TpxState;

} // end namespace gmx

#endif //GMXAPI_GROMACSFWD_H
