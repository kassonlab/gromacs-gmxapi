//
// Created by Eric Irrgang on 7/10/18.
//

#ifndef GMXAPI_SESSION_RESOURCES_H
#define GMXAPI_SESSION_RESOURCES_H

/*! \file
 * \brief Define interface to Session Resources for running gmxapi operations.
 */

namespace gmxapi {

/*!
 * \brief Handle to Session-provided resources.
 *
 * Session handle for workflow elements requiring resources provided through the Session.
 *
 * Provided during launch through gmx::IRestraintPotential::bindSession()
 *
 * No public interface yet. Use accompanying free functions.
 * \see gmxapi::getMdrunnerSignal()
 */
class SessionResources;

} // end namespace gmxapi

#endif //GMXAPI_SESSION_RESOURCES_H
