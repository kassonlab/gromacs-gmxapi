#ifndef GMXAPI_DATA_H
#define GMXAPI_DATA_H

/*! \file
 * \brief Classes of molecular data defined by the API.
 *
 * Other API classes may be specified as providers or consumers of these classes.
 * Objects of these types are proxy objects. If raw data must be accessed, an
 * appropriate Handle must be obtained, possibly forcing network data transfers.
 * ``const`` instances of proxy objects may only produce read-only data handles.
 * Writable handles obtain exclusive ownership and lock the data until the
 * handle is released. If a non-shared writeable copy of data is desired,
 * consider exportCopy() instead. A read-only handle may be obtained from a
 * non-const proxy by explicitly specifying the access level template parameter.
 *
 * Other data specifications can be reqested via template parameters as well.
 *
 * Example:
 *
 *     // create an empty (proxy) data object
 *     Positions r;
 *
 *     // Provide initialization data. Behavior may depend on execution context.
 *     r.initialize({{0,0,0},{0,1,0},{0,-1,0}});
 *
 *     // Get array-like access to data initialized with the contents of r at
 *     // the time the function is called, using the default floating point
 *     // precision with which the library was built. Other precision or data
 *     // layout can be requested with explicit template parameters.
 *     auto my_positions = r.exportCopy();
 *
 *     r.exportHandle();
 *
 *     try
 *     {
 *         // Lock the data structure
 *         writehandle = r.exportHandle();
 *     }
 *     catch (const gmxapi::NoAccess& e)
 *     {
 *         // Could not obtain lock on resource.
 *     }
 *
 */

namespace gmxapi
{

/*! \brief 3D spatial vectors.
 */
class Positions
{
    public:
        /*! \brief Make full numeric data available to the caller.
         *
         * Writable handles obtain exclusive ownership and lock the data until the
         * handle is released. If a non-shared writeable copy of data is desired,
         * consider exportCopy() instead.
         */
        template<class Scalar>
        const  exportHandle();
};


}

#endif //header guard
