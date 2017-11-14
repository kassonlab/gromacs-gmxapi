//
// Created by Eric Irrgang on 11/9/17.
//

#ifndef GMXAPI_MDMODULE_H
#define GMXAPI_MDMODULE_H
/*! \file
 * \brief Declare MDModule class for interacting with GROMACS library.
 *
 * \ingroup gmxapi_md
 */

#include <memory>

#include "gmxapi/gromacsfwd.h"



namespace gmxapi
{



/*!
 * \brief Base class for computational components of MD containers.
 *
 * Instances of this class and its descendents provide member functions that return
 * GROMACS library objects that are not defined in this API, but which should be
 * defined in the public part of the GROMACS library API. Forward declarations of the library
 * classes are in gmxapi/gromacsfwd.h so that basic API clients only need to compile
 * and link against the gmxapi target, but to extend the API requires GROMACS library
 * headers and possibly linking against `libgromacs`. Refer to the GROMACS developer
 * documentation for details.
 *
 * \ingroup gmxapi_md
 */
class MDModule
{
    public:
        virtual ~MDModule();
        virtual const char* name() { return "MDModule"; };

        /*!
         * \brief Allows module to provide a restraint implementation.
         *
         * To implement a restraint, override this function.
         * \return shared ownership of a restraint implementation or nullptr if not implemented.
         *
         * With future maturation, this interface will presumably be revised to something
         * more abstract, though I'm not sure what form that would take. We will probably
         * still need to have a general set of possible module types defined with the API,
         * in which case it does make sense to have clearly typed dispatching, and
         * `bool hasRestraint = module->getRestraint() != nullptr;` might the simplest thing.
         *
         * Implementing a restraint is explained in the GROMACS developer documentation,
         * which is currently built separately from the GMXAPI documentation.
         * Also, refer to the sample plugin in a repository hosted in the same
         * place this git repository is found.
         */
        virtual std::shared_ptr<::gmx::IRestraintPotential> getRestraint();
};


} // end namespace gmxapi

#endif //GMXAPI_MDMODULE_H
