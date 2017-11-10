//
// Created by Eric Irrgang on 11/9/17.
//

#ifndef GMXAPI_MDMODULE_H
#define GMXAPI_MDMODULE_H

#include "gmxapi/gromacsfwd.h"

namespace gmxapi
{



/*!
 * \brief Base class for computational components of MD containers.
 *
 *
 */
class MDModule
{
    public:
        virtual ~MDModule() = default;
        virtual const char* name() { return "MDModule"; };

        virtual std::shared_ptr<::gmx::IRestraintPotential> getRestraint()
        {
            return nullptr;
        };
};

} // end namespace gmxapi

#endif //GMXAPI_MDMODULE_H
