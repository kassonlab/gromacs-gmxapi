//
// Created by Eric Irrgang on 11/13/17.
//

#include "gmxapi/md/mdmodule.h"

namespace gmxapi
{

MDModule::~MDModule() = default;

std::shared_ptr<::gmx::IRestraintPotential> MDModule::getRestraint()
{
    return nullptr;
}

} // end namespace gmxapi
