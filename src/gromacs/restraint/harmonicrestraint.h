//
// Created by Eric Irrgang on 10/2/17.
//

#ifndef GMX_PULLING_ROUXRESTRAINT_H
#define GMX_PULLING_ROUXRESTRAINT_H

#include "restraintpotential.h"

namespace gmx
{

/*!
 * \brief Refine and apply a pair distance restraint.
 *
 * Applies a restraint to selected pair separation vectors to bias sampled configurations
 * towards an externally-defined distribution.
 *
 * The potential is a functional of selected pair distances
 * evaluated
 */
class HarmonicRestraint : public IRestraintPotential
{

};

} // end namespace gmx

#endif //GMX_PULLING_ROUXRESTRAINT_H
