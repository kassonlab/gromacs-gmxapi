//
// Created by Eric Irrgang on 10/30/17.
//

#ifndef GROMACS_RESTRAINTFUNCTOR_H
#define GROMACS_RESTRAINTFUNCTOR_H

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/real.h"

/*! \libinternal \file
 * \brief Provide a functor class to clarify calculation input, output, and point in time.
 */

struct pull_t;
struct t_mdatoms;
struct t_pbc;
struct t_commrec;


namespace gmx
{

namespace restraint
{

/*!
 * \brief Implement class for applying the restraint calculations.
 *
 * This class implements the stage of MD integration at which the contributions from restraints
 * (such as for structural refinement or free-energy calculation) are calculated.
 *
 * Example:
 *
 *     void MdIntegrator::applyConstraints()
 *     {
 *          auto calculator = RestraintFunctor(pull_work_, mdatoms_, pbc_, commRec_, time_, lambda_, position_);
 *          calculator.calculate(&energy_, &force_, &virial_, &work_);
 *     };
 */
class RestraintFunctor
{
    public:
        /*!
         * \brief Construct the functor with necessary input parameters.
         * \param pull
         * \param atoms
         * \param pbc
         * \param cr
         * \param t
         * \param lambda
         * \param x
         */
        RestraintFunctor(const pull_t& pull, const t_mdatoms& atoms, const t_pbc& pbc, const t_commrec& cr, double t, real lambda, const RVec& x);

        void calculate(real* energy, RVec* force, tensor* virial, real* work);

    private:

};

} // end namespace gmx::restraint

} // end namespace gmx

#endif //GROMACS_RESTRAINTFUNCTOR_H
