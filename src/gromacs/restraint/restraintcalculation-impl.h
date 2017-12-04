//
// Created by Eric Irrgang on 10/30/17.
//

#ifndef GROMACS_RESTRAINTCALCULATION_IMPL_H
#define GROMACS_RESTRAINTCALCULATION_IMPL_H


#include <memory>
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/real.h"
#include "restraintcalculation.h"

struct t_mdatoms;
struct t_pbc;


namespace gmx
{
namespace restraint
{

/*!
 * \brief State class for a set of restraints calculations.
 */
class Calculation : public ICalculation
{
    public:
        Calculation() = delete;
        Calculation(double time,
                            const t_commrec &commRec,
                            const t_mdatoms &atoms,
                            const t_pbc &pbc,
                            real lambda,
                            const rvec *positions,
                            pull_t *puller,
                            rvec *forces,
                            tensor virial);

        /*!
         * \brief Get energy calculated for current time.
         *
         * \implements ICalculation
         * \return kJ / mole
         */
        real energy() const noexcept override;

        /*!
         * \brief Get generalized work for free energy calculation.
         *
         * \implements ICalculation
         * \return dV/dlambda
         */
        real work() const noexcept override ;

        /*!
         * \brief Get current time.
         *
         * \return Simulation time value at which calculation was performed.
         */
        double time() const noexcept;
    private:

        double time_;
        real energy_;
        real work_;
};

}
}


#endif //GROMACS_RESTRAINTCALCULATION_IMPL_H
