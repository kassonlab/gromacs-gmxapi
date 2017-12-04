//
// Created by Eric Irrgang on 10/30/17.
//

#ifndef GROMACS_RESTRAINTCALCULATION_H
#define GROMACS_RESTRAINTCALCULATION_H

#include "gromacs/utility/real.h"

namespace gmx
{

namespace restraint
{

class ICalculation
{
    public:
        virtual ~ICalculation();

        /*!
         * \brief Get energy contribution for the given time.
         *
         * \return Energy contribution at time t in kJ / mole
         */
        virtual real energy() const noexcept  = 0;

        /*!
         * \brief Get change in energy per free energy coordinate lambda value.
         *
         * \return dVdl change in potential energy V per change in lambda value l
         * \ref gmx_enerdata_t
         */
        virtual real work() const noexcept  = 0;

};

}

}

#endif //GROMACS_RESTRAINTCALCULATION_H
