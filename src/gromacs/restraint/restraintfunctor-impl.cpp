//
// Created by Eric Irrgang on 10/30/17.
//

#include "restraintfunctor-impl.h"

namespace gmx
{

namespace restraint
{


RestraintFunctor::RestraintFunctor(const pull_t &pull,
                                   const t_mdatoms &atoms,
                                   const t_pbc &pbc,
                                   const t_commrec &cr,
                                   double t,
                                   real lambda,
                                   const RVec &x)
{

}

void RestraintFunctor::calculate(real *energy,
                                 RVec *force,
                                 tensor *virial,
                                 real *work)
{

}

} // end namespace gmx::restraint

} // end namespace gmx
