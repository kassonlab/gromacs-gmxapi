#ifndef GMXAPI_ATOMS_H
#define GMXAPI_ATOMS_H

#include <memory>

#include "gromacs/math/paddedvector.h"

namespace gmx
{
// PaddedRVecVector is an alias, not a class...
//     class PaddedRVecVector;
}

class t_state;

namespace gmxapi
{

/// Handle to a configuration of atoms.
class Atoms
{
    public:
        explicit Atoms(const Atoms &)   = default;
        Atoms &operator=(const Atoms &) = default;

        /// Get a new handle to this Atoms data
        std::unique_ptr<Atoms> handle();
        std::shared_ptr<gmx::PaddedRVecVector> x();

        explicit Atoms(const t_state &state);
    private:
        std::shared_ptr<gmx::PaddedRVecVector> x_;
        std::shared_ptr<gmx::PaddedRVecVector> v_;
};

}      //end namespace gmxapi

#endif // header guard
