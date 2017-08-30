#include "atoms.h"

#include "gromacs/compat/make_unique.h"
#include "gromacs/math/paddedvector.h"
#include "gromacs/mdtypes/state.h"

namespace gmxapi
{

template <class T>
bool hasVelocity(const T &)
{
    return false;
};

/// Whether state structure provides velocity.
template<>
bool hasVelocity(const t_state &state)
{
    return (state.flags & estV);
};

Atoms::Atoms(const t_state &state)
{
    using gmx::PaddedRVecVector;
    // What is the correct way to identify the size of PaddedRVecVector
    x_ = std::make_shared<PaddedRVecVector>(state.natoms + 1);
    if (x_ != nullptr)
    {
        std::copy(state.x.begin(), state.x.end(), x_->begin());
    }
    if (hasVelocity(state))
    {
        v_ = std::make_shared<PaddedRVecVector>(state.natoms + 1);
        if (v_ != nullptr)
        {
            std::copy(state.v.begin(), state.v.end(), v_->begin());
        }
    }
}

std::shared_ptr<PaddedRVecVector> Atoms::x()
{
    return x_;
}

std::unique_ptr<Atoms> Atoms::handle()
{
    return gmx::compat::make_unique<Atoms>(*this);
}


} // end namespace gmxapi
