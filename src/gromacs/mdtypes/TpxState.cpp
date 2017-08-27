//
// Created by Eric Irrgang on 8/27/17.
//

#include <string>

#include "gromacs/compat/make_unique.h"

#include "gromacs/fileio/tpxio.h"
#include "gromacs/utility/smalloc.h"
#include "TpxState.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/topology/topology.h"

namespace gmx
{

TpxState::TpxState() :
    inputrecInstance_{std::make_shared<t_inputrec>()},
    stateInstance_{std::make_shared<t_state>()},
    mtop_{nullptr},
    initialized_{false},
    dirty_{false}
{
    snew(mtop_, 1);
}

TpxState::~TpxState()
{
    if (mtop_ != nullptr)
    {
        sfree(mtop_);
    }
};

std::unique_ptr<TpxState> TpxState::initializeFromFile(const char* filename)
{
    std::string arg;
    arg = filename;
    return initializeFromFile(arg);
}

std::unique_ptr<TpxState> TpxState::initializeFromFile(const std::string &filename)
{
    auto newState = gmx::compat::make_unique<TpxState>();
    read_tpx_state(filename.c_str(), newState->inputrecInstance_.get(), newState->stateInstance_.get(), newState->mtop_);

    newState->initialized_ = true;
    return newState;
}

t_inputrec *TpxState::getRawInputrec()
{
    dirty_ = true;
    return inputrecInstance_.get();
}

gmx_mtop_t *TpxState::getRawMtop()
{
    dirty_ = true;
    return mtop_;
}

t_state *TpxState::getRawState()
{
    dirty_ = true;
    return stateInstance_.get();
}

bool TpxState::isInitialized() const
{
    return initialized_;
};



} // end namespace gmx