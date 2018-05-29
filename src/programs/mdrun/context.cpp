//
// Created by Eric Irrgang on 5/16/18.
//

#include "context.h"

#include "gromacs/mdlib/simulationsignal.h"

#include "runner.h"

namespace gmx
{
namespace md
{

Context::Context(const Mdrunner &runner)
{
    runner_ = &runner;
}

SimulationSignals* Context::simulationSignals() const
{
    return runner_->signals();
}


} // end namespace md
} // end namespace gmx
