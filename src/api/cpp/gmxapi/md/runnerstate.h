//
// Created by Eric Irrgang on 8/16/17.
//

#ifndef GROMACS_RUNNERPROXY_H
#define GROMACS_RUNNERPROXY_H

#include "gmxapi/runner.h"
#include <memory>
#include <gromacs/mdtypes/inputrec.h>
#include <gromacs/mdtypes/state.h>
#include <gromacs/topology/topology.h>

namespace gmxapi
{

class EmptyMDRunnerState : public IMDRunner
{
    public:
        EmptyMDRunnerState() = default;

        Status run() override;

        std::shared_ptr<IMDRunner> initialize(std::shared_ptr<Context> context) override;

        void registerMDBuilder(std::unique_ptr<MDBuilder> builder) override;

//        std::shared_ptr<IMDRunnerBuilder> builder() override;
};

/*!
 * \brief An MDRunner that has not yet started.
 *
 * Accumulates configuration information that can be used to launch a gmx::Mdrunner.
 */
class UninitializedMDRunnerState : public IMDRunner
{
    public:
        ~UninitializedMDRunnerState() override;

        // Disallow copy
        UninitializedMDRunnerState(const UninitializedMDRunnerState&) = delete;
        UninitializedMDRunnerState& operator=(const UninitializedMDRunnerState&) = delete;

        // Allow move
        UninitializedMDRunnerState(UninitializedMDRunnerState&&) noexcept = default;
        UninitializedMDRunnerState& operator=(UninitializedMDRunnerState&&) noexcept = default;

        Status run() override;

        std::shared_ptr<IMDRunner> initialize(std::shared_ptr<Context> context) override;

        void registerMDBuilder(std::unique_ptr<MDBuilder> builder) override;

        class Builder
        {
            public:
                Builder();
                ~Builder();
                Builder(const Builder&) = delete;
                Builder(Builder&&) noexcept = default;
                Builder& operator=(const Builder&) = delete;
                Builder& operator=(Builder&&) noexcept = default;

                Builder& mdEngine(std::shared_ptr<MDEngine> md);
                Builder& inputRecord(std::shared_ptr<t_inputrec> inputRecord);
                Builder& state(std::shared_ptr<t_state> state);
                Builder& topology(std::shared_ptr<gmx_mtop_t> topology);
                std::unique_ptr<UninitializedMDRunnerState> build();
            private:
                std::unique_ptr<UninitializedMDRunnerState> runner_;
        };

    private:
        UninitializedMDRunnerState();
        /// Private implementation class
        class Impl;
        /// pointer to implementation
        std::unique_ptr<Impl> impl_;
};

/*!
 * \brief Handle to an active gmx::Mdrunner
 */
class RunningMDRunnerState : public IMDRunner
{
    public:
        ~RunningMDRunnerState() override;

        Status run() override;

        std::shared_ptr<IMDRunner> initialize(std::shared_ptr<Context> context) override;

        void registerMDBuilder(std::unique_ptr<MDBuilder> builder) override;

        // Todo: we can just template some of this: class Builder : public RunnerBuilder<RunningMDRunnerState>
        class Builder
        {
            public:
                Builder();
                ~Builder();
                Builder(const Builder&) = delete;
                Builder(Builder&&) noexcept = default;
                Builder& operator=(const Builder&) = delete;
                Builder& operator=(Builder&&) noexcept = default;

                std::unique_ptr<RunningMDRunnerState> build();
            private:
                std::unique_ptr<RunningMDRunnerState> runner_;
        };
    private:
        RunningMDRunnerState();
        /// Private implementation class
        class Impl;
        /// pointer to implementation
        std::unique_ptr<Impl> impl_;};


} // end namespace gmxapi
#endif //GROMACS_RUNNERPROXY_H
