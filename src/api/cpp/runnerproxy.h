//
// Created by Eric Irrgang on 8/16/17.
//

#ifndef GROMACS_RUNNERPROXY_H
#define GROMACS_RUNNERPROXY_H

#include "gmxapi/runner.h"
#include "runner-impl.h"

namespace gmxapi
{

/*!
 * \brief Base class for SingleNodeRunner states
 */
class SingleNodeRunnerProxy::State
{
    public:
        virtual ~State() = default;
        virtual std::shared_ptr<IRunnerBuilder> builder() = 0;
};

class RunnerImplState : public SingleNodeRunnerProxy::State, public IMDRunner
{
    private:
        // This is the soon-to-be-replaced cut-and-paste filename-based mdrunner...
        std::shared_ptr<gmxapi::RunnerImpl> impl_;
        std::shared_ptr<MDProxy> mdproxy_;
        std::weak_ptr<SingleNodeRunnerProxy> owner_;

    public:
        explicit RunnerImplState(std::string filename);
        ~RunnerImplState() override ;
        explicit RunnerImplState(std::shared_ptr<SingleNodeRunnerProxy> owner, std::shared_ptr<MDProxy> mdproxy);
        std::shared_ptr<IRunnerBuilder> builder() override;
        gmxapi::Status run() override;
        gmxapi::Status run(long int nsteps) override;

        virtual void registerMDBuilder(std::unique_ptr<MDBuilder> builder) override;

};


class RunnerImplBuilder : public IRunnerBuilder
{
    private:
        std::shared_ptr<MDProxy> md_;
        std::shared_ptr<SingleNodeRunnerProxy> owner_;
    public:
        virtual ~RunnerImplBuilder() override = default;
        explicit RunnerImplBuilder(std::shared_ptr<SingleNodeRunnerProxy> owner, std::shared_ptr<MDProxy> md);

        virtual std::shared_ptr<IMDRunner> build() override;
};

} // end namespace gmxapi
#endif //GROMACS_RUNNERPROXY_H
