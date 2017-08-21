#ifndef GMXAPI_RUNNER_H
#define GMXAPI_RUNNER_H

/*! \file
 * \brief Public C++ API for Gromacs module runners.
 *
 * \ingroup gmxapi
 */

#include <memory>
#include <string>

#include "gmxapi/gmxapi.h"
#include "gmxapi/system.h"

namespace gmxapi
{
class MDProxy;

/// Proxy to an API Gromacs runner.
/*! Can be converted to a concrete runner when work and compute
 * environment are available. Provides a factory (or factory methods or builder(s)) to create an appropriate
 * implementation.
 * \ingroup gmxapi
 */
class RunnerProxy
{
public:
    virtual ~RunnerProxy() = default;
};

class MDBuilder;

/*!
 * \brief Interface provided by runners for MD tasks.
 *
 * A runner implements this interface in order to bind to an MD task. A caller passes an IMDRunner pointer to the
 * bind() method of an MDProxy object. The MDProxy object then provides the runner with a builder for an MD task by
 * calling IMDRunner::registerMDBuilder(std::unique_ptr<ModuleBuilder>).
 *
 * The caller of bind() should guarantee the lifetime of IMDRunner through the subsequent call to registerMDBuilder().
 * The call to bind() should be made when the caller is ready to construct an MD task, such that the state of the
 * MDProxy at the time of the call is appropriate for the MD object to be built. The likely use case is for the call
 * to be made during the builder of a runner that is about to execute.
 *
 * Todo: Check or enforce assumptions: the base class for the builder can guarantee that registerMDBuilder()
 * will be called before destruction.
 *
 * Example:
 * \code
 * class MyRunner : public IMDRunner
 * {
 * private:
 *      std::shared_ptr<MDProxy> md_;
 *      std::unique_ptr<ModuleBuilder> mdbuilder_;
 * public:
 *      IMDRunner* actualRunner() {return this;};
 *      virtual void registerMDBuilder(std::unique_ptr<ModuleBuilder> builder) override
 *      {
 *          mdbuilder_ = std::move(builder);
 *      };
 *
 *      void run()
 *      {
 *          md_->bind(actualRunner);
 *          auto actualMd = mdbuilder_.build();
 *          gmx::MDRunner runner{};
 *          runner.setMD(actualMD);
 *          runner.run();
 *      };
 * };
 *
 * class MyMDProxy : public MDState
 * {
 * public:
 *      std::unique_ptr<ModuleBuilder> builder();
 *      virtual void bind(IMDRunner* runner)
 *      {
 *          runner->registerMDBuilder(builder());
 *      };
 * };
 *
 * \endcode
 */
class IMDRunner
{
    public:
        IMDRunner() = default;
        IMDRunner(const IMDRunner&) = default;
        IMDRunner& operator=(const IMDRunner&) = default;
        virtual ~IMDRunner() = default;
        virtual void registerMDBuilder(std::unique_ptr<MDBuilder> builder) = 0;
        // todo: reconsider interfaces
        virtual gmxapi::Status run() = 0;
        virtual gmxapi::Status run(long int) { throw gmxapi::Exception(); };
};

/*!
 * \brief Get a builder for a concrete runner.
 *
 * A class provides this interface to allow a computational object to be created and launched. The
 * object returned provides a build() method from which to get a runnable object. Other interface
 * features TBD.
 *
 * In general, the class providing this interface will bind to a concrete task before returning from
 * build(). \see registerMDBuiler()
 */
class IRunnerBuilder
{
    public:
        virtual ~IRunnerBuilder() = default;
        /// Build a runner. Return a handle to something that can be run.
        virtual std::shared_ptr<IMDRunner> build() = 0;
};

/// \cond internal
/*! \brief Runner for trivial graphs.
 *
 * The requirements of running only a single Gromacs tool once are substantially
 * different than those of a runner for a chain or data graph. This interface
 * allows classes to offer only a simple implementation.
 *
 * Implementations for object states will depend on execution context and, possibly,
 * on the module to be run.
 * \ingroup gmxapi
 */
class SingleNodeRunnerProxy : public RunnerProxy, public std::enable_shared_from_this<SingleNodeRunnerProxy>
{
    public:
        class State;
        SingleNodeRunnerProxy();
        virtual ~SingleNodeRunnerProxy();

        explicit SingleNodeRunnerProxy(std::shared_ptr<MDProxy> md);

        std::shared_ptr<IRunnerBuilder> builder();

        void setState(std::shared_ptr<State> state);

    private:
        /// bound task, if any
        std::shared_ptr<MDProxy> module_;
        std::shared_ptr<State> state_;
};
/// \endcond



}; // namespace gmxapi


#endif // header guard
