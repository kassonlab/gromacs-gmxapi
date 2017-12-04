/*! \file
 * \brief Implementation details of gmxapi::Context
 *
 */

#include "gmxapi/context.h"

#include <memory>
#include <utility>

#include <cassert>
#include "programs/mdrun/runner.h"
#include "gromacs/mdtypes/tpxstate.h"

#include "gmxapi/gmxapi.h"
#include "gmxapi/session.h"
#include "gmxapi/status.h"

#include "workflow.h"
#include "workflow-impl.h"
#include "session-impl.h"
#include "gromacs/compat/make_unique.h"

namespace gmxapi
{

/*! \brief Temporary shim until proper messaging.
 *
 */
class warn
{
    public:
        /*! \brief Create a warning message.
         *
         * \param message must outlive warn instance...
         */
        explicit warn(const char* message) :
            message_ {message}
        { };
        /// pointer to string managed somewhere else.
        const char* message_;
};

/*!
 * \brief Context implementation base class.
 *
 * Execution contexts have a uniform interface specified by the API. Implementations for
 * particular execution environments can specialize / derive from this base.
 *
 * \todo Separate interface and implementation.
 */
class ContextImpl
{
    public:
        ContextImpl();
        virtual ~ContextImpl() = default;

        /*!
         * \brief Get a reference to the current status object.
         *
         * \return shared ownership of the current status.
         */
        std::shared_ptr<const Status> status() const noexcept;

        /*!
         * \brief Translate the workflow to the execution context and launch.
         *
         * \param work workflow graph
         * \return ownership of a new session
         *
         * \todo This probably makes more sense as a free function.
         */
        std::shared_ptr<Session> launch(std::shared_ptr<ContextImpl> context, const Workflow& work);

        /*!
         * \brief Status of the last operation in the local context.
         *
         * This pointer should always be valid while the context handle exists and
         * client code can extend the life of the object. We use a shared_ptr because
         * it may be expensive or dangerous to copy the object when it is most needed.
         */
        std::shared_ptr<Status> status_;
        std::weak_ptr<Session> session_;
};

ContextImpl::ContextImpl() :
    status_{std::make_shared<Status>(true)},
    session_{}
{
    assert(status_->success());
    assert(session_.expired());
}

std::shared_ptr<const Status> ContextImpl::status() const noexcept
{
    return status_;
}

std::shared_ptr<Session> ContextImpl::launch(std::shared_ptr<ContextImpl> context, const Workflow &work)
{
    // Assume failure until proven otherwise.
    assert(status_ != nullptr);
    *status_ = false;

    std::shared_ptr<Session> session{nullptr};

    // This implementation can only run one workflow at a time.
    if (session_.expired())
    {
        // Check workflow spec, build graph for current context, launch and return new session.
        // \todo This is specific to the session implementation...
        auto mdNode = work.getNode("MD");
        std::string filename{};
        if (mdNode != nullptr)
        {
            filename = mdNode->params();
        }
        auto newMdRunner = gmx::compat::make_unique<gmx::Mdrunner>();
        if (!filename.empty())
        {
            auto tpxState = gmx::TpxState::initializeFromFile(filename);
            newMdRunner->setTpx(std::move(tpxState));
            newMdRunner->initFromAPI();
        }

        {
            auto newSession = SessionImpl::create(std::move(context),
                                                  std::move(newMdRunner));
            session = std::make_shared<Session>(std::move(newSession));
        }


//        for (auto&& node : work)
//        {
//            // If we can do something with the node, do it. If the spec is bad, error.
//        }

//        session = std::make_shared<>();
    }
    // \todo Make some note about the unsuccessful launch.

    if (session != nullptr)
    {
        // Update weak reference.
        session_ = session;
        // Set successful status.
        *status_ = true;
    }
    return session;
}

// In 0.0.3 there is only one Context type
Context::Context() :
    impl_ {gmx::compat::make_unique<ContextImpl>()}
{
    assert(impl_ != nullptr);
}

std::shared_ptr<Session> Context::launch(const Workflow& work)
{
    return impl_->launch(impl_, work);
}

Context::Context(std::shared_ptr<ContextImpl> &&impl) :
    impl_{std::move(impl)}
{
    assert(impl_ != nullptr);
}

Context::~Context() = default;

std::unique_ptr<Context> defaultContext()
{
    auto impl = gmx::compat::make_unique<ContextImpl>();
    auto context = gmx::compat::make_unique<Context>(std::move(impl));
    return context;
}

} // end namespace gmxapi
