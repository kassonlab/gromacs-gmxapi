//
// Created by Eric Irrgang on 11/29/17.
//

#include "gmxapi/session.h"

#include <cassert>

#include "gromacs/compat/make_unique.h"
#include "gromacs/mdlib/sighandler.h"
#include "gromacs/restraint/restraintpotential.h"
#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/init.h"

#include "gmxapi/context.h"
#include "gmxapi/exceptions.h"
#include "gmxapi/status.h"
#include "gmxapi/md/mdmodule.h"

#include "mdsignals-impl.h"
#include "session-impl.h"
#include "sessionresources-impl.h"

namespace gmxapi
{

class MpiContextManager
{
    public:
        MpiContextManager()
        {
            gmx::init(nullptr, nullptr);
#ifdef GMX_MPI
#if GMX_MPI
            assert(gmx_mpi_initialized());
#endif
#endif
        };

        ~MpiContextManager()
        {
            gmx::finalize();
        }

        // Disallow copying
        MpiContextManager(const MpiContextManager&) = delete;
        MpiContextManager& operator=(const MpiContextManager&) = delete;

        // Trivial move
        MpiContextManager(MpiContextManager&&) noexcept = default;
        MpiContextManager& operator=(MpiContextManager&&) noexcept = default;
};

/*!
 * \brief Check if a an object can be considered "open".
 *
 * This should be generalized to an API idiom.
 *
 * \tparam T type that can be open or closed.
 * \param object something that has a concept of "open" or "closed."
 * \return true if open, false if closed, compiler error if non-sensical.
 */
template<class T>
bool isOpen(const T& object);
//{
//    (void) object;
//    static_assert(false, "Compiler could not find open/close concept for the given object.");
//    return false;
//}

template<>
bool isOpen<SessionImpl>(const SessionImpl& object)
{
    return object.isOpen();
}

template<>
bool isOpen<Session>(const Session& object)
{
    return object.isOpen();
}


SignalManager::SignalManager(gmx::Mdrunner *runner) : runner_(runner)
{

}

SignalManager::~SignalManager() = default;

bool SessionImpl::isOpen() const noexcept
{
    return status_ != nullptr;
}

Status SessionImpl::status() const noexcept
{
    return *status_;
}

std::unique_ptr<Status> SessionImpl::close()
{
    // When the Session is closed, we need to know that the MD output has been finalized, which currently requires
    // gmx::MDrunner::~MDrunner() to be called.
    runner_.reset();
    std::unique_ptr<Status> status{nullptr};
    status.swap(status_);
    assert(status_ == nullptr);
    return status;
}

Status SessionImpl::run() noexcept
{
    // Status is failure until proven otherwise.
    Status status{false};
    assert(runner_ != nullptr);
    auto rc = runner_->mdrunner();
    if (rc == 0)
    {
        status = true;
    }
    return status;
}

std::unique_ptr<SessionImpl> SessionImpl::create(std::shared_ptr<ContextImpl> context,
                                                 std::unique_ptr<gmx::Mdrunner> runner)
{
    std::unique_ptr<SessionImpl> impl{new SessionImpl(std::move(context), std::move(runner))};
    return impl;
}

SessionImpl::SessionImpl(std::shared_ptr<ContextImpl> context,
                         std::unique_ptr<gmx::Mdrunner> runner) :
    status_{gmx::compat::make_unique<Status>(true)},
    context_{std::make_shared<Context>(std::move(context))},
    mpiContextManager_{gmx::compat::make_unique<MpiContextManager>()},
    runner_{std::move(runner)},
    signal_{gmx::compat::make_unique<SignalManager>(runner_.get())}
{
    assert(status_ != nullptr);
    assert(context_ != nullptr);
    assert(mpiContextManager_ != nullptr);
    assert(runner_ != nullptr);
    assert(signal_ != nullptr);
    // For the libgromacs context, a session should explicitly reset global variables that could
    // have been set in a previous simulation during the same process.
    gmx_sighandler_reset();
}

Status SessionImpl::setRestraint(std::shared_ptr<gmxapi::MDModule> module)
{
    assert(runner_ != nullptr);
    Status status{false};

    if (module != nullptr)
    {
        auto restraint = module->getRestraint();
        if (restraint != nullptr)
        {
            auto sessionResources = createResources(module);
            if (!sessionResources)
            {
                status = false;
            }
            else
            {
                runner_->addPullPotential(restraint, module->name());
                status = true;
            }
        }
    }
    return status;
}

gmx::Mdrunner *SessionImpl::getRunner()
{
    gmx::Mdrunner * runner{nullptr};
    if (runner_)
    {
        runner = runner_.get();
    }
    return runner;
}

gmxapi::SessionResources *SessionImpl::getResources(const std::string &name) const noexcept
{
    gmxapi::SessionResources * resources{nullptr};
    try
    {
        resources = resources_.at(name).get();
    }
    catch (const std::out_of_range& e)
    {
        // named operation does not have any resources registered.
    };

    return resources;
}

gmxapi::SessionResources *SessionImpl::createResources(std::shared_ptr<gmxapi::MDModule> module) noexcept
{
    // check if resources already exist for this module
    // If not, create resources and return handle.
    // Return nullptr for any failure.
    gmxapi::SessionResources * resources{nullptr};
    if (resources_.find(module->name()) == resources_.end())
    {
        auto resourcesInstance = gmx::compat::make_unique<SessionResources>(this, module->name());
        resources_.emplace(std::make_pair(module->name(), std::move(resourcesInstance)));
        resources = resources_.at(module->name()).get();
        // This should be more dynamic.
        getSignalManager()->addSignaller(module->name());
        auto restraint = module->getRestraint();
        if (restraint)
        {
            restraint->bindSession(resources);
        }
    };
    return resources;
}

SignalManager *SessionImpl::getSignalManager()
{
    SignalManager* ptr{nullptr};
    if (isOpen())
    {
        ptr = signal_.get();
    }
    return ptr;
}

Session::Session(std::unique_ptr<SessionImpl>&& impl) noexcept :
    impl_{std::move(impl)}
{
    assert(impl_ != nullptr);
    assert(impl == nullptr);
    assert(impl_->isOpen());
}

Status Session::run() noexcept
{
    assert(impl_ != nullptr);

    Status status{impl_->run()};
    return status;
}

Status Session::close()
{
    assert(impl_ != nullptr);

    Status status{false};
    if (isOpen())
    {
        // \todo catch exceptions when we know what they might be
        auto status_ptr = impl_->close();
        if (status_ptr != nullptr)
        {
            status = *status_ptr;
        }
        // what to do if we get nullptr?
    }

    return status;
}

Session::~Session()
{
    assert(impl_ != nullptr);
    if (isOpen())
    {
        try
        {
            impl_->close();
        }
        catch (const std::exception&)
        {
            // \todo find some exception-safe things to do with this via the Context interface.
        }
    }

}

bool Session::isOpen() const noexcept
{
    assert(impl_ != nullptr);
    const auto result = impl_->isOpen();
    return result;
}

Status setSessionRestraint(Session *session,
                           std::shared_ptr<gmxapi::MDModule> module)
{
    auto status = gmxapi::Status(false);

    if (session != nullptr && module != nullptr)
    {
        auto sessionImpl = session->getRaw();

        assert(sessionImpl);
        status = sessionImpl->setRestraint(std::move(module));
    }
    return status;
}

SessionImpl *Session::getRaw() const noexcept
{
    return impl_.get();
}

std::shared_ptr<Session> launchSession(Context* context, const Workflow& work) noexcept
{
    auto session = context->launch(work);
    return session;
}

SessionImpl::~SessionImpl() = default;

SessionResources::SessionResources(gmxapi::SessionImpl *session,
                                   std::string name) :
    sessionImpl_{session},
    name_{std::move(name)}
{
}

SessionResources::~SessionResources() = default;

const std::string SessionResources::name() const
{
    return name_;
}

Signal SessionResources::getMdrunnerSignal(md::signals signal)
{
    //// while there is only one choice...
//    if (signal == md::signals::STOP)
//    {
    if(signal != md::signals::STOP)
    {
        throw gmxapi::NotImplementedError("This signaller only handles stop signals.");
    };

    // Get a signalling proxy for the caller.
    auto signalManager = sessionImpl_->getSignalManager();
    if(signalManager == nullptr)
    {
        throw gmxapi::ProtocolError("Client requested access to a signaller that is not available.");
    };
    auto functor = signalManager->getSignal(name_, signal);

    return functor;
}

} // end namespace gmxapi
