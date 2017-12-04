//
// Created by Eric Irrgang on 11/29/17.
//

#include "gmxapi/session.h"

#include <cassert>
#include "gromacs/compat/make_unique.h"

#include "gmxapi/context.h"
#include "gmxapi/status.h"

#include "session-impl.h"

namespace gmxapi
{

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

bool SessionImpl::isOpen() const noexcept
{
    return status_ == nullptr;
}

Status SessionImpl::status() const noexcept
{
    return *status_;
}

std::unique_ptr<Status> SessionImpl::close()
{
    std::unique_ptr<Status> status{nullptr};
    status.swap(status_);
    return status;
}

Status SessionImpl::run() noexcept
{
    // Status is failure until proven otherwise.
    Status status{false};
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
    status_{gmx::compat::make_unique<Status>()},
    context_{std::make_shared<Context>(std::move(context))},
    runner_{std::move(runner)}
{
    assert(status_ != nullptr);
    assert(context_ != nullptr);
    assert(runner_ != nullptr);
}

Session::Session(std::unique_ptr<SessionImpl>&& impl) noexcept :
    impl_{std::move(impl)}
{
    assert(impl_ != nullptr);
    assert(impl == nullptr);
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

std::shared_ptr<Session> launchSession(Context* context, const Workflow& work) noexcept
{
    auto session = context->launch(work);
    return session;
}

} // end namespace gmxapi
