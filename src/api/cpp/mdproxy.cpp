//
// Created by Eric Irrgang on 7/31/17.
//

#include <memory>
#include <functional>

#include "gmxapi/gmxapi.h"
#include "gmxapi/md.h"

#include "gromacs/compat/make_unique.h"
#include "md-impl.h"
#include "runnerproxy.h"

namespace gmxapi
{

MDState::~MDState() = default;

std::string MDState::info()
{
    return "Generic MDState object";
}

std::unique_ptr<MDBuilder> MDState::builder()
{
//    std::unique_ptr<MDBuilder> builder = std::make_unique<MDEngineBuilder>();
//    return builder;
    return nullptr;
}

MDStateFromMDInput::MDStateFromMDInput() : input_{nullptr} {}

MDStateFromMDInput::MDStateFromMDInput(std::unique_ptr<MDInput> input) :
        MDStateFromMDInput{std::move(input), std::string()}
{}

MDStateFromMDInput::MDStateFromMDInput(std::unique_ptr<MDInput> input, std::string metadata) :
    input_{std::move(input)},
    metadata_{metadata}
{}

MDStateFromMDInput::~MDStateFromMDInput() = default;


MDProxy::MDProxy() = default;
MDProxy::~MDProxy() = default;
MDProxy::MDProxy(MDProxy &&proxy) noexcept = default;
MDProxy::MDProxy(const MDProxy &proxy) = default;
MDProxy& MDProxy::operator=(const MDProxy&) = default;
MDProxy& MDProxy::operator=(MDProxy&&) noexcept = default;

void MDProxy::setState(std::shared_ptr<MDState> state)
{
    state_ = state;
}

std::string MDProxy::info()
{
    return state_->info();
}

std::unique_ptr<MDBuilder> MDStateFromMDInput::builder()
{
    return nullptr;
}

std::string MDStateFromMDInput::info()
{
    if (input_ == nullptr)
    {
        return "uninitialized MDStateFromMDInput";
    }
    else
    {
        std::string output("MDStateFromMDInput initialized");
        if (!metadata_.empty())
        {
            output.append(" with metadata: ");
            output.append(metadata_);
        }
        return output;
    }
}

void MDProxy::bind(IMDRunner* runner)
{
    auto builder = state_->builder();
    runner->registerMDBuilder(std::move(builder));
}

} // end namespace gmxapi