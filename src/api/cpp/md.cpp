#include <memory>
#include <cassert>
#include <iostream>

#include "gmxapi/gmxapi.h"
#include "gmxapi/md.h"
#include "md-impl.h"
#include "gmxapi/md/mdmodule.h"

#include "gromacs/compat/make_unique.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/utility/keyvaluetree.h"

namespace gmxapi
{

class MDWorkSpec::Impl
{
    public:
        static std::unique_ptr<Impl> create();

        std::vector<std::shared_ptr<gmxapi::MDModule>> modules{};
};

std::unique_ptr<MDWorkSpec::Impl> MDWorkSpec::Impl::create()
{
    auto newImpl = gmx::compat::make_unique<MDWorkSpec::Impl>();
    assert(newImpl != nullptr);
    assert(newImpl->modules.empty());
    return newImpl;
}

MDWorkSpec::MDWorkSpec() :
    impl_{Impl::create()}
{
    assert(impl_ != nullptr);
}

void MDWorkSpec::addModule(std::shared_ptr<gmxapi::MDModule> module)
{
    assert(impl_ != nullptr);
    std::cout << "Adding module " << module->name() << " to work specification" << std::endl;
    impl_->modules.emplace_back(std::move(module));
}

std::vector<std::shared_ptr<gmxapi::MDModule>> &MDWorkSpec::getModules()
{
    assert(impl_ != nullptr);
    return impl_->modules;
}

MDWorkSpec::~MDWorkSpec() = default;

std::shared_ptr<::gmxapi::MDWorkSpec> MDHolder::getSpec()
{
    assert(impl_ != nullptr);
    assert(impl_->spec_ != nullptr);
    return impl_->spec_;
}

std::shared_ptr<const ::gmxapi::MDWorkSpec> MDHolder::getSpec() const
{
    assert(impl_ != nullptr);
    assert(impl_->spec_ != nullptr);
    return impl_->spec_;
}

MDHolder::MDHolder() :
    MDHolder{std::make_shared<MDWorkSpec>()}
{
    assert(impl_ != nullptr);
    assert(impl_->spec_ != nullptr);
}

MDHolder::MDHolder(std::shared_ptr<MDWorkSpec> spec) :
    name_{},
    impl_{std::make_shared<MDHolder::Impl>(std::move(spec))}
{
    assert(impl_ != nullptr);
    assert(impl_->spec_ != nullptr);
}

MDHolder::Impl::Impl(std::shared_ptr<MDWorkSpec>&& spec) :
    spec_{spec}
{
    assert(spec_ != nullptr);
}

MDHolder::MDHolder(std::string name) :
    MDHolder{}
{
    name_ = name;
    assert(impl_ != nullptr);
    assert(impl_->spec_ != nullptr);
}

std::string MDHolder::name() const
{
    return name_;
}


} //end namespace gmxapi
