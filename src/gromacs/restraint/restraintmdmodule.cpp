//
// Created by Eric Irrgang on 11/10/17.
//

#include <gromacs/mdtypes/iforceprovider.h>
#include "gromacs/compat/make_unique.h"
#include "restraintmdmodule.h"
#include "restraintmdmodule-impl.h"

gmx::RestraintMDModuleImpl::~RestraintMDModuleImpl() = default;

gmx::RestraintMDModuleImpl::RestraintMDModuleImpl(std::shared_ptr<gmx::IRestraintPotential> restraint) :
    forceProvider_{::gmx::compat::make_unique<RestraintForceProvider>(restraint)},
    outputProvider_{::gmx::compat::make_unique<RestraintOutputProvider>()},
    optionProvider_{::gmx::compat::make_unique<RestraintOptionProvider>()}
{

}

gmx::IMdpOptionProvider *gmx::RestraintMDModuleImpl::mdpOptionProvider()
{
    assert(optionProvider_ != nullptr);
    return optionProvider_.get();
}

gmx::IMDOutputProvider *gmx::RestraintMDModuleImpl::outputProvider()
{
    assert(outputProvider_ != nullptr);
    return outputProvider_.get();
}

void gmx::RestraintMDModuleImpl::initForceProviders(ForceProviders *forceProviders)
{
    assert(forceProvider_ != nullptr);
    assert(forceProviders != nullptr);
    forceProviders->addForceProvider(forceProvider_.get());
}


// Needs to be defined after implementation type is complete in order to have unique_ptr member.
gmx::RestraintMDModule::~RestraintMDModule() = default;



gmx::IMdpOptionProvider *gmx::RestraintMDModule::mdpOptionProvider()
{
    assert(impl_ != nullptr);
    return impl_->mdpOptionProvider();
}

gmx::IMDOutputProvider *gmx::RestraintMDModule::outputProvider()
{
    assert(impl_ != nullptr);
    return impl_->outputProvider();
}

void gmx::RestraintMDModule::initForceProviders(ForceProviders *forceProviders)
{
    assert(impl_ != nullptr);
    impl_->initForceProviders(forceProviders);
}

std::unique_ptr<gmx::RestraintMDModule>
gmx::RestraintMDModule::create(std::shared_ptr<gmx::IRestraintPotential> restraint)
{
    auto implementation = ::gmx::compat::make_unique<RestraintMDModuleImpl>(std::move(restraint));
    std::unique_ptr<gmx::RestraintMDModule> newModule{new RestraintMDModule(std::move(implementation))};
    return newModule;
}

// private constructor to implement static create() method.
gmx::RestraintMDModule::RestraintMDModule(std::unique_ptr<RestraintMDModuleImpl> restraint) :
    impl_{std::move(restraint)}
{
}
