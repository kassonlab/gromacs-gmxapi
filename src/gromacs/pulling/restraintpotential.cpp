//
// Created by Eric Irrgang on 9/23/17.
//

#include "restraintpotential.h"

#include <cassert>

#include <vector>

#include "gromacs/compat/make_unique.h"

class PotentialContainer::Impl
{
    public:
        std::vector<std::shared_ptr<gmx::RestraintPotential>> pullers_;
};

void PotentialContainer::addPotential(std::shared_ptr<gmx::RestraintPotential> puller) noexcept
{
    assert(impl_ != nullptr);
    impl_->pullers_.emplace_back(std::move(puller));
}

PotentialContainer::PotentialContainer() :
        impl_{gmx::compat::make_unique<PotentialContainer::Impl>()}
{}

PotentialContainer::~PotentialContainer() = default;

PotentialContainer &PotentialContainer::operator=(PotentialContainer &&) noexcept = default;

PotentialContainer::PotentialContainer(PotentialContainer&&) noexcept = default;

gmx::LegacyPullingPack::LegacyPullingPack(pull_t *pullWorkPointer) : pullWorkPointer_(pullWorkPointer)
{};

gmx::LegacyPullingPack::LegacyPullingPack(const gmx::LegacyPullingPack &source) : LegacyPullingPack(source.pullWorkPointer_)
{

}

gmx::LegacyPullingPack &gmx::LegacyPullingPack::operator=(const gmx::LegacyPullingPack &source)
{
    this->pullWorkPointer_ = source.pullWorkPointer_;
    return *this;
}

gmx::LegacyPullingPack::LegacyPullingPack(gmx::LegacyPullingPack &&old) noexcept : LegacyPullingPack(old.pullWorkPointer_)
{
    old.pullWorkPointer_ = nullptr;
}

gmx::LegacyPullingPack &gmx::LegacyPullingPack::operator=(gmx::LegacyPullingPack &&old) noexcept
{
    this->pullWorkPointer_ = old.pullWorkPointer_;
    old.pullWorkPointer_ = nullptr;
    return *this;
}
