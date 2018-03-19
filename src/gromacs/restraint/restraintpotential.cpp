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
        std::vector<std::shared_ptr<gmx::IRestraintPotential>> pullers_;
};

void PotentialContainer::addPotential(std::shared_ptr<gmx::IRestraintPotential> puller) noexcept
{
    assert(impl_ != nullptr);
    impl_->pullers_.emplace_back(std::move(puller));
}

PotentialContainer::RestraintIterator PotentialContainer::begin()
{
    return impl_->pullers_.begin();
}

PotentialContainer::RestraintIterator PotentialContainer::end()
{
    return impl_->pullers_.end();
}

PotentialContainer::PotentialContainer() :
        impl_{gmx::compat::make_unique<PotentialContainer::Impl>()}
{}

//template<typename T>
//std::function<gmx::PotentialPointData(const gmx::Vector &,
//                                           const gmx::Vector &,
//                                           gmx::Time)> gmx::RestraintPotential<T>::getEvaluator()
//{
//    return nullptr;
//}

PotentialContainer::~PotentialContainer() = default;

PotentialContainer &PotentialContainer::operator=(PotentialContainer &&) noexcept = default;

PotentialContainer::PotentialContainer(PotentialContainer&&) noexcept = default;

namespace gmx
{

LegacyPuller::LegacyPuller(pull_t *pullWorkPointer) : pullWorkPointer_(pullWorkPointer)
{};

LegacyPuller::LegacyPuller(const LegacyPuller &source) : LegacyPuller(source.pullWorkPointer_)
{

}

LegacyPuller &LegacyPuller::operator=(const LegacyPuller &source)
{
    this->pullWorkPointer_ = source.pullWorkPointer_;
    return *this;
}

LegacyPuller::LegacyPuller(LegacyPuller &&old) noexcept : LegacyPuller(old.pullWorkPointer_)
{
    old.pullWorkPointer_ = nullptr;
}

LegacyPuller &LegacyPuller::operator=(LegacyPuller &&old) noexcept
{
    this->pullWorkPointer_ = old.pullWorkPointer_;
    old.pullWorkPointer_ = nullptr;
    return *this;
}

PotentialPointData LegacyPuller::calculate(Vector r1,
                                          Vector r2,
                                          Time t)
{
    (void)(r1);
    (void)(r2);
    (void)(t);
    return {};
}

struct pull_t *LegacyPuller::getRaw()
{
    return pullWorkPointer_;
}

const struct pull_t *LegacyPuller::getRaw() const
{
    return pullWorkPointer_;
}

void IRestraintPotential::update(gmx::Vector v,
                                 gmx::Vector v0,
                                 double t)
{
    // By default, an IRestraintPotential has a null updater.
    (void)(v);
    (void)(v0);
    (void)(t);
}
} // end namespace gmx
