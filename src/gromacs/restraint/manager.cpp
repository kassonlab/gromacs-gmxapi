//
// Created by Eric Irrgang on 10/25/17.
//

/*! \internal
 * \brief Implement the restraint manager.
 *
 * \ingroup module_restraint
 */

#include "manager.h"

#include <cassert>
#include <memory>
#include <gromacs/utility/basedefinitions.h>
#include <gromacs/pulling/pull.h>
#include "gromacs/compat/make_unique.h"

namespace gmx
{

class LegacyPuller;

namespace restraint
{

// Initialize static members
std::shared_ptr<Manager> Manager::instance_{nullptr};
std::mutex Manager::initializationMutex_{};

/*!
 * \brief Implementation class for restraint manager.
 */
class ManagerImpl
{
    public:
        void addLegacy(std::shared_ptr<LegacyPuller> puller, std::string name);
        std::shared_ptr<LegacyPuller> getLegacy();
    private:
        std::shared_ptr<LegacyPuller> puller_;
};

void ManagerImpl::addLegacy(std::shared_ptr<LegacyPuller> puller,
                            std::string name)
{
    puller_ = std::move(puller);
    (void)name;
}

std::shared_ptr<LegacyPuller> ManagerImpl::getLegacy()
{
    return puller_;
}

Manager::Manager() : impl_(gmx::compat::make_unique<ManagerImpl>()) {};

Manager::~Manager() = default;

std::shared_ptr<Manager> Manager::instance()
{
    std::lock_guard<std::mutex> lock(initializationMutex_);
    if (instance_ == nullptr)
    {
        // What do we want to do if `new` throws?
        instance_ = std::shared_ptr<Manager>(new Manager);
    }
    assert(instance_ != nullptr);
    return instance_;
}

void Manager::print(gmx_int64_t step,
                           double time)
{
    assert(impl_ != nullptr);
    if (auto puller = impl_->getLegacy())
    {
        if (auto pull_work = puller->getRaw())
        {
            pull_print_output(pull_work, step, time);
        }
    }
}

void Manager::finish()
{
    finish_pull(impl_->getLegacy()->getRaw());
}

pull_t *Manager::getRaw()
{
    assert(impl_ != nullptr);
    pull_t* pull_work{nullptr};
    if (auto puller = impl_->getLegacy())
    {
        pull_work = puller->getRaw();
    }
    return pull_work;
}

void Manager::makeLocalGroups(t_commrec *cr,
                              t_mdatoms *mdatoms)
{
    assert(impl_ != nullptr);
    if(auto puller = impl_->getLegacy())
    {
        if(auto pull_work = puller->getRaw())
        {
            dd_make_local_pull_groups(cr, impl_->getLegacy()->getRaw(), mdatoms);
        }
    }
}

bool Manager::contributesEnergy()
{
    assert(impl_ != nullptr);
    bool energetic{false};
    if (auto puller = impl_->getLegacy())
    {
        if (auto pull_work = puller->getRaw())
        {
            energetic = bool(pull_have_potential(impl_->getLegacy()->getRaw()));
        }
    }
    return energetic;
}

void Manager::clearConstraintForces()
{
    assert(impl_ != nullptr);
    if (auto puller = impl_->getLegacy())
    {
        auto pull_work = impl_->getLegacy()->getRaw();
        assert (pull_work != nullptr);
        if (pull_have_constraint(pull_work))
        {
            clear_pull_forces(pull_work);
        }
    }

}

void Manager::add(std::shared_ptr<LegacyPuller> puller, std::string name)
{
    assert(impl_ != nullptr);
    impl_->addLegacy(std::move(puller), std::move(name));
}

void Manager::add(std::shared_ptr<gmx::IRestraintPotential> puller,
                  std::string name)
{
    (void)puller;
    (void)name;
}

} // end namespace restraint
} // end namespace gmx
