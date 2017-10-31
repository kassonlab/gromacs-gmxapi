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
#include <map>

#include "restraintfunctor-impl.h"

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/pulling/pull.h"
#include "gromacs/compat/make_unique.h"
#include "restraintcalculation-impl.h"

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
        std::tuple<double, real> energy(double time) const noexcept;
        std::tuple<double, real> work(double time) const noexcept;

        void currentTime(double currentTime)
        {
            currentTime_ = currentTime;
        }

        void previousTime(double previousTime)
        {
            previousTime_ = previousTime;
        }

        void atoms(const t_mdatoms *atoms)
        {
            atoms_ = atoms;
        }

        void pbc(const t_pbc *pbc)
        {
            pbc_ = pbc;
        }

        void lambda(real lambda)
        {
            lambda_ = lambda;
        }

        void positions(rvec const *positions)
        {
            positions_ = positions;
        }

        void forces(rvec *forces)
        {
            forces_ = forces;
        }

        void virial(tensor *virial)
        {
            virial_ = virial;
        }

        std::shared_ptr<ICalculation> calculate(double t);

    private:
        double currentTime_{0};
        double previousTime_{0};

        const t_mdatoms* atoms_{nullptr};
        const t_pbc* pbc_{nullptr};
        real lambda_{0};
        const rvec* positions_{nullptr};
        rvec* forces_{nullptr};
        tensor* virial_{nullptr};

        std::shared_ptr<LegacyPuller> puller_;
};

std::shared_ptr<ICalculation> ManagerImpl::calculate(double t)
{
    auto calculation = std::make_shared<Calculation>(t, *atoms_, *pbc_, lambda_, positions_, forces_, virial_);
    if (calculation != nullptr)
    {
        previousTime_ = currentTime_;
        currentTime_ = t;
    }
    return calculation;
}

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
    assert(impl_ != nullptr);
    auto puller = impl_->getLegacy();
    if (puller)
    {
        finish_pull(impl_->getLegacy()->getRaw());
    }
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
            dd_make_local_pull_groups(cr, pull_work, mdatoms);
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
            energetic = bool(pull_have_potential(pull_work));
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

std::shared_ptr<ICalculation> Manager::calculate(double t)
{
    // Disambiguate the input and output parameters of the old pull_potential by
    // constructing a functor and then applying it.
//    auto restraints = ::gmx::restraint::RestraintFunctor(t, const t_mdatoms& md, const t_pbc& pbc, real lambda);
//    restraints.calculate();

    assert(impl_ != nullptr);
    auto calculation = impl_->calculate(t);
    return calculation;
}

void Manager::setAtomsSource(const t_mdatoms &atoms)
{
    impl_->atoms(&atoms);
}

void Manager::setBoundaryConditionsSource(const t_pbc &pbc)
{
    impl_->pbc(&pbc);
}

void Manager::setPositionsSource(const rvec &x)
{
    impl_->positions(&x);
}

void Manager::setForceOwner(rvec *f)
{
    impl_->forces(f);
}

void Manager::setVirialOwner(tensor *virial_force)
{
    impl_->virial(virial_force);
}

void Manager::setLambdaSource(real lambda)
{
    impl_->lambda(lambda);
}

//void Manager::add(std::shared_ptr<gmx::IRestraintPotential> puller,
//                  std::string name)
//{
//    (void)puller;
//    (void)name;
//}

} // end namespace restraint
} // end namespace gmx
