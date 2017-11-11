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
#include <iostream>

#include "restraintfunctor-impl.h"

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

        void add(std::shared_ptr<::gmx::IRestraintPotential> restraint, std::string name);

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

        void virial(tensor virial)
        {
            virial_ = virial;
        }

        void communicator(const t_commrec& commRec)
        {
            communicator_ = &commRec;
        }

        std::shared_ptr<ICalculation> calculate(double t);

        std::shared_ptr<LegacyPuller> puller_;
        std::shared_ptr<::gmx::IRestraintPotential> restraint_;

    private:
        double currentTime_{0};
        double previousTime_{0};

        const t_mdatoms* atoms_{nullptr};
        const t_pbc* pbc_{nullptr};
        real lambda_{0};
        const rvec* positions_{nullptr};
        rvec* forces_{nullptr};
        rvec* virial_{nullptr};

        const t_commrec* communicator_{nullptr};

};

std::shared_ptr<ICalculation> ManagerImpl::calculate(double t)
{
    assert(communicator_ != nullptr);
    assert(atoms_ != nullptr);
    assert(pbc_ != nullptr);
    assert(getLegacy() != nullptr);
    assert(getLegacy()->getRaw() != nullptr);
    assert(forces_ != nullptr);
    assert(virial_ != nullptr);

    auto calculation = std::make_shared<Calculation>(t, *communicator_, *atoms_, *pbc_, lambda_, positions_, getLegacy()->getRaw(), forces_, virial_);
    if (calculation != nullptr)
    {
        previousTime_ = currentTime_;
        currentTime_ = t;
    }
    std::cout << "Manager calculating at time " << t << std::endl;
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

void ManagerImpl::add(std::shared_ptr<::gmx::IRestraintPotential> restraint, std::string name)
{

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
    if (impl_->restraint_ != nullptr)
    {
        energetic = true;
    };
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
    assert(impl_ != nullptr);
    impl_->atoms(&atoms);
}

void Manager::setBoundaryConditionsSource(const t_pbc &pbc)
{
    assert(impl_ != nullptr);
    impl_->pbc(&pbc);
}

void Manager::setPositionsSource(const rvec &x)
{
    assert(impl_ != nullptr);
    impl_->positions(&x);
}

void Manager::setForceOwner(rvec *f)
{
    assert(impl_ != nullptr);
    impl_->forces(f);
}

void Manager::setVirialOwner(tensor virial_force)
{
    assert(impl_ != nullptr);
    impl_->virial(virial_force);
}

void Manager::setLambdaSource(real lambda)
{
    assert(impl_ != nullptr);
    impl_->lambda(lambda);
}

void Manager::setCommunicator(const t_commrec &commRec)
{
    assert(impl_ != nullptr);
    impl_->communicator(commRec);
}

void Manager::add(std::shared_ptr<LegacyPuller> puller, std::string name)
{
    assert(impl_ != nullptr);
    impl_->addLegacy(std::move(puller), std::move(name));
}

void Manager::add(std::shared_ptr<gmx::IRestraintPotential> puller,
                  std::string name)
{
    assert(impl_ != nullptr);
    impl_->add(std::move(puller), name);
}

//void Manager::add(std::shared_ptr<gmx::IRestraintPotential> puller,
//                  std::string name)
//{
//    (void)puller;
//    (void)name;
//}

} // end namespace restraint
} // end namespace gmx
