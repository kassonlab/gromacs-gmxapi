//
// Created by Eric Irrgang on 10/25/17.
//

#ifndef GROMACS_RESTRAINT_MANAGER_H
#define GROMACS_RESTRAINT_MANAGER_H

/*! \libinternal \file
 * \brief Declare the Manager for restraint potentials.
 *
 * \ingroup module_restraint
 */

#include <memory>
#include <mutex>
#include <string>
#include "gromacs/utility/basedefinitions.h"
#include "restraintpotential.h"

struct t_commrec;
struct t_mdatoms;
struct pull_t;

namespace gmx
{
class LegacyPuller;

/*!
 * \brief Implementation details for MD restraints
 */
namespace restraint
{

class ManagerImpl;

/*!
 * \brief Manage the Restraint potentials available for Molecular Dynamics.
 *
 * Until further factoring of the MD integrators and force calculations, we use a singleton
 * to reduce coupling between rapidly changing GROMACS components. Ultimately, this manager
 * should either not be necessary or can be used in more tightly scoped instances.
 */
class Manager final
{
    public:
        ~Manager();
        /// Get a shared reference to the global manager.
        static std::shared_ptr<Manager> instance();

        Manager(const Manager&) = delete;
        Manager& operator=(const Manager&) = delete;
        Manager(Manager&&) = delete;
        Manager& operator=(Manager&&) = delete;

        void add(std::shared_ptr<LegacyPuller> puller, std::string name);
        void add(std::shared_ptr<gmx::IRestraintPotential> puller, std::string name);

        void print(gmx_int64_t step, double time);

        void finish();

        pull_t* getRaw();

        /*!
         * \brief Callback used when domain decomposition is able to provide local pull groups.
         */
        void makeLocalGroups(t_commrec* cr, t_mdatoms* mdatoms);

        /*!
         * \brief Whether managed restraints affect calculated potential energy.
         *
         * \return true if restraints are present that are relevant to potential energy calculation.
         */
        bool contributesEnergy();

        /*!
         * \brief Clear forces provided by constraints.
         *
         * If constraints are present, reset their force contribution. If they are not, do nothing.
         */
        void clearConstraintForces();

    private:
        Manager();

        static std::mutex initializationMutex_;
        static std::shared_ptr<Manager> instance_;

        std::unique_ptr<ManagerImpl> impl_;
};


} // end namespace restraint

} // end namespace gmx

#endif //GROMACS_RESTRAINT_MANAGER_H
