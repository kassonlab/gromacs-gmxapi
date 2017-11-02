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
class ICalculation;

/*!
 * \brief Manage the Restraint potentials available for Molecular Dynamics.
 *
 * Until further factoring of the MD integrators and force calculations, we use a singleton
 * to reduce coupling between rapidly changing GROMACS components. Ultimately, this manager
 * should either not be necessary or can be used in more tightly scoped instances.
 *
 * The manager takes ownership of the "pull groups" (or atomic selections) and of
 * the various restraints and constraints applied for a given simulation.
 *
 * Calling code provides the manager with a means to access the various required input data
 * to be used when restraints are computed.
 *
 * \internal
 * When calculate(t) is called, the various input and output data sources are provided to
 * a CalculationBuilder to produce a Calculation object for the point in time, t.
 * Constructing the Calculation object triggers updates to the simulation state force array
 * and virial tensor. After construction, the Calculation object can be queried for calculated
 * data such as energy or pulling work.
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
//        void add(std::shared_ptr<gmx::IRestraintPotential> puller, std::string name);

        /*!
         * \brief Provide restraints with a source of atom information.
         *
         * Until additional factoring and clarification, restraints manager has no way to know if
         * the previously used pointer is valid or who to ask for it, so a pointer must be provided
         * to the working data structure for atom data. At the very least, there could be an object
         * from which the restraints manager _requests_ a fresh handle, but right now the raw pointer
         * is all we have.
         *
         * \param atoms various quick-reference data for atoms in system.
         */
        void setAtomsSource(const t_mdatoms& atoms);

        /*!
         * \brief Provide restraints with a source of information on periodic boundary conditions.
         *
         * Until additional factoring and clarification, restraints manager has no way to know if
         * the previously used pointer is valid or who to ask for it, so a pointer must be provided
         * to the structure providing more complete box information. At the very least, there could be an object
         * from which the restraints manager _requests_ a fresh handle, but right now the raw pointer
         * is all we have.
         *
         * \param pbc periodic boundary conditions and simulation box data
         */
        void setBoundaryConditionsSource(const t_pbc& pbc);

        /*!
         * \brief Provide restraints with the current communicator handle.
         *
         * Until additional factoring and clarification, restraints manager has no way to know if
         * the previously used pointer is valid or who to ask for it, so a raw pointer must be provided
         * to the communcation record. At the very least, there could be an object
         * from which the restraints manager _requests_ a fresh handle with access controls, but right now the raw pointer
         * is all we have.
         *
         * \param commRec communications record
         */
        void setCommunicator(const t_commrec& commRec);

        /*!
         * \brief Provide restraints with a source of atomic coordinates.
         *
         * Until additional factoring and clarification, restraints manager has no way to know if
         * the previously used pointer is valid or who to ask for it, so a pointer must be provided.
         * At the very least, there could be an object
         * from which the restraints manager _requests_ a fresh handle, but right now the raw pointer
         * is all we have.
         * \param x positions
         */
        void setPositionsSource(const rvec& x);

        /*!
         * \brief Provide restraints with a connection to the integrator force data.
         *
         * Until additional factoring and clarification, restraints manager has no way to know if
         * the previously used pointer is valid or who to ask for it, so a pointer must be provided.
         * At the very least, there could be an object
         * from which the restraints manager _requests_ a fresh handle, but right now the raw pointer
         * is all we have.
         * \param f read/write handle to force data.
         */
        void setForceOwner(rvec* f);

        /*!
         * \brief Provide restraints with a connection to the integrator virial data.
         *
         * Until additional factoring and clarification, restraints manager has no way to know if
         * the previously used pointer is valid or who to ask for it, so a pointer must be provided. At the very least, there could be an object
         * from which the restraints manager _requests_ a fresh handle, but right now the raw pointer
         * is all we have.
         * \param virial_force read/write handle with which to correct for the virial contribution.
         */
        void setVirialOwner(tensor virial_force);

        /*!
         * \brief Provide restraints with a source of current lambda value.
         *
         * Until additional factoring and clarification, restraints manager has no way to retrieve the current value on its own or who to ask for it, so the current value must be provided in addition to the current time value. At the very least, there could be an object
         * from which the restraints manager _requests_ a fresh handle, but right now the raw pointer
         * is all we have.
         * \param lambda
         */
        void setLambdaSource(real lambda);

        /*!
         * \brief (Re)calculate the restraint forces for time t.
         *
         * Refreshes the values available to the accessors for energy, force,
         * and virial.
         * \param t Simulation time (coordinate value) in picoseconds.
         *
         * Note: this call could be made implicit or could be triggered by some
         * other signal or work schedule.
         *
         * Fulfills the role of pull_potential() defined in pull.cpp
         */
        std::shared_ptr<ICalculation> calculate(double t);


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
