//
// Created by Eric Irrgang on 9/23/17.
//

#ifndef GMX_PULLING_RESTRAINTPOTENTIAL_H
#define GMX_PULLING_RESTRAINTPOTENTIAL_H

/*! \file
 * \defgroup modules_pulling
 * \brief Declare generic interface for restraint implementations.
 *
 */

#include "vectortype.h"

#include <memory>
#include <vector>

#include "gromacs/math/vectypes.h"
//#include "gromacs/mdtypes/pull-params.h"
//#include "gromacs/pulling/pull_internal.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

struct gmx_mtop_t;
struct gmx_output_env_t;
struct pull_params_t;
struct pull_t;
struct t_commrec;
struct t_filenm;
struct t_inputrec;
struct t_mdatoms;
struct t_pbc;

namespace gmx
{

using detail::vec3;

/*!
 * \brief Typed unitless time in GROMACS.
 *
 * It may be helpful to explicitly specify the units of time.
 *
 * Example:
 *
 *     gmx_time time;
 *     // float t = time; // error: no implicit conversion because units are ambiguous.
 *     using gmx::time::ps;
 *     float t = time * ps;
 *
 * TODO: If this persists, it will be moved from this header.
 */
struct gmx_time
{
    using float_type = real;
    float_type t;
};

/*!
 * \brief Structure to hold the results of IRestraintPotential::evaluate().
 */
class PotentialPointData
{
    public:
        vec3<real> force{};
        real energy{0.0};
};

/*!
 * \brief Interface for Restraint potentials.
 *
 * Derive from this interface class to implement a restraint potential. The derived class
 * must implement the evaluate() member function that produces a PotentialPointData instance.
 *
 * For a set of \f$n\f$ coordinates, generate a force field according to a
 * scalar potential that is a fun. \f$F_i = - \nabla_{q_i} \Phi (q_0, q_1, ... q_n; t)\f$
 *
 * Potentials implemented with these classes may be long ranged and are appropriate
 * for only a small number of particles to avoid substantial performance impact.
 * The atom indices appearing in the potential are provided by the user before the
 * simulation is run. The potential function is evaluated with a time argument
 * and can be updated during the simulation.
 */
class IRestraintPotential
{
    public:
        virtual ~IRestraintPotential() = default;

        /*!
         * \brief Calculate a force vector according to two input positions at a given time.
         *
         * If not overridden by derived class, returns a zero vector.
         * \param r1 position of first site
         * \param r2 position of second site
         * \param t simulation time in picoseconds
         * \return force vector and potential energy to be applied by calling code.
         */
        virtual PotentialPointData evaluate(vec3<real> r1,
                              vec3<real> r2,
                              double t) = 0;
};

/*!
 * \brief Encapsulate the old pulling schemes.
 *
 * Give the old and new pulling code the same management and wrap up the old set of
 * schemes into a single pulling class. If the input record indicates a pulling protocol,
 * the RestraintPotential may use the associated resources (pull_t, pull_params_t, ...).
 */
class LegacyPuller : public IRestraintPotential
{
    public:
        PotentialPointData evaluate(vec3<real> r1,
                                    vec3<real> r2,
                                    double t) override;

        ~LegacyPuller() = default;

        /*!
         * \brief Copy the metadata, but manage the same internals.
         *
         * \param source object to copy
         */
        LegacyPuller(const LegacyPuller& source);

        /*!
         * \brief Copy the metadata, but manage the same internals.
         *
         * \param source object to copy
         */
        LegacyPuller& operator=(const LegacyPuller& source);

        /*!
         * \brief Move ownership of pull code manager.
         *
         * \param old object to move from.
         */
        LegacyPuller(LegacyPuller&& old) noexcept;

        /*!
         * \brief Move ownership of pull code manager.
         *
         * \param old move ownership from rhs of operation.
         */
        LegacyPuller& operator=(LegacyPuller&& old) noexcept;

        /*!
         * \brief Construct manager for legacy pulling code.
         *
         * \param pullWorkPointer pointer to struct created by init_pull() and managed by caller.
         */
        explicit LegacyPuller(pull_t* pullWorkPointer);

    private:
        // borrow access to pull_t owned by calling code.
        struct pull_t* pullWorkPointer_;
};

} // end namespace gmx

/*!
 * \brief Hold an arbitrary number of gmx::IRestraintPotential objects.
 */
class PotentialContainer
{
    public:
        /// Default constructor
        PotentialContainer();
        /// Copy constructor disabled.
        PotentialContainer(const PotentialContainer&) = delete;
        /// Move construction allowed.
        PotentialContainer(PotentialContainer&&) noexcept;
        /// Copy assignment disabled.
        PotentialContainer& operator=(const PotentialContainer&) = delete;
        /// Allow move assignment.
        PotentialContainer& operator=(PotentialContainer&&) noexcept;
        /// non-virtual destructor.
        ~PotentialContainer();

        /*!
         * \brief Add a pulling potential to the managed list.
         *
         * \param puller GROMACS-provided or custom pulling potential
         */
        void addPotential(std::shared_ptr<gmx::IRestraintPotential> puller) noexcept;

        using RestraintIterator = std::vector<std::shared_ptr<gmx::IRestraintPotential>>::iterator;
        RestraintIterator begin();
        RestraintIterator end();

    private:
        /// Private implementation class.
        class Impl;
        /// Implementation object.
        std::unique_ptr<Impl> impl_;
};


#endif //GMX_PULLING_PULLPOTENTIAL_H
