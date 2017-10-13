//
// Created by Eric Irrgang on 9/23/17.
//

#ifndef GMX_PULLING_RESTRAINTPOTENTIAL_H
#define GMX_PULLING_RESTRAINTPOTENTIAL_H

#include "vectortype.h"

#include <memory>

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
 * \brief Base class for providers of pull_potential()
 *
 * For a set of \f$n\f$ coordinates, generate a force field according to a
 * scalar potential that is a fun. \f$F_i = - \grad_{q_i} \Phi (q_0, q_1, ... q_n; t)\f$
 *
 * Potentials implemented with these classes may be long ranged and are appropriate
 * for only a small number of particles to avoid substantial performance impact.
 * The atom indices appearing in the potential are provided by the user before the
 * simulation is run. The potential function is evaluated with a time argument
 * and can be updated during the simulation.
 *
 */
class RestraintPotential
{
    public:
        virtual ~RestraintPotential() = default;

        /*!
         * \brief Calculate a force vector according to two input positions.
         *
         * If not overridden by derived class, returns a zero vector.
         * \param r1 position of first site
         * \param r2 position of second site
         * \return force vector to be applied by calling code.
         */
        virtual vec3<real> calculateForce(vec3<real> r1, vec3<real> r2);
};

/*!
 * \brief Encapsulate the old pulling schemes.
 *
 * Give the old and new pulling code the same management and wrap up the old set of
 * schemes into a single pulling class. If the input record indicates a pulling protocol,
 * the RestraintPotential may use the associated resources (pull_t, pull_params_t, ...).
 */
class LegacyPullingPack final : public RestraintPotential
{
    public:
        ~LegacyPullingPack() override = default;

        /*!
         * \brief Copy the metadata, but manage the same internals.
         *
         * \param source object to copy
         */
        LegacyPullingPack(const LegacyPullingPack& source);

        /*!
         * \brief Copy the metadata, but manage the same internals.
         *
         * \param source object to copy
         */
        LegacyPullingPack& operator=(const LegacyPullingPack& source);

        /*!
         * \brief Move ownership of pull code manager.
         *
         * \param old object to move from.
         */
        LegacyPullingPack(LegacyPullingPack&& old) noexcept;

        /*!
         * \brief Move ownership of pull code manager.
         *
         * \param old move ownership from rhs of operation.
         */
        LegacyPullingPack& operator=(LegacyPullingPack&& old) noexcept;

        /*!
         * \brief Construct manager for legacy pulling code.
         *
         * \param pullWorkPointer pointer to struct created by init_pull() and managed by caller.
         */
        explicit LegacyPullingPack(pull_t* pullWorkPointer);

    private:
        struct pull_t* pullWorkPointer_;
};

} // end namespace gmx

/*!
 * \brief Hold an arbitrary number of PullPotential objects.
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
        void addPotential(std::shared_ptr<gmx::RestraintPotential> puller) noexcept;
    private:
        /// Private implementation class.
        class Impl;
        /// Implementation object.
        std::unique_ptr<Impl> impl_;
};


#endif //GMX_PULLING_PULLPOTENTIAL_H
