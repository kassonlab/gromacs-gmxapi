//
// Created by Eric Irrgang on 9/23/17.
//

#ifndef GMX_PULLING_RESTRAINTPOTENTIAL_H
#define GMX_PULLING_RESTRAINTPOTENTIAL_H

/*!
 * \defgroup module_restraint MD restraints
 * \ingroup group_mdrun
 * Apply restraints during MD integration.
 */
/*! \file
 * \brief Declare generic interface for restraint implementations.
 *
 * \ingroup module_restraint
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
 * \ingroup module_restraint
 */
struct gmx_time
{
    using float_type = real;
    float_type t;
};

/*!
 * \brief Structure to hold the results of IRestraintPotential::evaluate().
 *
 * \ingroup module_restraint
 */
class PotentialPointData
{
    public:

        using detail::vec3;

        /*!
         * \brief Force vector calculated for first position.
         */
        vec3<real> force;
        /*!
         * \brief Potential energy calculated for this interaction.
         */
        real energy;


        /*!
         * \brief Initialize a new data structure.
         */
        PotentialPointData() : PotentialPointData{vec3<real>(), real(0.0)}
        {};

        /*!
         * \brief Initialize from an argument list
         *
         * \param f Force vector.
         * \param e Energy value.
         *
         * Note that if force was calculated as a scalar, it needs to be multiplied by a unit vector in the direction to which it should be applied.
         * If this calculation is in a subclass of gmx::RestraintPotential,
         * you should be able to use the make_force_vec() helper function.
         */
        PotentialPointData(const vec3<real>& f, const real e) :
                force{f},
                energy{e}
        {};
};



/*!
 * \brief Interface for Restraint potentials.
 *
 * Derive from this interface class to implement a restraint potential. The derived class
 * must implement the evaluate() member function that produces a PotentialPointData instance.
 * For various convenience functions and to decouple the internal library
 * interface from implementations of potentials, it is expected that
 * restraints will be programmed by subclassing gmx::RestraintPotential<>
 * rather than gmx::IRestraintPotential.
 *
 * For a set of \f$n\f$ coordinates, generate a force field according to a
 * scalar potential that is a fun. \f$F_i = - \nabla_{q_i} \Phi (q_0, q_1, ... q_n; t)\f$
 *
 * Potentials implemented with these classes may be long ranged and are appropriate
 * for only a small number of particles to avoid substantial performance impact.
 *
 * The indices in the above equation refer to the input and output sites specified before the simulation is started.
 * In the current interface, the evaluate() virtual function allows an
 * implementer to calculate the energy and force acting at the first of a pair of sites. The equal and opposite force is applied at the second site.
 *
 * The potential function is evaluated with a time argument
 * and can be updated during the simulation.
 * For non-time-varying potentials, the time argument may still be useful for
 * internal optimizations, such as managing cached values.
 *
 * In the simplest and most common case, pairs of sites (atoms or groups)
 * are specified by the user and then, during integration, GROMACS provides
 * the current positions of each pair for the restraint potential to be evaluated.
 * In such a case, the potential can be implemented by overriding...
 *
 * \ingroup module_restraint
 */
class IRestraintPotential
{
    public:
        using detail::vec3;

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
 * \brief Allow easy implementation of restraint potentials
 *
 * \tparam T implementation type
 *
 * Implementation of type T uses RestraintPotential<T> as a mix-in. E.g.
 *
 *     class MyRestraint : public RestraintPotential<T>
 *     {
 *         // ...
 *     }
 *
 * \ingroup module_restraint
 */
template<class T>
class RestraintPotential : public IRestraintPotential
{
    public:

        /// template interface
        PotentialPointData calculate(vec3 <real> r1,
                                     vec3 <real> r2,
                                     double t);
        /*! \cond internal
         * \brief IRestraintPotential interface
         *
         * Implements the interface for MD force evaluation. Though the definition
         * is necessarily in the template header, this interface is less stable
         * and restraint potentials should be coded using the public interface as described.
         */
        PotentialPointData evaluate(vec3<real> r1,
                                    vec3<real> r2,
                                    double t) override;
        /*! \endcond */

    private:
        /*! \cond internal
         * \brief Allow construction only by the specific class T
         *
         * This avoids potentially hair-pulling behavior, say, in the event of
         * typos that use RestraintPotential<T> where RestraintPotential<U> was intended, such as the following.
         *
         *     class RestraintB : public RestraintPotential<A> {
         *     //...
         *     };
         * \{
         */
        friend class T;
        RestraintPotential() = default;
        /*! \} \endcond */
};

/*!
 * The obvious ways to resolve which method to call would be with SFINAE or some sort of static dispatch.
 *
 * \ingroup module_restraint
 */
template<class T>
PotentialPointData RestraintPotential<T>::evaluate(vec3<real> r1,
                                                vec3<real> r2,
                                                double t)
{
    return static_cast<T&>(*this).calculate(r1, r2, t);
}

template<class T>
PotentialPointData RestraintPotential<T>::calculate(vec3 <real> r1,
                                                    vec3 <real> r2,
                                                    double t)
{
    (void)(r1);
    (void)(r2);
    (void)(t);
    data_ = PotentialPointData{};
}

template <class T>
vec3<real> RestraintPotential<T>::force(double t) const
{
    (void)(t);
    return data_.force;
}

template <class T>
real RestraintPotential<T>::energy(double t) const
{
    (void)(t);
    return data_.energy;
}

/*!
 * \brief Encapsulate the old pulling schemes.
 *
 * Give the old and new pulling code the same management and wrap up the old set of
 * schemes into a single pulling class. If the input record indicates a pulling protocol,
 * the RestraintPotential may use the associated resources (pull_t, pull_params_t, ...).
 *
 * \ingroup module_restraint
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
 *
 * \ingroup module_restraint
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
