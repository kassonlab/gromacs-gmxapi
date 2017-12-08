//
// Created by Eric Irrgang on 11/10/17.
//

#ifndef GROMACS_RESTRAINTMDMODULE_IMPL_H
#define GROMACS_RESTRAINTMDMODULE_IMPL_H

/*! \internal \file
 * \brief Implementation details for RestraintMDModule
 *
 *
 * \ingroup module_restraint
 */

#include "restraintpotential.h"

#include <iostream>
#include <gromacs/mdtypes/commrec.h>
#include <gromacs/domdec/domdec_struct.h>
#include <gromacs/domdec/ga2la.h>
#include "gromacs/mdtypes/imdpoptionprovider.h"
#include "gromacs/mdtypes/imdoutputprovider.h"
#include "gromacs/mdtypes/iforceprovider.h"
#include "gromacs/mdtypes/imdmodule.h"
#include "gromacs/utility/arrayref.h"

namespace gmx
{

namespace
{
/*! \internal
 * \brief Abstraction for a restraint interaction site.
 *
 * A restraint may operate on a single atom or some other entity, such as a selection of atoms.
 * The Restraint implementation is very independent from how coordinates are provided or what they mean.
 *
 * First implementation can only represent single atoms in a global context.
 */
class Site
{
    public:
        /*! \brief Construct from global atom indices
         *
         * \param globalIndex Atom index in the global state (as input to the simulation)
         */
        explicit Site(unsigned long int globalIndex) : index_{globalIndex} {};

        unsigned long int index() const {return index_;};

    private:
        unsigned long int index_;
};

RVec sitePosition(const t_commrec *commRec,
                                     const Site& site,
                                     const rvec* const localPositions)
{
    int localIndex = static_cast<int>(site.index());
    RVec position{0, 0, 0};
    if (commRec->dd != nullptr)
    {
        // Get global-to-local indexing structure
        auto crossRef = commRec->dd->ga2la;
        assert(crossRef != nullptr);
        if (ga2la_get_home(crossRef,
                           static_cast<int>(site.index()),
                           &localIndex))
        {
            position = localPositions[localIndex];
        }
    }
    else
    {
        position = localPositions[localIndex];
    }
    return position;
};

} //end anonymous namespace

/*!
 * \brief Concrete MdpOptionProvider for Restraints
 */
class RestraintOptionProvider : public gmx::IMdpOptionProvider
{
    public:
        void initMdpTransform(IKeyValueTreeTransformRules *transform) override
        {
            (void)(transform);
        }

        void initMdpOptions(IOptionsContainerWithSections *options) override
        {
            (void)(options);
        }

        void buildMdpOutput(KeyValueTreeObjectBuilder *builder) const override
        {
            (void)(builder);
        }
};

/*!
 * \brief MDOutputProvider concrete class for Restraints
 */
class RestraintOutputProvider : public gmx::IMDOutputProvider
{
    public:
        void initOutput(FILE *fplog,
                        int nfile,
                        const t_filenm *fnm,
                        bool bAppendFiles,
                        const gmx_output_env_t *oenv) override
        {
            IMDOutputProvider::initOutput(fplog,
                                          nfile,
                                          fnm,
                                          bAppendFiles,
                                          oenv);
        }

        void finishOutput() override
        {
            IMDOutputProvider::finishOutput();
        }
};

/*!
 * \brief Provide IForceProvider for RestraintMDModuleImpl
 *
 * Adapter class from IForceProvider to IRestraintPotential.
 * Objects of this type are uniquely owned by instances of RestraintMDModuleImpl. The object will
 * dispatch calls to IForceProvider->calculateForces() to the functor managed by RestraintMDModuleImpl.
 */
class RestraintForceProvider : public gmx::IForceProvider
{
    public:
        /*!
         * \brief Can only be constructed when initialized from a restraint.
         */
        RestraintForceProvider() = delete;

        /*!
         * \brief RAII construction with an IRestraintPotential
         *
         * Note, this object must outlive the pointer that will be provided to ForceProviders.
         * \param restraint
         */
        explicit RestraintForceProvider(std::shared_ptr<gmx::IRestraintPotential> restraint,
                                        unsigned long int site1,
                                        unsigned long int site2);

        /*!
         * \brief Implement the IForceProvider interface.
         *
         * Update the force array with restraint contribution(s) for local atoms.
         *
         * RestraintForceProvider is implemented with the assumption that few restraints apply to many atoms. That is,
         * the number of restraints affecting a large number of atoms is small, though there may be several restraints
         * that apply to few atoms each. Under this assumption, it is considered computationally inexpensive to iterate
         * over restraints in an outer loop and iterate over atoms within each restraint. This would be an invalid assumption
         * if, say, several restraints applied to an entire membrane or the entire solvent group.
         *
         * If the assumption causes performance problems, we can look for a good way to reduce from several restraints
         * in a single pass or a very lightweight way to determine whether a given restraint applies to a given atom.
         * There is also the notion in the pulling code of a limited number of "pull groups" used by the "pull coordinates".
         * The right optimization will depend on how the code is being used, but I expect allocating and reusing even
         * large arrays for lookup tables and calculation staging areas will be effective.
         *
         * Call the evaluator(s) for the restraints for the configured
         * \param cr
         * \param mdatoms
         * \param box
         * \param t
         * \param x
         * \param force
         */
        void calculateForces(const t_commrec          *cr,
                             const t_mdatoms          *mdatoms,
                             const matrix              box,
                             double                    t,
                             const rvec               *x,
                             gmx::ArrayRef<gmx::RVec>  force)

        override
        {
            using gmx::detail::make_vec3;

            assert(restraint_ != nullptr);
            auto a = site1_.index();
            auto b = site2_.index();

            const RVec r1 = sitePosition(cr, site1_, x);
            const RVec r2 = sitePosition(cr, site2_, x);

            auto result = restraint_->evaluate(make_vec3<real>(r1[0], r1[1], r1[2]),
                                               make_vec3<real>(r2[0], r2[1], r2[3]),
                                               t);

            size_t aLocal{a};
            // Set forces using index a if no domain decomposition, otherwise set with local index if available.
            if ((cr->dd == nullptr) || ga2la_get_home(cr->dd->ga2la,
                               static_cast<int>(a),
                               reinterpret_cast<int*>(&aLocal)))
            {
                force[aLocal][0] += result.force.x;
                force[aLocal][1] += result.force.y;
                force[aLocal][2] += result.force.z;
            }

            // Note: Currently calculateForces is called once per restraint and each restraint
            // applies to a pair of atoms. Future optimizations may consolidate multiple restraints
            // with possibly duplicated sites, in which case we may prefer to iterate over non-frozen
            // sites to apply forces without explicitly expressing pairwise symmetry as in the
            // following logic.
            size_t bLocal{b};
            if ((cr->dd == nullptr) || ga2la_get_home(cr->dd->ga2la,
                               static_cast<int>(b),
                               reinterpret_cast<int*>(&bLocal)))
            {
                force[b][0] -= result.force.x;
                force[b][1] -= result.force.y;
                force[b][2] -= result.force.z;
            }

            std::cout << "Evaluated restraint forces on atoms at " << make_vec3<real>(r1[0], r1[1], r1[2]) << " and " << make_vec3<real>(r2[0], r2[1], r2[3]) << ": " << result.force << "\n";

        }
    private:
        std::shared_ptr<gmx::IRestraintPotential> restraint_;
        Site site1_;
        Site site2_;
};

RestraintForceProvider::RestraintForceProvider(std::shared_ptr<gmx::IRestraintPotential> restraint,
                                               unsigned long int site1,
                                               unsigned long int site2) :
    restraint_{std::move(restraint)},
    site1_{site1},
    site2_{site2}
{
    assert(restraint_ != nullptr);
}

class RestraintMDModuleImpl final: public gmx::IMDModule
{
    public:
        RestraintMDModuleImpl() = delete;
        RestraintMDModuleImpl(std::shared_ptr<gmx::IRestraintPotential>,
                              unsigned long int site1,
                              unsigned long int site2);

        // Does default move work right (exception safe) with multiple unique_ptr members?
        RestraintMDModuleImpl(RestraintMDModuleImpl&&) noexcept = default;
        RestraintMDModuleImpl& operator=(RestraintMDModuleImpl&&) noexcept = default;

        virtual ~RestraintMDModuleImpl();

        IMdpOptionProvider *mdpOptionProvider() override;

        IMDOutputProvider *outputProvider() override;

        /*!
         * \brief
         *
         * \param forceProviders force module manager in the force record that will call this.
         *
         * The calling code must ensure that this object stays alive as long as forceProviders needs
         * the RestraintForceProvider, since forceProviders can't. Typically that is the duration of a do_md() call.
         */
        void initForceProviders(ForceProviders *forceProviders) override;

        std::unique_ptr<RestraintForceProvider> forceProvider_;
        std::unique_ptr<RestraintOutputProvider> outputProvider_;
        std::unique_ptr<RestraintOptionProvider> optionProvider_;
};

} // end namespace gmx

#endif //GROMACS_RESTRAINTMDMODULE_IMPL_H
