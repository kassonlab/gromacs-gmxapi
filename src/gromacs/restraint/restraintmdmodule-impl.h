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
#include <gromacs/gmxlib/network.h>
#include <gromacs/pbcutil/pbc.h>
#include <gromacs/mdtypes/mdatom.h>
#include "gromacs/math/vec.h"
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
            const auto site1 = static_cast<size_t>(site1_.index());
            const auto site2 = static_cast<size_t>(site2_.index());

            // Cooperatively get Cartesian coordinates for center of mass of each site
            RVec r1{0,0,0};
            RVec r2{0,0,0};
            if (cr->dd != nullptr)
            {
                // Get global-to-local indexing structure
                auto crossRef = cr->dd->ga2la;
                assert(crossRef != nullptr);
                int localIndex{-1};
                if (ga2la_get_home(crossRef,
                                   static_cast<int>(site1),
                                   &localIndex))
                {
                    assert(localIndex < mdatoms->homenr);
                    assert(localIndex >= 0);
                    // If atom is local, get its location
                    copy_rvec(x[localIndex], r1);
                }
                else
                {
                    // leave position == [0,0,0]
                }
                if (ga2la_get_home(crossRef,
                                   static_cast<int>(site2),
                                   &localIndex))
                {
                    assert(localIndex < mdatoms->homenr);
                    assert(localIndex >= 0);
                    // If atom is local, get its location
                    copy_rvec(x[localIndex], r2);
                }
                else
                {
                    // leave position == [0,0,0]
                }
            }
            else
            {
                // No DD so all atoms are local.
                copy_rvec(x[site1], r1);
                copy_rvec(x[site2], r2);
            }
            // r1 and r2 are now correct if local or [0,0,0] if not local.

            if (cr != nullptr && DOMAINDECOMP(cr))
            {
                // our quick-and-dirty short-term solution just relies on getting non-zero positions on exactly one ranks
                std::array<double, 6> buffer {{static_cast<double>(r1[0]),
                                                  static_cast<double>(r1[1]),
                                                  static_cast<double>(r1[2]),
                                                  static_cast<double>(r2[0]),
                                                  static_cast<double>(r2[1]),
                                                  static_cast<double>(r2[2]),
                                              }};
                // \todo This definitely needs some abstraction and checks.
                gmx_sumd(6, buffer.data(), cr);
                assert((r1[0] == 0) || (buffer[0] == r1[0]));
                assert((r1[1] == 0) || (buffer[1] == r1[1]));
                assert((r1[2] == 0) || (buffer[2] == r1[2]));
                assert((r2[0] == 0) || (buffer[3] == r2[0]));
                assert((r2[1] == 0) || (buffer[4] == r2[1]));
                assert((r2[2] == 0) || (buffer[5] == r2[2]));
                r1[0] = static_cast<real>(buffer[0]);
                r1[1] = static_cast<real>(buffer[1]);
                r1[2] = static_cast<real>(buffer[2]);
                r2[0] = static_cast<real>(buffer[3]);
                r2[1] = static_cast<real>(buffer[4]);
                r2[2] = static_cast<real>(buffer[5]);
            }

            // Apply minimum image convention to get into the same coordinate system.
            // \todo allow reference site for distances greater than half a box length.
            // \todo make conditional on what the implemented potential wants.
            {
                assert(check_box(-1, box) == nullptr);
                RVec dx{0,0,0};
                t_pbc pbc{};
                set_pbc(&pbc, -1, box);
                pbc_dx(&pbc, r2, r1, dx);
                rvec_add(r1, dx, r2);
            }

            auto result = restraint_->evaluate(make_vec3<real>(r1[0], r1[1], r1[2]),
                                               make_vec3<real>(r2[0], r2[1], r2[2]),
                                               t);

            size_t aLocal{site1};
            // Set forces using index a if no domain decomposition, otherwise set with local index if available.
            if ((cr->dd == nullptr) || ga2la_get_home(cr->dd->ga2la,
                               static_cast<int>(site1),
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
            size_t bLocal{site2};
            if ((cr->dd == nullptr) || ga2la_get_home(cr->dd->ga2la,
                               static_cast<int>(site2),
                               reinterpret_cast<int*>(&bLocal)))
            {
                force[bLocal][0] -= result.force.x;
                force[bLocal][1] -= result.force.y;
                force[bLocal][2] -= result.force.z;
            }

            if (int(t*1000) % 100 == 0)
            {
                if ((cr->dd == nullptr) || MASTER(cr))
                {
                    std::cout << "Evaluated restraint forces on atoms at " << make_vec3<real>(r1[0], r1[1], r1[2]) << " and " << make_vec3<real>(r2[0], r2[1], r2[2]) << ": " << result.force << ". rank,time: " << cr->rank_pp_intranode << "," << t << "\n";
                }
            }

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

/*! \internal
 * \brief IMDModule implementation for RestraintMDModule.
 *
 * Provides IMDModule interface.
 *
 * \ingroup module_restraint
 */
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
