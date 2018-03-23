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
#include <mutex>

#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/domdec/ga2la.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/iforceprovider.h"
#include "gromacs/mdtypes/imdmodule.h"
#include "gromacs/mdtypes/imdoutputprovider.h"
#include "gromacs/mdtypes/imdpoptionprovider.h"
#include "gromacs/pbcutil/pbc.h"
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
        explicit Site(unsigned long int globalIndex) :
            index_{globalIndex},
            t_{0.},
            r_{0, 0, 0}
        {};

        // Explicitly define copies. Implicit definition is not possible because of the mutex member,
        // and a copy constructor is necessary to use Site in a std::vector. Really, we should make
        // a copy point to the same implementation object to reuse its cache.
        Site(const Site& site) :
            index_{site.index_},
            t_{site.t_},
            r_{site.r_},
            cacheMutex_{}
        {}
        // Assignment doesn't make sense because it implies that a site's meaning is fuzzy.
        // If the definition of a site is changing, just make a new site. There's nothing to be gained
        // by reusing one or by creating it uninitialized.
        Site& operator=(const Site& site) = delete;

        unsigned long int index() const {return index_;};

        /*!
         * \brief Get the position of this site at time t.
         *
         * By providing the current time, we can cache results in order to use them once per timestep.
         * In the long term, we would prefer to also allow client code to preregister interest in a
         * position at a given time, or issue "futures".
         *
         * \param cr Communications record.
         * \param nx Number of locally available atoms (size of local atom data arrays)
         * \param x Array of locally available atom coordinates.
         * \param t the current time.
         * \return position vector.
         */
        RVec position(const t_commrec *cr, size_t nx, const rvec* x, double t)
        {
            // I think I'm overlooking a better way to do this, but I think we should move to an issuer
            // of futures instead of explicitly locking, so it might not be worth looking into. I'm trying
            // to avoid locking the mutex if possible, but if we wait for a lock, we might perform the
            // cache update several times, and yet it would be nice not to have to make any lock at all
            // if the cache is already up-to-date. The other thing we can do is to require that someone
            // trigger the Sites update at the global atom position update(s) and make it an error to
            // request site positions for the wrong time, since it is likely that all threads will request
            // site position data at almost exactly the same time. Finally, we could just use a barrier,
            // but then we have to make assumptions or keep track of which threads are interested.
            // Also, we should manage all sites at once instead of one at a time.
            if(t_ <= t) // Do we need to update the cache?
            {
                std::lock_guard<std::mutex> cacheLock(cacheMutex_);
                // Now that we've got the lock, check again whether the cache was updated while we were waiting.
                if (t_ <= t)
                {
                    gmx::RVec r{0,0,0};
                    if (cr != nullptr && DOMAINDECOMP(cr)) // Domain decomposition
                    {
                        // Get global-to-local indexing structure
                        auto crossRef = cr->dd->ga2la;
                        assert(crossRef != nullptr);
                        int localIndex{-1};
                        if (ga2la_get_home(crossRef,
                                           static_cast<int>(index_),
                                           &localIndex))
                        {
                            assert(localIndex < static_cast<decltype(localIndex)>(nx));
                            assert(localIndex >= 0);
                            // If atom is local, get its location
                            copy_rvec(x[localIndex], r);
                        }
                        else
                        {
                            // Nothing to contribute on this rank. Leave position == [0,0,0].
                        }
                        // AllReduce across the ranks of the simulation.
                        std::array<double , 3> buffer{{r[0], r[1], r[2]}};
                        // \todo This definitely needs some abstraction and checks.
                        // This should be an all-reduce sum, which gmx_sumd appears to be.
                        gmx_sumd(3, buffer.data(), cr);
                        r[0] = static_cast<real>(buffer[0]);
                        r[1] = static_cast<real>(buffer[1]);
                        r[2] = static_cast<real>(buffer[2]);

                    } // end domain decomposition branch
                    else
                    {
                        // No DD so all atoms are local.
                        copy_rvec(x[index_], r);
                    }
                    // Update cache and cache status.
                    copy_rvec(r, r_);
                    t_ = t;
                } // end inner check if (t_ <= t)
                // else the cache was updated while we were waiting for the lock. We're done here.
                // release lock at end of block.
            } // end outer check if (t_ <= t)

            return r_;
        }

    private:
        /// global index of the single-atom site.
        const unsigned long int index_;
        volatile double t_;
        RVec r_;
        std::mutex cacheMutex_;
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
                                        const std::vector<unsigned long int>& sites);

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
         * Call the evaluator(s) for the restraints for the configured sites. Forces are applied to atoms in the first
         * and last site listed. Intermediate sites are used as reference coordinates when the relevant vector between
         * sites is on the order of half a box length or otherwise ambiguous in the case of periodic boundary conditions.
         *
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
            using gmx::detail::vec3;
            using gmx::detail::make_vec3;

            assert(restraint_ != nullptr);

            assert(mdatoms->homenr >= 0);

            assert(check_box(-1, box) == nullptr);
            t_pbc pbc{};
            set_pbc(&pbc, -1, box);

            // Cooperatively get Cartesian coordinates for center of mass of each site
            RVec r1 = sites_[0].position(cr,
                                      static_cast<size_t>(mdatoms->homenr),
                                      x,
                                      t);
            RVec r2{r1[0],r1[1],r1[2]};
            rvec dr{0,0,0};
            // Build r2 by following a path of difference vectors that are each presumed to be less than
            // a half-box apart, in case we are battling periodic boundary conditions along the lines of
            // a big molecule in a small box.
            for (size_t i = 0; i < sites_.size() - 1; ++i)
            {
                RVec a = sites_[i].position(cr,
                                            static_cast<size_t>(mdatoms->homenr),
                                            x,
                                            t);
                RVec b = sites_[i+1].position(cr,
                                              static_cast<size_t>(mdatoms->homenr),
                                              x,
                                              t);
                pbc_dx(&pbc, b, a, dr);
                r2[0] += dr[0];
                r2[1] += dr[1];
                r2[2] += dr[2];
            }
            // In the case of a single-atom site, r1 and r2 are now correct if local or [0,0,0] if not local.


            // Master rank update call-back. This needs to be moved to a discrete place in the
            // time step to avoid extraneous barriers. The code would be prettier with "futures"...
            if ((cr->dd == nullptr) || MASTER(cr))
            {
                restraint_->update(make_vec3<real>(r1[0],
                                                   r1[1],
                                                   r1[2]),
                                   make_vec3<real>(r2[0],
                                                   r2[1],
                                                   r2[2]),
                                   t);
            }
            // All ranks wait for the update to finish.
            // tMPI ranks are depending on structures that may have just been updated.
            if (cr != nullptr && DOMAINDECOMP(cr))
            {
                // Note: this assumes that all ranks are hitting this line, which is not generally true.
                // I need to find the right subcommunicator. What I really want is a _scoped_ communicator...
                gmx_barrier(cr);
            }

            auto result = restraint_->evaluate(make_vec3<real>(r1[0], r1[1], r1[2]),
                                               make_vec3<real>(r2[0], r2[1], r2[2]),
                                               t);

            int site1{static_cast<int>(sites_.front().index())};
            int aLocal{site1};
            // Set forces using index a if no domain decomposition, otherwise set with local index if available.
            if ((cr->dd == nullptr) || ga2la_get_home(cr->dd->ga2la,
                               site1,
                               &aLocal))
            {
                force[static_cast<size_t>(aLocal)][0] += result.force.x;
                force[static_cast<size_t>(aLocal)][1] += result.force.y;
                force[static_cast<size_t>(aLocal)][2] += result.force.z;
            }

            // Note: Currently calculateForces is called once per restraint and each restraint
            // applies to a pair of atoms. Future optimizations may consolidate multiple restraints
            // with possibly duplicated sites, in which case we may prefer to iterate over non-frozen
            // sites to apply forces without explicitly expressing pairwise symmetry as in the
            // following logic.
            int site2{static_cast<int>(sites_.back().index())};
            int bLocal{site2};
            if ((cr->dd == nullptr) || ga2la_get_home(cr->dd->ga2la,
                               site2,
                               &bLocal))
            {
                force[static_cast<size_t>(bLocal)][0] -= result.force.x;
                force[static_cast<size_t>(bLocal)][1] -= result.force.y;
                force[static_cast<size_t>(bLocal)][2] -= result.force.z;
            }

            if (int(t*1000) % 100 == 0)
            {
                if ((cr->dd == nullptr) || MASTER(cr))
                {
                    std::cout << "Evaluated restraint forces on sites at " << make_vec3<real>(r1[0], r1[1], r1[2]) << " and " << make_vec3<real>(r2[0], r2[1], r2[2]) << ": " << result.force << ". rank,time: " << cr->rank_pp_intranode << "," << t << "\n";
                }
            }

        }
    private:
        std::shared_ptr<gmx::IRestraintPotential> restraint_;
        std::vector<Site> sites_;
};

RestraintForceProvider::RestraintForceProvider(std::shared_ptr<gmx::IRestraintPotential> restraint,
                                               const std::vector<unsigned long int>& sites) :
    restraint_{std::move(restraint)},
    sites_{}
{
    assert(restraint_ != nullptr);
    assert(sites_.empty());
    for (auto&& site : sites)
    {
        sites_.emplace_back(site);
    }
    assert(sites_.size() >= 2);
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
                              const std::vector<unsigned long int>& sites);

        // Does default move work right (exception safe) with multiple unique_ptr members?
        RestraintMDModuleImpl(RestraintMDModuleImpl&&) noexcept = default;
        RestraintMDModuleImpl& operator=(RestraintMDModuleImpl&&) noexcept = default;

        ~RestraintMDModuleImpl() override;

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
