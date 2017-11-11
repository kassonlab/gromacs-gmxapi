//
// Created by Eric Irrgang on 11/10/17.
//

#ifndef GROMACS_RESTRAINTMDMODULE_IMPL_H
#define GROMACS_RESTRAINTMDMODULE_IMPL_H

#include <iostream>
#include "gromacs/mdtypes/imdpoptionprovider.h"
#include "gromacs/mdtypes/imdoutputprovider.h"
#include "gromacs/mdtypes/iforceprovider.h"
#include "gromacs/mdtypes/imdmodule.h"
#include "gromacs/utility/arrayref.h"

namespace gmx
{


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
        RestraintForceProvider() = delete;
        explicit RestraintForceProvider(std::shared_ptr<gmx::IRestraintPotential> restraint);

        void calculateForces(const t_commrec          *cr,
                             const t_mdatoms          *mdatoms,
                             const matrix              box,
                             double                    t,
                             const rvec               *x,
                             gmx::ArrayRef<gmx::RVec>  force)

        override
        {
            assert(restraint_ != nullptr);
            size_t a = 1;
            size_t b = 3;

            const gmx::detail::vec3<real> r1{gmx::detail::make_vec3<real>(x[a][0], x[a][1], x[a][2])};
            const gmx::detail::vec3<real> r2{gmx::detail::make_vec3<real>(x[b][0], x[b][1], x[b][2])};

            auto result = restraint_->evaluate(r1, r2, t);

            force[a][0] += result.force.x;
            force[a][1] += result.force.y;
            force[a][2] += result.force.z;

            force[b][0] -= result.force.x;
            force[b][1] -= result.force.y;
            force[b][2] -= result.force.z;

            std::cout << "Evaluated restraint forces on atoms at " << r1 << " and " << r2 << ": " << result.force << "\n";

        }
    private:
        std::shared_ptr<gmx::IRestraintPotential> restraint_;
};

RestraintForceProvider::RestraintForceProvider(std::shared_ptr<gmx::IRestraintPotential> restraint) :
    restraint_{std::move(restraint)}
{
    assert(restraint_ != nullptr);
}

class RestraintMDModuleImpl final: public gmx::IMDModule
{
    public:
        RestraintMDModuleImpl() = delete;
        explicit RestraintMDModuleImpl(std::shared_ptr<gmx::IRestraintPotential>);

        RestraintMDModuleImpl(RestraintMDModuleImpl&&) = default;
        RestraintMDModuleImpl& operator=(RestraintMDModuleImpl&&) = default;

        virtual ~RestraintMDModuleImpl();

        IMdpOptionProvider *mdpOptionProvider() override;

        IMDOutputProvider *outputProvider() override;

        /*!
         * \brief
         *
         * \param forceProviders force module manager in from the force record that will call this.
         *
         * The calling code must ensure that this object stays alive as long as forceProviders needs
         * this. Typically that is the duration of a do_md() call.
         */
        void initForceProviders(ForceProviders *forceProviders) override;

        std::unique_ptr<RestraintForceProvider> forceProvider_;
        std::unique_ptr<RestraintOutputProvider> outputProvider_;
        std::unique_ptr<RestraintOptionProvider> optionProvider_;
};

} // end namespace gmx

#endif //GROMACS_RESTRAINTMDMODULE_IMPL_H
