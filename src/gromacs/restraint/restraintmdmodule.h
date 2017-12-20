//
// Created by Eric Irrgang on 11/10/17.
//

#ifndef GROMACS_RESTRAINTMDMODULE_H
#define GROMACS_RESTRAINTMDMODULE_H

/*! \libinternal \file
 * \brief Library interface for RestraintMDModule
 * 
 * \ingroup module_restraint
 */

#include "gromacs/mdtypes/imdmodule.h"
#include "restraintpotential.h"

namespace gmx
{

// Forward declaration to allow opaque pointer to library internal class.
class RestraintMDModuleImpl;

/*!
 * \brief MDModule wrapper for Restraint implementations.
 *
 * Shares ownership of an object implementing the IRestraintPotential interface.
 * Provides the IMDModule interfaces.
 */
class RestraintMDModule final : public gmx::IMDModule
{
    public:
        RestraintMDModule() = delete;

        ~RestraintMDModule() override ;

        /*!
         * \brief Wrap a restraint potential as an MDModule
         *
         * \param restraint shared ownership of an object for calculating restraint forces.
         * \return new wrapper object sharing ownership of restraint.
         *
         * Consumers of the interfaces provided by an IMDModule do not extend the lifetime
         * of the interface objects returned by mdpOptionProvider(), outputProvider(), or
         * registered via initForceProviders(). Calling code must keep this object alive
         * as long as those interfaces are needed (probably the duration of an MD run).
         */
        static std::unique_ptr<RestraintMDModule>
        create(std::shared_ptr<gmx::IRestraintPotential> restraint, unsigned long int site1, unsigned long int site2);

        /*!
         * \brief Implement IMDModule interface
         *
         * \return see gmx::IMDModule::outputProvider.
         */
        IMdpOptionProvider *mdpOptionProvider() override;

        /*!
         * \brief Implement IMDModule interface
         *
         * \return see gmx::IMDModule::outputProvider.
         */
        IMDOutputProvider *outputProvider() override;

        /*!
         * \brief Implement IMDModule interface.
         *
         * See gmx::IMDModule::initForceProviders()
         * \param forceProviders manager in the force record.
         */
        void initForceProviders(ForceProviders *forceProviders) override;

    private:
        /*!
         * \brief Private constructor used by static create() method.
         */
        explicit RestraintMDModule(std::unique_ptr<RestraintMDModuleImpl> restraint);
        /*!
         * \brief Private implementation opaque pointer.
         */
        std::unique_ptr<RestraintMDModuleImpl> impl_;
};

} // end namespace gmx

#endif //GROMACS_RESTRAINTMDMODULE_H
