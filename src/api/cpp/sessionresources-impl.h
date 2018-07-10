//
// Created by Eric Irrgang on 7/10/18.
//

#ifndef GMXAPI_SESSION_RESOURCES_IMPL_H
#define GMXAPI_SESSION_RESOURCES_IMPL_H

#include <string>

#include "gmxapi/md/mdsignals.h"

namespace gmxapi {

class SessionResources
{
    public:
        /*!
         * \brief Construct a resources object for the named operation.
         *
         * \param session implementation object backing these resources.
         * \param name Unique name of workflow operation.
         */
        SessionResources(SessionImpl* session, std::string name);
        SessionResources() = delete;

        // Objects of this type should only exist in their Session container.
        // If necessary, ownership can be transferred by owning through a unique_ptr handle.
        SessionResources(const SessionResources&) = delete;
        SessionResources& operator=(const SessionResources&) = delete;
        SessionResources(SessionResources&&) = delete;
        SessionResources& operator=(SessionResources&&) = delete;

        ~SessionResources();

        /*!
         * \brief Get the name of the gmxapi operation for which these resources exist.
         * \return workflow element name
         */
        const std::string name() const ;

        Signal getMdrunnerSignal(md::signals signal);
    private:
        /*!
         * \brief pointer to the session owning these resources
         */
        SessionImpl* sessionImpl_{nullptr};

        /*!
         * \brief name of the associated gmxapi operation
         */
        std::string name_;
};

} // end namespace gmxapi

#endif //GMXAPI_SESSION_RESOURCES_IMPL_H
