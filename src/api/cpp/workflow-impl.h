//
// Created by Eric Irrgang on 11/29/17.
//

#ifndef GROMACS_WORKFLOW_IMPL_H
#define GROMACS_WORKFLOW_IMPL_H

/*! \internal \file
 * \brief Implementation details for Workflow infrastructure.
 *
 * \ingroup gmxapi
 */

#include "gmxapi/exceptions.h"

namespace gmxapi
{

class WorkflowKeyError: public BasicException<WorkflowKeyError>
{
    public:
        using BasicException::BasicException;
};


class MDNodeSpecification : public NodeSpecification
{
    public:
        explicit MDNodeSpecification(std::string filename);

        std::unique_ptr<NodeSpecification> clone() override;

    private:
        std::string tprfilename_;
};


} // end namespace gmxapi

#endif //GROMACS_WORKFLOW_IMPL_H
