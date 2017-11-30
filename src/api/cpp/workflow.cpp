//
// Created by Eric Irrgang on 11/27/17.
//

#include "workflow.h"
#include "workflow-impl.h"
#include "gmxapi/exceptions.h"
#include "gmxapi/status.h"
#include "gromacs/compat/make_unique.h"
#include <cassert>

namespace gmxapi
{


NodeSpecification::~NodeSpecification() = default;


std::unique_ptr<NodeSpecification> MDNodeSpecification::clone()
{
    assert(!tprfilename_.empty());
    std::unique_ptr<NodeSpecification> node{nullptr};
    node = gmx::compat::make_unique<MDNodeSpecification>(tprfilename_);
    return node;
}

MDNodeSpecification::MDNodeSpecification(std::string filename) :
    tprfilename_{std::move(filename)}
{
    assert(!tprfilename_.empty());
}

Status Workflow::addNode(std::unique_ptr<NodeSpecification> &&spec) noexcept
{
    (void)spec;
    Status status{false};
    return status;
}

std::unique_ptr<Workflow> Workflow::create(const std::string &filename)
{
    std::string name{"MD"};
    auto spec = gmx::compat::make_unique<MDNodeSpecification>(filename);
    Workflow::Impl graph;
    graph.emplace(std::make_pair(name, std::move(spec)));
    auto workflow = gmx::compat::make_unique<Workflow>(std::move(graph));
    return workflow;
}

std::unique_ptr<NodeSpecification> Workflow::getNode(const NodeKey &key) const noexcept
{
    const Impl &graph = graph_;
    assert((graph.count(key) == 0) || (graph.count(key) == 1));
    auto const iter = graph.find(key);
    // Can return a null NodeSpecification if key is not found...
//    if (iter == graph.end())
//    {
//        auto const error = WorkflowKeyError(std::move(std::string(key)));
//        throw error;
//    }
    std::unique_ptr<NodeSpecification> node{nullptr};
    if (iter != graph.end())
    {
        node = iter->second->clone();
    }
    return node;
}

Workflow::Workflow(Workflow::Impl &&impl) :
    graph_{std::forward<Workflow::Impl>(impl)}
{}


} // end namespace gmxapi
