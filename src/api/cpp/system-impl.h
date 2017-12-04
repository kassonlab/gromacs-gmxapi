#ifndef GMXAPI_SYSTEM_IMPL_H
#define GMXAPI_SYSTEM_IMPL_H

#include "gmxapi/system.h"

namespace gmxapi
{

class Workflow;

class System::Impl final
{
    public:
        Impl();
        ~Impl();

        Impl(Impl&&) noexcept = default;
        Impl& operator=(Impl&&) noexcept = default;

        explicit Impl(std::unique_ptr<gmxapi::Workflow>&& workflow) noexcept;

        Status status() const;

        std::shared_ptr<Session> launch(std::shared_ptr<Context> context);
        std::shared_ptr<Session> launch();

        Status setRestraint(std::shared_ptr<gmxapi::MDModule> module);
        std::shared_ptr<MDWorkSpec> getSpec();

    private:
        std::shared_ptr<Context> context_;
        std::shared_ptr<Workflow> workflow_;
        // \todo merge Workflow and MDWorkSpec
        std::shared_ptr<gmxapi::MDWorkSpec> spec_;
        std::unique_ptr<Status> status_;
};

}      // end namespace gmxapi

#endif // header guard
