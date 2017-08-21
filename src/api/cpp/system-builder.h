#ifndef GMXAPI_SYSTEM_BUILDER_H
#define GMXAPI_SYSTEM_BUILDER_H

#include "gmxapi/system.h"
//#include "gmxapi/md.h"

namespace gmxapi
{

class MDInput;

class System::Builder
{
    public:
        // Construct a trivial builder that just uses a TPR file.
        explicit Builder();

        // Allow an appropriate default Context to be determined and configured.
        Builder &defaultContext(const MDInput &inputrec);

        // Use the information in the input record to configure an appropriate runner.
        Builder &runner(const MDInput &inputrec);

        Builder &structure(const MDInput &inputrec);

        Builder &mdEngine(const MDInput &inputrec);

        Builder &topology(const MDInput &inputrec);

        // Pass ownership of the assembled System.
        std::unique_ptr<System> build();
    private:
        std::unique_ptr<System> system_;
};

}      // end namespace gmxapi

#endif // header guard
