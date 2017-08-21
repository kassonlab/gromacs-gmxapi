#ifndef GMXAPI_RUNNER_IMPL_H
#define GMXAPI_RUNNER_IMPL_H

#include <cstdio>

#include <functional>
#include <iostream>
#include <memory>
#include <string>

#include <gromacs/commandline/filenm.h>
#include <gromacs/ewald/pme-internal.h>
#include <gromacs/hardware/hw_info.h>
#include <gromacs/mdlib/integrator.h>
#include <gromacs/mdrunutility/mdmodules.h>
#include "gromacs/timing/wallcycle.h"
#include <gromacs/utility/basedefinitions.h>
#include <gromacs/utility/loggerbuilder.h>
#include <gromacs/utility/real.h>

struct t_commrec;
struct gmx_hw_info_t;
struct gmx_pme_t;
struct t_mdatoms;
struct gmx_vsite_t;
struct t_fcdata;
struct ReplicaExchangeParameters;
struct t_nrnb;

namespace gmx
{
class IntegratorParams;
struct DomDecParams;
}

namespace gmxapi
{
class MDInput;

/*!
 * \brief Reimplementation of gmx::mdrunner while requirements are sorted out.
 *
 * Todo: migrate to gmx::MDRunner class.
 *
 * Objects of this class are instantiated as execution is launched using input previously provided at higher API levels.
 *
 * This implementation is basically a rough copy-paste to allow multiple run() calls during a single execution. The
 * constructor includes a lot of options-processing code that can be migrated out as library infrastructure evolves.
 */
class RunnerImpl

{
    public:
        // Default initializations are delegated to the default constructor, but a better practice is probably to
        // delegate the default constructor to the most complete constructor.
        RunnerImpl();

        // Prepare the simulation and advance to the point right before calling the integrator.
        explicit RunnerImpl(const std::string &filename);

        // Call the integrator with current parameters
        int run();

        // Set number of steps and call integrator
        int run(unsigned int numSteps);

        // Tear down routine
        int close();

        // Clean up a la mdrunner() and gmx_mdrun.
        ~RunnerImpl();

        // Get a copy of the current positions in the local state structure, which may be null.
        std::shared_ptr < std::vector < std::array<real, 3>>> getX() const;

    private:
        unsigned int long                          flags_;
        std::unique_ptr<MDInput>                   input_;
        FILE                                     * fplog_;
        std::unique_ptr<ReplicaExchangeParameters> replExParams_;
        std::unique_ptr<gmx::LoggerOwner>          logOwner_;
        std::unique_ptr<gmx_hw_info_t>             hardwareInfo_;
        std::unique_ptr < gmx_pme_t, std::function < void(gmx_pme_t*)>> pmeData_;
        t_commrec                                 *commRec_; // Note unmanaged opaque pointer
        gmx_bool                                   doAppendFiles_;
        bool                                doMembed_;
        real                                ewaldcoeff_q_;
        real                                ewaldcoeff_lj_;
        gmx::MDModules                      mdModules_;
        gmx_output_env_t                   *oenv_; // hidden definition (opaque pointer) with no deleter
        bool                                verbose_;
        int                                 nstglobalcomm_;
        std::unique_ptr<gmx::DomDecParams>  ddParams_;
        std::unique_ptr<ObservablesHistory> observablesHistory_;
        int                                 imdport_;
        int                                 nstepout_;
        real                                cpt_period_;
        std::unique_ptr<t_mdatoms>          mdAtoms_;
        std::unique_ptr<gmx_vsite_t>        vSite_;
        std::unique_ptr<t_fcdata>           forceCalcData_;
        std::unique_ptr<t_nrnb>             nrNonBonded_;
        gmx_wallcycle_t                     wallCycle_;
        std::unique_ptr<t_forcerec>         forceRecord_;
        std::unique_ptr < gmx_membed_t, std::function < void(gmx_membed_t*)>> membed_;
        real                                maxHours_;
        gmx_walltime_accounting_t           walltimeAccounting_; // Note this is actually an opaque pointer
        bool                                initialized_;

        std::vector<t_filenm>               fnm_ {{
                                                      { efTPR, nullptr,      nullptr,       ffREAD },
                                                      { efTRN, "-o",      nullptr,       ffWRITE },
                                                      { efCOMPRESSED, "-x", nullptr,     ffOPTWR },
                                                      { efCPT, "-cpi",    nullptr,       ffOPTRD | ffALLOW_MISSING },
                                                      { efCPT, "-cpo",    nullptr,       ffOPTWR },
                                                      { efSTO, "-c",      "confout",  ffWRITE },
                                                      { efEDR, "-e",      "ener",     ffWRITE },
                                                      { efLOG, "-g",      "md",       ffWRITE },
                                                      { efXVG, "-dhdl",   "dhdl",     ffOPTWR },
                                                      { efXVG, "-field",  "field",    ffOPTWR },
                                                      { efXVG, "-table",  "table",    ffOPTRD },
                                                      { efXVG, "-tablep", "tablep",   ffOPTRD },
                                                      { efXVG, "-tableb", "table",    ffOPTRDMULT },
                                                      { efTRX, "-rerun",  "rerun",    ffOPTRD },
                                                      { efXVG, "-tpi",    "tpi",      ffOPTWR },
                                                      { efXVG, "-tpid",   "tpidist",  ffOPTWR },
                                                      { efEDI, "-ei",     "sam",      ffOPTRD },
                                                      { efXVG, "-eo",     "edsam",    ffOPTWR },
                                                      { efXVG, "-devout", "deviatie", ffOPTWR },
                                                      { efXVG, "-runav",  "runaver",  ffOPTWR },
                                                      { efXVG, "-px",     "pullx",    ffOPTWR },
                                                      { efXVG, "-pf",     "pullf",    ffOPTWR },
                                                      { efXVG, "-ro",     "rotation", ffOPTWR },
                                                      { efLOG, "-ra",     "rotangles", ffOPTWR },
                                                      { efLOG, "-rs",     "rotslabs", ffOPTWR },
                                                      { efLOG, "-rt",     "rottorque", ffOPTWR },
                                                      { efMTX, "-mtx",    "nm",       ffOPTWR },
                                                      { efRND, "-multidir", nullptr,      ffOPTRDMULT},
                                                      { efDAT, "-membed", "membed",   ffOPTRD },
                                                      { efTOP, "-mp",     "membed",   ffOPTRD },
                                                      { efNDX, "-mn",     "membed",   ffOPTRD },
                                                      { efXVG, "-if",     "imdforces", ffOPTWR },
                                                      { efXVG, "-swap",   "swapions", ffOPTWR }
                                                  }};
};


}      // end namesace gmxapi

#endif // header guard
