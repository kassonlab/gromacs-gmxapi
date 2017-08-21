#include "gmxpre.h"

#include "config.h"

#include <assert.h>
#include <signal.h>
#include <stdlib.h>
#include <string.h>

#include <algorithm>

#include "gmxapi/runner.h"

#include "api/cpp/md-impl.h"
#include "api/cpp/runner-impl.h"
#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/ewald/pme.h"
#include "gromacs/ewald/pme-internal.h"
#include "gromacs/fileio/checkpoint.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/hardware/cpuinfo.h"
#include "gromacs/hardware/detecthardware.h"
#include "gromacs/hardware/hardwareassign.h"
#include "gromacs/hardware/hw_info.h"
#include "gromacs/listed-forces/disre.h"
#include "gromacs/listed-forces/orires.h"
#include "gromacs/math/calculate-ewald-splitting-coefficient.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/calc_verletbuf.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/mdlib/force.h"
#include "gromacs/mdlib/forcerec.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/integrator.h"
#include "gromacs/mdlib/main.h"
#include "gromacs/mdlib/md_support.h"
#include "gromacs/mdlib/mdatoms.h"
#include "gromacs/mdlib/mdrun.h"
#include "gromacs/mdlib/minimize.h"
#include "gromacs/mdlib/nbnxn_search.h"
#include "gromacs/mdlib/qmmm.h"
#include "gromacs/mdlib/sighandler.h"
#include "gromacs/mdlib/sim_util.h"
#include "gromacs/mdlib/tpi.h"
#include "gromacs/mdrunutility/handlerestart.h"
#include "gromacs/mdrunutility/mdmodules.h"
#include "gromacs/mdrunutility/threadaffinity.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/edsamhistory.h"
#include "gromacs/mdtypes/energyhistory.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/observableshistory.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/mdtypes/swaphistory.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pulling/pull.h"
#include "gromacs/pulling/pull_rotation.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/filestream.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/init.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/loggerbuilder.h"
#include "gromacs/utility/pleasecite.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"
#include "programs/mdrun/deform.h"
#include "programs/mdrun/membed.h"
#include "programs/mdrun/repl_ex.h"
#include "programs/mdrun/resource-division.h"

#ifdef GMX_FAHCORE
#include "programs/mdrun/corewrap.h"
#endif

#include "gromacs/compat/make_unique.h"
#include "gromacs/mdlib/integrator.h"
using gmx::compat::make_unique;

namespace gmx
{

//! MD simulations
integrator_t do_md;

}      // namespace gmx

namespace gmxapi
{
using gmx::integrator_t;
using gmx::do_md;
using gmx::do_steep;
using gmx::do_cg;
using gmx::do_nm;
using gmx::do_lbfgs;
using gmx::do_tpi;
using gmx::IntegratorParams;
using gmx::APIError;
using gmx::NotImplementedError;

/*! \brief Return whether either of the command-line parameters that
 *  will trigger a multi-simulation is set */
static inline bool is_multisim_option_set(int argc, const char *const argv[])
{
    for (int i = 0; i < argc; ++i)
    {
        if (strcmp(argv[i], "-multi") == 0 || strcmp(argv[i], "-multidir") == 0)
        {
            return true;
        }
    }
    return false;
}

int mdrunner(gmx_hw_opt_t *hw_opt,
             FILE *fplog, t_commrec *cr, int nfile,
             const t_filenm fnm[], const gmx_output_env_t *oenv, gmx_bool bVerbose,
             int nstglobalcomm,
             ivec ddxyz, int dd_rank_order, int npme, real rdd, real rconstr,
             const char *dddlb_opt, real dlb_scale,
             const char *ddcsx, const char *ddcsy, const char *ddcsz,
             const char *nbpu_opt, int nstlist_cmdline,
             gmx_int64_t nsteps_cmdline, int nstepout, int resetstep,
             int gmx_unused nmultisim,
             const ReplicaExchangeParameters &replExParams,
             real pforce, real cpt_period, real max_hours,
             int imdport, unsigned long Flags);

struct mdrunner_arglist
{
    /// \cond
    gmx_hw_opt_t                     hw_opt;
    FILE                            *fplog;
    t_commrec                       *cr;
    int                              nfile;
    const t_filenm                  *fnm;
    const gmx_output_env_t          *oenv;
    gmx_bool                         bVerbose;
    int                              nstglobalcomm;
    ivec                             ddxyz;
    int                              dd_rank_order;
    int                              npme;
    real                             rdd;
    real                             rconstr;
    const char                      *dddlb_opt;
    real                             dlb_scale;
    const char                      *ddcsx;
    const char                      *ddcsy;
    const char                      *ddcsz;
    const char                      *nbpu_opt;
    int                              nstlist_cmdline;
    gmx_int64_t                      nsteps_cmdline;
    int                              nstepout;
    int                              resetstep;
    int                              nmultisim;
    const ReplicaExchangeParameters *replExParams;
    real                             pforce;
    real                             cpt_period;
    real                             max_hours;
    int                              imdport;
    unsigned long                    Flags;
    /// \endcond
};

mdrunner_arglist make_mdrunner_arglist(gmx_hw_opt_t *hw_opt,
                                       FILE *fplog, t_commrec *cr, int nfile,
                                       const t_filenm fnm[], const gmx_output_env_t *oenv, gmx_bool bVerbose,
                                       int nstglobalcomm, ivec ddxyz, int dd_rank_order, int npme,
                                       real rdd, real rconstr, const char *dddlb_opt, real dlb_scale,
                                       const char *ddcsx, const char *ddcsy, const char *ddcsz,
                                       const char *nbpu_opt, int nstlist_cmdline,
                                       gmx_int64_t nsteps_cmdline, int nstepout, int resetstep,
                                       int nmultisim,
                                       const ReplicaExchangeParameters &replExParams,
                                       real pforce, real cpt_period, real max_hours,
                                       int imdport, unsigned long Flags);
// int mdrunner(mdrunner_arglist args);

//! First step used in pressure scaling
gmx_int64_t         deform_init_init_step_tpx;
//! Initial box for pressure scaling
matrix              deform_init_box_tpx;
//! MPI variable for use in pressure scaling
tMPI_Thread_mutex_t deform_init_box_mutex = TMPI_THREAD_MUTEX_INITIALIZER;

#if GMX_THREAD_MPI
/* The minimum number of atoms per tMPI thread. With fewer atoms than this,
 * the number of threads will get lowered.
 */
#define MIN_ATOMS_PER_MPI_THREAD    90
#define MIN_ATOMS_PER_GPU           900

/* The function used for spawning threads. Extracts the mdrunner()
   arguments from its one argument and calls mdrunner(), after making
   a commrec. */
// called by mdrunner_start_threads
static void mdrunner_start_fn(void *arg)
{
    try
    {
        struct mdrunner_arglist *mda = (struct mdrunner_arglist*)arg;
        struct mdrunner_arglist  mc  = *mda; /* copy the arg list to make sure
                                                that it's thread-local. This doesn't
                                                copy pointed-to items, of course,
                                                but those are all const. */
        t_commrec *cr;                       /* we need a local version of this */
        FILE      *fplog = nullptr;
        t_filenm  *fnm;

        fnm = dup_tfn(mc.nfile, mc.fnm);

        cr = reinitialize_commrec_for_this_thread(mc.cr);

        if (MASTER(cr))
        {
            fplog = mc.fplog;
        }

        mdrunner(&mc.hw_opt, fplog, cr, mc.nfile, fnm, mc.oenv,
                 mc.bVerbose, mc.nstglobalcomm,
                 mc.ddxyz, mc.dd_rank_order, mc.npme, mc.rdd,
                 mc.rconstr, mc.dddlb_opt, mc.dlb_scale,
                 mc.ddcsx, mc.ddcsy, mc.ddcsz,
                 mc.nbpu_opt, mc.nstlist_cmdline,
                 mc.nsteps_cmdline, mc.nstepout, mc.resetstep,
                 mc.nmultisim, *mc.replExParams, mc.pforce,
                 mc.cpt_period, mc.max_hours, mc.imdport, mc.Flags);
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
}


/* called by mdrunner() to start a specific number of threads (including
   the main thread) for thread-parallel runs. This in turn calls mdrunner()
   for each thread.
   All options besides nthreads are the same as for mdrunner(). */
static t_commrec *mdrunner_start_threads(gmx_hw_opt_t *hw_opt,
                                         FILE *fplog, t_commrec *cr, int nfile,
                                         const t_filenm fnm[], const gmx_output_env_t *oenv, gmx_bool bVerbose,
                                         int nstglobalcomm,
                                         ivec ddxyz, int dd_rank_order, int npme,
                                         real rdd, real rconstr,
                                         const char *dddlb_opt, real dlb_scale,
                                         const char *ddcsx, const char *ddcsy, const char *ddcsz,
                                         const char *nbpu_opt, int nstlist_cmdline,
                                         gmx_int64_t nsteps_cmdline,
                                         int nstepout, int resetstep,
                                         int nmultisim,
                                         const ReplicaExchangeParameters &replExParams,
                                         real pforce, real cpt_period, real max_hours,
                                         unsigned long Flags)
{
    int                      ret;
    struct mdrunner_arglist *mda;
    t_commrec               *crn; /* the new commrec */
    t_filenm                *fnmn;

    /* first check whether we even need to start tMPI */
    if (hw_opt->nthreads_tmpi < 2)
    {
        return cr;
    }

    /* a few small, one-time, almost unavoidable memory leaks: */
    snew(mda, 1);
    fnmn = dup_tfn(nfile, fnm);

    /* fill the data structure to pass as void pointer to thread start fn */
    /* hw_opt contains pointers, which should all be NULL at this stage */
    mda->hw_opt          = *hw_opt;
    mda->fplog           = fplog;
    mda->cr              = cr;
    mda->nfile           = nfile;
    mda->fnm             = fnmn;
    mda->oenv            = oenv;
    mda->bVerbose        = bVerbose;
    mda->nstglobalcomm   = nstglobalcomm;
    mda->ddxyz[XX]       = ddxyz[XX];
    mda->ddxyz[YY]       = ddxyz[YY];
    mda->ddxyz[ZZ]       = ddxyz[ZZ];
    mda->dd_rank_order   = dd_rank_order;
    mda->npme            = npme;
    mda->rdd             = rdd;
    mda->rconstr         = rconstr;
    mda->dddlb_opt       = dddlb_opt;
    mda->dlb_scale       = dlb_scale;
    mda->ddcsx           = ddcsx;
    mda->ddcsy           = ddcsy;
    mda->ddcsz           = ddcsz;
    mda->nbpu_opt        = nbpu_opt;
    mda->nstlist_cmdline = nstlist_cmdline;
    mda->nsteps_cmdline  = nsteps_cmdline;
    mda->nstepout        = nstepout;
    mda->resetstep       = resetstep;
    mda->nmultisim       = nmultisim;
    mda->replExParams    = &replExParams;
    mda->pforce          = pforce;
    mda->cpt_period      = cpt_period;
    mda->max_hours       = max_hours;
    mda->Flags           = Flags;

    /* now spawn new threads that start mdrunner_start_fn(), while
       the main thread returns, we set thread affinity later */
    ret = tMPI_Init_fn(TRUE, hw_opt->nthreads_tmpi, TMPI_AFFINITY_NONE,
                       mdrunner_start_fn, (void*)(mda) );
    if (ret != TMPI_SUCCESS)
    {
        return nullptr;
    }

    crn = reinitialize_commrec_for_this_thread(cr);
    return crn;
}

#endif /* GMX_THREAD_MPI */


/*! \brief Cost of non-bonded kernels
 *
 * We determine the extra cost of the non-bonded kernels compared to
 * a reference nstlist value of 10 (which is the default in grompp).
 */
static const int    nbnxnReferenceNstlist = 10;
//! The values to try when switching
const int           nstlist_try[] = { 20, 25, 40 };
//! Number of elements in the neighborsearch list trials.
#define NNSTL  sizeof(nstlist_try)/sizeof(nstlist_try[0])
/* Increase nstlist until the non-bonded cost increases more than listfac_ok,
 * but never more than listfac_max.
 * A standard (protein+)water system at 300K with PME ewald_rtol=1e-5
 * needs 1.28 at rcoulomb=0.9 and 1.24 at rcoulomb=1.0 to get to nstlist=40.
 * Note that both CPU and GPU factors are conservative. Performance should
 * not go down due to this tuning, except with a relatively slow GPU.
 * On the other hand, at medium/high parallelization or with fast GPUs
 * nstlist will not be increased enough to reach optimal performance.
 */
/* CPU: pair-search is about a factor 1.5 slower than the non-bonded kernel */
//! Max OK performance ratio beween force calc and neighbor searching
static const float  nbnxn_cpu_listfac_ok    = 1.05;
//! Too high performance ratio beween force calc and neighbor searching
static const float  nbnxn_cpu_listfac_max   = 1.09;
/* CPU: pair-search is about a factor 2-3 slower than the non-bonded kernel */
//! Max OK performance ratio beween force calc and neighbor searching
static const float  nbnxn_knl_listfac_ok    = 1.22;
//! Too high performance ratio beween force calc and neighbor searching
static const float  nbnxn_knl_listfac_max   = 1.3;
/* GPU: pair-search is a factor 1.5-3 slower than the non-bonded kernel */
//! Max OK performance ratio beween force calc and neighbor searching
static const float  nbnxn_gpu_listfac_ok    = 1.20;
//! Too high performance ratio beween force calc and neighbor searching
static const float  nbnxn_gpu_listfac_max   = 1.30;

/*! \brief Try to increase nstlist when using the Verlet cut-off scheme */
// ir->nstlist and ir->rlist are set during this function call.
static void increase_nstlist(FILE *fp, t_commrec *cr,
                             t_inputrec *ir, int nstlist_cmdline,
                             const gmx_mtop_t *mtop, matrix box,
                             gmx_bool bGPU, const gmx::CpuInfo &cpuinfo)
{
    float                  listfac_ok, listfac_max;
    int                    nstlist_orig, nstlist_prev;
    verletbuf_list_setup_t ls;
    real                   rlistWithReferenceNstlist, rlist_inc, rlist_ok, rlist_max;
    real                   rlist_new, rlist_prev;
    size_t                 nstlist_ind = 0;
    gmx_bool               bBox, bDD, bCont;
    const char            *nstl_gpu = "\nFor optimal performance with a GPU nstlist (now %d) should be larger.\nThe optimum depends on your CPU and GPU resources.\nYou might want to try several nstlist values.\n";
    const char            *nve_err  = "Can not increase nstlist because an NVE ensemble is used";
    const char            *vbd_err  = "Can not increase nstlist because verlet-buffer-tolerance is not set or used";
    const char            *box_err  = "Can not increase nstlist because the box is too small";
    const char            *dd_err   = "Can not increase nstlist because of domain decomposition limitations";
    char                   buf[STRLEN];

    if (nstlist_cmdline <= 0)
    {
        if (ir->nstlist == 1)
        {
            /* The user probably set nstlist=1 for a reason,
             * don't mess with the settings.
             */
            return;
        }

        if (fp != nullptr && bGPU && ir->nstlist < nstlist_try[0])
        {
            fprintf(fp, nstl_gpu, ir->nstlist);
        }
        nstlist_ind = 0;
        while (nstlist_ind < NNSTL && ir->nstlist >= nstlist_try[nstlist_ind])
        {
            nstlist_ind++;
        }
        if (nstlist_ind == NNSTL)
        {
            /* There are no larger nstlist value to try */
            return;
        }
    }

    if (EI_MD(ir->eI) && ir->etc == etcNO)
    {
        if (MASTER(cr))
        {
            fprintf(stderr, "%s\n", nve_err);
        }
        if (fp != nullptr)
        {
            fprintf(fp, "%s\n", nve_err);
        }

        return;
    }

    if (ir->verletbuf_tol == 0 && bGPU)
    {
        gmx_fatal(FARGS, "You are using an old tpr file with a GPU, please generate a new tpr file with an up to date version of grompp");
    }

    if (ir->verletbuf_tol < 0)
    {
        if (MASTER(cr))
        {
            fprintf(stderr, "%s\n", vbd_err);
        }
        if (fp != nullptr)
        {
            fprintf(fp, "%s\n", vbd_err);
        }

        return;
    }

    if (bGPU)
    {
        listfac_ok  = nbnxn_gpu_listfac_ok;
        listfac_max = nbnxn_gpu_listfac_max;
    }
    else if (cpuinfo.feature(gmx::CpuInfo::Feature::X86_Avx512ER))
    {
        listfac_ok  = nbnxn_knl_listfac_ok;
        listfac_max = nbnxn_knl_listfac_max;
    }
    else
    {
        listfac_ok  = nbnxn_cpu_listfac_ok;
        listfac_max = nbnxn_cpu_listfac_max;
    }

    nstlist_orig = ir->nstlist;
    if (nstlist_cmdline > 0)
    {
        if (fp)
        {
            sprintf(buf, "Getting nstlist=%d from command line option",
                    nstlist_cmdline);
        }
        //
        // This is where *ir is an output param...
        //
        ir->nstlist = nstlist_cmdline;
    }

    verletbuf_get_list_setup(TRUE, bGPU, &ls);

    /* Allow rlist to make the list a given factor larger than the list
     * would be with the reference value for nstlist (10).
     */
    nstlist_prev = ir->nstlist;
    ir->nstlist  = nbnxnReferenceNstlist;
    calc_verlet_buffer_size(mtop, det(box), ir, -1, &ls, nullptr,
                            &rlistWithReferenceNstlist);
    ir->nstlist  = nstlist_prev;

    /* Determine the pair list size increase due to zero interactions */
    rlist_inc = nbnxn_get_rlist_effective_inc(ls.cluster_size_j,
                                              mtop->natoms/det(box));
    rlist_ok  = (rlistWithReferenceNstlist + rlist_inc)*std::cbrt(listfac_ok) - rlist_inc;
    rlist_max = (rlistWithReferenceNstlist + rlist_inc)*std::cbrt(listfac_max) - rlist_inc;
    if (debug)
    {
        fprintf(debug, "nstlist tuning: rlist_inc %.3f rlist_ok %.3f rlist_max %.3f\n",
                rlist_inc, rlist_ok, rlist_max);
    }

    nstlist_prev = nstlist_orig;
    rlist_prev   = ir->rlist;
    do
    {
        if (nstlist_cmdline <= 0)
        {
            ir->nstlist = nstlist_try[nstlist_ind];
        }

        /* Set the pair-list buffer size in ir */
        calc_verlet_buffer_size(mtop, det(box), ir, -1, &ls, nullptr, &rlist_new);

        /* Does rlist fit in the box? */
        bBox = (gmx::square(rlist_new) < max_cutoff2(ir->ePBC, box));
        bDD  = TRUE;
        if (bBox && DOMAINDECOMP(cr))
        {
            /* Check if rlist fits in the domain decomposition */
            if (inputrec2nboundeddim(ir) < DIM)
            {
                gmx_incons("Changing nstlist with domain decomposition and unbounded dimensions is not implemented yet");
            }
            t_state state_tmp;
            copy_mat(box, state_tmp.box);
            bDD = change_dd_cutoff(cr, &state_tmp, ir, rlist_new);
        }

        if (debug)
        {
            fprintf(debug, "nstlist %d rlist %.3f bBox %d bDD %d\n",
                    ir->nstlist, rlist_new, bBox, bDD);
        }

        bCont = FALSE;

        if (nstlist_cmdline <= 0)
        {
            if (bBox && bDD && rlist_new <= rlist_max)
            {
                /* Increase nstlist */
                nstlist_prev = ir->nstlist;
                rlist_prev   = rlist_new;
                bCont        = (nstlist_ind+1 < NNSTL && rlist_new < rlist_ok);
            }
            else
            {
                /* Stick with the previous nstlist */
                ir->nstlist = nstlist_prev;
                rlist_new   = rlist_prev;
                bBox        = TRUE;
                bDD         = TRUE;
            }
        }

        nstlist_ind++;
    }
    while (bCont);

    if (!bBox || !bDD)
    {
        gmx_warning(!bBox ? box_err : dd_err);
        if (fp != nullptr)
        {
            fprintf(fp, "\n%s\n", bBox ? box_err : dd_err);
        }
        ir->nstlist = nstlist_orig;
    }
    else if (ir->nstlist != nstlist_orig || rlist_new != ir->rlist)
    {
        sprintf(buf, "Changing nstlist from %d to %d, rlist from %g to %g",
                nstlist_orig, ir->nstlist,
                ir->rlist, rlist_new);
        if (MASTER(cr))
        {
            fprintf(stderr, "%s\n\n", buf);
        }
        if (fp != nullptr)
        {
            fprintf(fp, "%s\n\n", buf);
        }
        ir->rlist     = rlist_new;
    }
}

/*! \brief Initialize variables for Verlet scheme simulation */
static void prepare_verlet_scheme(FILE                           *fplog,
                                  t_commrec                      *cr,
                                  t_inputrec                     *ir,
                                  int                             nstlist_cmdline,
                                  const gmx_mtop_t               *mtop,
                                  matrix                          box,
                                  gmx_bool                        bUseGPU,
                                  const gmx::CpuInfo             &cpuinfo)
{
    /* For NVE simulations, we will retain the initial list buffer */
    if (EI_DYNAMICS(ir->eI) &&
        ir->verletbuf_tol > 0 &&
        !(EI_MD(ir->eI) && ir->etc == etcNO))
    {
        /* Update the Verlet buffer size for the current run setup */
        verletbuf_list_setup_t ls;
        real                   rlist_new;

        /* Here we assume SIMD-enabled kernels are being used. But as currently
         * calc_verlet_buffer_size gives the same results for 4x8 and 4x4
         * and 4x2 gives a larger buffer than 4x4, this is ok.
         */
        verletbuf_get_list_setup(TRUE, bUseGPU, &ls);

        calc_verlet_buffer_size(mtop, det(box), ir, -1, &ls, nullptr, &rlist_new);

        if (rlist_new != ir->rlist)
        {
            if (fplog != nullptr)
            {
                fprintf(fplog, "\nChanging rlist from %g to %g for non-bonded %dx%d atom kernels\n\n",
                        ir->rlist, rlist_new,
                        ls.cluster_size_i, ls.cluster_size_j);
            }
            ir->rlist     = rlist_new;
        }
    }

    if (nstlist_cmdline > 0 && (!EI_DYNAMICS(ir->eI) || ir->verletbuf_tol <= 0))
    {
        gmx_fatal(FARGS, "Can not set nstlist without %s",
                  !EI_DYNAMICS(ir->eI) ? "dynamics" : "verlet-buffer-tolerance");
    }

    if (EI_DYNAMICS(ir->eI))
    {
        /* Set or try nstlist values */
        increase_nstlist(fplog, cr, ir, nstlist_cmdline, mtop, box, bUseGPU, cpuinfo);
    }
}

/*! \brief Override the nslist value in inputrec
 *
 * with value passed on the command line (if any)
 */
static void override_nsteps_cmdline(const gmx::MDLogger &mdlog,
                                    gmx_int64_t          nsteps_cmdline,
                                    t_inputrec          *ir)
{
    assert(ir);

    /* override with anything else than the default -2 */
    if (nsteps_cmdline > -2)
    {
        char sbuf_steps[STEPSTRSIZE];
        char sbuf_msg[STRLEN];

        ir->nsteps = nsteps_cmdline;
        if (EI_DYNAMICS(ir->eI) && nsteps_cmdline != -1)
        {
            sprintf(sbuf_msg, "Overriding nsteps with value passed on the command line: %s steps, %.3g ps",
                    gmx_step_str(nsteps_cmdline, sbuf_steps),
                    fabs(nsteps_cmdline*ir->delta_t));
        }
        else
        {
            sprintf(sbuf_msg, "Overriding nsteps with value passed on the command line: %s steps",
                    gmx_step_str(nsteps_cmdline, sbuf_steps));
        }

        GMX_LOG(mdlog.warning).asParagraph().appendText(sbuf_msg);
    }
    else if (nsteps_cmdline < -2)
    {
        gmx_fatal(FARGS, "Invalid nsteps value passed on the command line: %d",
                  nsteps_cmdline);
    }
    /* Do nothing if nsteps_cmdline == -2 */
}

//! \brief Return the correct integrator function.
static integrator_t *my_integrator(unsigned int ei)
{
    switch (ei)
    {
        case eiMD:
        case eiBD:
        case eiSD1:
        case eiVV:
        case eiVVAK:
            if (!EI_DYNAMICS(ei))
            {
                GMX_THROW(APIError("do_md integrator would be called for a non-dynamical integrator"));
            }
            return do_md;
        case eiSteep:
            return do_steep;
        case eiCG:
            return do_cg;
        case eiNM:
            return do_nm;
        case eiLBFGS:
            return do_lbfgs;
        case eiTPI:
        case eiTPIC:
            if (!EI_TPI(ei))
            {
                GMX_THROW(APIError("do_tpi integrator would be called for a non-TPI integrator"));
            }
            return do_tpi;
        case eiSD2_REMOVED:
            GMX_THROW(NotImplementedError("SD2 integrator has been removed"));
        default:
            GMX_THROW(APIError("Non existing integrator selected"));
    }
}

//! Initializes the logger for mdrun.
static gmx::LoggerOwner buildLogger(FILE *fplog, const t_commrec *cr)
{
    gmx::LoggerBuilder builder;
    if (fplog != nullptr)
    {
        builder.addTargetFile(gmx::MDLogger::LogLevel::Info, fplog);
    }
    if (cr == nullptr || SIMMASTER(cr))
    {
        builder.addTargetStream(gmx::MDLogger::LogLevel::Warning,
                                &gmx::TextOutputFile::standardError());
    }
    return builder.build();
}


/// \endcond

/// \cond internal
int mdrunner(mdrunner_arglist args)
{
    return mdrunner(&args.hw_opt, args.fplog, args.cr, args.nfile, args.fnm, args.oenv,
                    args.bVerbose, args.nstglobalcomm,
                    args.ddxyz, args.dd_rank_order, args.npme, args.rdd,
                    args.rconstr, args.dddlb_opt, args.dlb_scale,
                    args.ddcsx, args.ddcsy, args.ddcsz,
                    args.nbpu_opt, args.nstlist_cmdline,
                    args.nsteps_cmdline, args.nstepout, args.resetstep,
                    args.nmultisim, *args.replExParams, args.pforce,
                    args.cpt_period, args.max_hours, args.imdport, args.Flags);
}
/// \endcond

// Already documented in runner.h.
/*! \cond internal
 * \brief Make parameter structure for driver routine.
 */
mdrunner_arglist make_mdrunner_arglist(gmx_hw_opt_t *hw_opt,
                                       FILE *fplog, t_commrec *cr, int nfile,
                                       const t_filenm fnm[], const gmx_output_env_t *oenv, gmx_bool bVerbose,
                                       int nstglobalcomm, ivec ddxyz, int dd_rank_order, int npme,
                                       real rdd, real rconstr, const char *dddlb_opt, real dlb_scale,
                                       const char *ddcsx, const char *ddcsy, const char *ddcsz,
                                       const char *nbpu_opt, int nstlist_cmdline,
                                       gmx_int64_t nsteps_cmdline, int nstepout, int resetstep,
                                       int nmultisim,
                                       const ReplicaExchangeParameters &replExParams,
                                       real pforce, real cpt_period, real max_hours,
                                       int imdport, unsigned long Flags)
{
    mdrunner_arglist retval;   // Create object to be returned.
    retval.hw_opt          = *hw_opt;
    retval.fplog           = fplog;
    retval.cr              = cr;
    retval.nfile           = nfile;
    retval.fnm             = fnm;
    retval.oenv            = oenv;
    retval.bVerbose        = bVerbose;
    retval.nstglobalcomm   = nstglobalcomm;
    retval.ddxyz[XX]       = ddxyz[XX];
    retval.ddxyz[YY]       = ddxyz[YY];
    retval.ddxyz[ZZ]       = ddxyz[ZZ];
    retval.dd_rank_order   = dd_rank_order;
    retval.npme            = npme;
    retval.rdd             = rdd;
    retval.rconstr         = rconstr;
    retval.dddlb_opt       = dddlb_opt;
    retval.dlb_scale       = dlb_scale;
    retval.ddcsx           = ddcsx;
    retval.ddcsy           = ddcsy;
    retval.ddcsz           = ddcsz;
    retval.nbpu_opt        = nbpu_opt;
    retval.nstlist_cmdline = nstlist_cmdline;
    retval.nsteps_cmdline  = nsteps_cmdline;
    retval.nstepout        = nstepout;
    retval.resetstep       = resetstep;
    retval.nmultisim       = nmultisim;
    retval.replExParams    = &replExParams;
    retval.pforce          = pforce;
    retval.cpt_period      = cpt_period;
    retval.max_hours       = max_hours;
    retval.imdport         = imdport;
    retval.Flags           = Flags;
    return retval;
}
/// \endcond


/*! \cond internal
 *
 * \brief Implementation of the runner. Overloads declared function in header file.
 *
 * \param[in] hw_opt   Hardware detection structure
 * \param[in] fplog    File pointer for log file
 * \param[in] cr       Communication data
 * \param[in] nfile    Number of files
 * \param[in] fnm      Array of filenames and file properties
 * \param[in] oenv     Output variables for storing xvg files etc.
 * \param[in] bVerbose Verbose output or not
 * \param[in] nstglobalcomm Number of steps between global communication
 * \param[in] ddxyz    Division of sub-boxes over processors for
 *                     use in domain decomposition parallellization
 * \param[in] dd_rank_order Ordering of the PP and PME ranks
 * \param[in] npme     The number of separate PME ranks requested, -1 = auto
 * \param[in] rdd      The maximum distance for bonded interactions with DD (nm)
 * \param[in] rconstr  Maximum distance for P-LINCS (nm)
 * \param[in] dddlb_opt File name for debugging
 * \param[in] dlb_scale File name for debugging
 * \param[in] ddcsx     File name for debugging
 * \param[in] ddcsy     File name for debugging
 * \param[in] ddcsz     File name for debugging
 * \param[in] nbpu_opt  Type of nonbonded processing unit
 * \param[in] nstlist_cmdline  Override neighbor search frequency
 * \param[in] nsteps_cmdline   Override number of simulation steps
 * \param[in] nstepout     How often to write to the console
 * \param[in] resetstep    Reset the step counter
 * \param[in] nmultisim    Number of parallel simulations to run
 * \param[in] replExParams Parameters for the replica exchange algorithm
 * \param[in] pforce       Minimum force for printing (for debugging)
 * \param[in] cpt_period    How often to checkpoint the simulation
 * \param[in] max_hours     Maximume length of the simulation (wall time)
 * \param[in] imdport       Interactive MD port (socket)
 * \param[in] Flags         More command line options
 */
int mdrunner(gmx_hw_opt_t *hw_opt,
             FILE *fplog, t_commrec *cr, int nfile,
             const t_filenm fnm[], const gmx_output_env_t *oenv, gmx_bool bVerbose,
             int nstglobalcomm,
             ivec ddxyz, int dd_rank_order, int npme, real rdd, real rconstr,
             const char *dddlb_opt, real dlb_scale,
             const char *ddcsx, const char *ddcsy, const char *ddcsz,
             const char *nbpu_opt, int nstlist_cmdline,
             gmx_int64_t nsteps_cmdline, int nstepout, int resetstep,
             int gmx_unused nmultisim,
             const ReplicaExchangeParameters &replExParams,
             real pforce, real cpt_period, real max_hours,
             int imdport, unsigned long Flags)
{
    gmx_bool                  bForceUseGPU, bTryUseGPU, bRerunMD;
    matrix                    box;
    gmx_ddbox_t               ddbox = {0};
    int                       npme_major, npme_minor;
    t_nrnb                   *nrnb;
    gmx_mtop_t               *mtop          = nullptr;
    t_mdatoms                *mdatoms       = nullptr;
    t_forcerec               *fr            = nullptr;
    t_fcdata                 *fcd           = nullptr;
    real                      ewaldcoeff_q  = 0;
    real                      ewaldcoeff_lj = 0;
    struct gmx_pme_t        **pmedata       = nullptr;
    gmx_vsite_t              *vsite         = nullptr;
    gmx_constr_t              constr;
    int                       nChargePerturbed = -1, nTypePerturbed = 0, status;
    gmx_wallcycle_t           wcycle;
    gmx_walltime_accounting_t walltime_accounting = nullptr;
    int                       rc;
    gmx_int64_t               reset_counters;
    int                       nthreads_pme = 1;
    gmx_membed_t *            membed       = nullptr;
    gmx_hw_info_t            *hwinfo       = nullptr;
    /* The master rank decides early on bUseGPU and broadcasts this later */
    gmx_bool                  bUseGPU            = FALSE;

    /* CAUTION: threads may be started later on in this function, so
       cr doesn't reflect the final parallel state right now */
    gmx::MDModules mdModules;
    t_inputrec     inputrecInstance;
    t_inputrec    *inputrec = &inputrecInstance;
    snew(mtop, 1);

    if (Flags & MD_APPENDFILES)
    {
        fplog = nullptr;
    }

    bool doMembed = opt2bSet("-membed", nfile, fnm);
    bRerunMD     = (Flags & MD_RERUN);

    /* Handle GPU-related user options. Later, we check consistency
     * with things like whether support is compiled, or tMPI thread
     * count. */
    bForceUseGPU = (strncmp(nbpu_opt, "gpu", 3) == 0);
    bTryUseGPU   = (strncmp(nbpu_opt, "auto", 4) == 0) || bForceUseGPU;
    gmx_parse_gpu_ids(&hw_opt->gpu_opt);

    // Here we assume that SIMMASTER(cr) does not change even after the
    // threads are started.
    gmx::LoggerOwner logOwner(buildLogger(fplog, cr));
    gmx::MDLogger    mdlog(logOwner.logger());

    /* Detect hardware, gather information. This is an operation that is
     * global for this process (MPI rank). */
    hwinfo = gmx_detect_hardware(mdlog, cr, bTryUseGPU);

    gmx_print_detected_hardware(fplog, cr, mdlog, hwinfo);

    if (fplog != nullptr)
    {
        /* Print references after all software/hardware printing */
        please_cite(fplog, "Abraham2015");
        please_cite(fplog, "Pall2015");
        please_cite(fplog, "Pronk2013");
        please_cite(fplog, "Hess2008b");
        please_cite(fplog, "Spoel2005a");
        please_cite(fplog, "Lindahl2001a");
        please_cite(fplog, "Berendsen95a");
    }

    std::unique_ptr<t_state> stateInstance = std::unique_ptr<t_state>(new t_state);
    t_state *                state         = stateInstance.get();

    if (SIMMASTER(cr))
    {
        /* Read (nearly) all data required for the simulation */
        read_tpx_state(ftp2fn(efTPR, nfile, fnm), inputrec, state, mtop);

        if (inputrec->cutoff_scheme == ecutsVERLET)
        {
            /* Here the master rank decides if all ranks will use GPUs */
            bUseGPU = (hwinfo->gpu_info.n_dev_compatible > 0 ||
                       getenv("GMX_EMULATE_GPU") != nullptr);

            /* TODO add GPU kernels for this and replace this check by:
             * (bUseGPU && (ir->vdwtype == evdwPME &&
             *               ir->ljpme_combination_rule == eljpmeLB))
             * update the message text and the content of nbnxn_acceleration_supported.
             */
            if (bUseGPU &&
                !nbnxn_gpu_acceleration_supported(mdlog, inputrec, bRerunMD))
            {
                /* Fallback message printed by nbnxn_acceleration_supported */
                if (bForceUseGPU)
                {
                    gmx_fatal(FARGS, "GPU acceleration requested, but not supported with the given input settings");
                }
                bUseGPU = FALSE;
            }

            prepare_verlet_scheme(fplog, cr,
                                  inputrec, nstlist_cmdline, mtop, state->box,
                                  bUseGPU, *hwinfo->cpuInfo);
        }
        else
        {
            if (nstlist_cmdline > 0)
            {
                gmx_fatal(FARGS, "Can not set nstlist with the group cut-off scheme");
            }

            if (hwinfo->gpu_info.n_dev_compatible > 0)
            {
                GMX_LOG(mdlog.warning).asParagraph().appendText(
                        "NOTE: GPU(s) found, but the current simulation can not use GPUs\n"
                        "      To use a GPU, set the mdp option: cutoff-scheme = Verlet");
            }

            if (bForceUseGPU)
            {
                gmx_fatal(FARGS, "GPU requested, but can't be used without cutoff-scheme=Verlet");
            }

#if GMX_TARGET_BGQ
            md_print_warn(commRec_, fplog,
                          "NOTE: There is no SIMD implementation of the group scheme kernels on\n"
                          "      BlueGene/Q. You will observe better performance from using the\n"
                          "      Verlet cut-off scheme.\n");
#endif
        }
    }

    /* Check and update the hardware options for internal consistency */
    check_and_update_hw_opt_1(hw_opt, cr, npme);

    /* Early check for externally set process affinity. */
    gmx_check_thread_affinity_set(mdlog, cr,
                                  hw_opt, hwinfo->nthreads_hw_avail, FALSE);

#if GMX_THREAD_MPI
    if (SIMMASTER(cr))
    {
        if (npme > 0 && hw_opt->nthreads_tmpi <= 0)
        {
            gmx_fatal(FARGS, "You need to explicitly specify the number of MPI threads (-ntmpi) when using separate PME ranks");
        }

        /* Since the master knows the cut-off scheme, update hw_opt for this.
         * This is done later for normal MPI and also once more with tMPI
         * for all tMPI ranks.
         */
        check_and_update_hw_opt_2(hw_opt, inputrec->cutoff_scheme);

        /* NOW the threads will be started: */
        hw_opt->nthreads_tmpi = get_nthreads_mpi(hwinfo,
                                                 hw_opt,
                                                 inputrec, mtop,
                                                 mdlog, bUseGPU,
                                                 doMembed);

        if (hw_opt->nthreads_tmpi > 1)
        {
            t_commrec *cr_old       = cr;
            /* now start the threads. */
            cr = mdrunner_start_threads(hw_opt, fplog, cr_old, nfile, fnm,
                                        oenv, bVerbose, nstglobalcomm,
                                        ddxyz, dd_rank_order, npme, rdd, rconstr,
                                        dddlb_opt, dlb_scale, ddcsx, ddcsy, ddcsz,
                                        nbpu_opt, nstlist_cmdline,
                                        nsteps_cmdline, nstepout, resetstep, nmultisim,
                                        replExParams, pforce,
                                        cpt_period, max_hours,
                                        Flags);
            /* the main thread continues here with a new cr. We don't deallocate
               the old cr because other threads may still be reading it. */
            if (cr == nullptr)
            {
                gmx_comm("Failed to spawn threads");
            }
        }
    }
#endif
    /* END OF CAUTION: cr is now reliable */

    if (PAR(cr))
    {
        /* now broadcast everything to the non-master nodes/threads: */
        init_parallel(cr, inputrec, mtop);

        /* The master rank decided on the use of GPUs,
         * broadcast this information to all ranks.
         */
        gmx_bcast_sim(sizeof(bUseGPU), &bUseGPU, cr);
    }
    // TODO: Error handling
    mdModules.assignOptionsToModules(*inputrec->params, nullptr);

    if (fplog != nullptr)
    {
        pr_inputrec(fplog, 0, "Input Parameters", inputrec, FALSE);
        fprintf(fplog, "\n");
    }

    /* now make sure the state is initialized and propagated */
    set_state_entries(state, inputrec);

    /* A parallel command line option consistency check that we can
       only do after any threads have started. */
    if (!PAR(cr) &&
        (ddxyz[XX] > 1 || ddxyz[YY] > 1 || ddxyz[ZZ] > 1 || npme > 0))
    {
        gmx_fatal(FARGS,
                  "The -dd or -npme option request a parallel simulation, "
#if !GMX_MPI
                  "but %s was compiled without threads or MPI enabled"
#else
#if GMX_THREAD_MPI
                  "but the number of MPI-threads (option -ntmpi) is not set or is 1"
#else
                  "but %s was not started through mpirun/mpiexec or only one rank was requested through mpirun/mpiexec"
#endif
#endif
                  , output_env_get_program_display_name(oenv)
                  );
    }

    if (bRerunMD &&
        (EI_ENERGY_MINIMIZATION(inputrec->eI) || eiNM == inputrec->eI))
    {
        gmx_fatal(FARGS, "The .mdp file specified an energy mininization or normal mode algorithm, and these are not compatible with mdrun -rerun");
    }

    if (can_use_allvsall(inputrec, TRUE, cr, fplog) && DOMAINDECOMP(cr))
    {
        gmx_fatal(FARGS, "All-vs-all loops do not work with domain decomposition, use a single MPI rank");
    }

    if (!(EEL_PME(inputrec->coulombtype) || EVDW_PME(inputrec->vdwtype)))
    {
        if (npme > 0)
        {
            gmx_fatal_collective(FARGS, cr->mpi_comm_mysim, MASTER(cr),
                                 "PME-only ranks are requested, but the system does not use PME for electrostatics or LJ");
        }

        npme = 0;
    }

    if (bUseGPU && npme < 0)
    {
        /* With GPUs we don't automatically use PME-only ranks. PME ranks can
         * improve performance with many threads per GPU, since our OpenMP
         * scaling is bad, but it's difficult to automate the setup.
         */
        npme = 0;
    }

#ifdef GMX_FAHCORE
    if (MASTER(commRec_))
    {
        fcRegisterSteps(inputrec->nsteps, inputrec->init_step);
    }
#endif

    /* NMR restraints must be initialized before load_checkpoint,
     * since with time averaging the history is added to t_state.
     * For proper consistency check we therefore need to extend
     * t_state here.
     * So the PME-only nodes (if present) will also initialize
     * the distance restraints.
     */
    snew(fcd, 1);

    /* This needs to be called before read_checkpoint to extend the state */
    init_disres(fplog, mtop, inputrec, cr, fcd, state, replExParams.exchangeInterval > 0);

    init_orires(fplog, mtop, as_rvec_array(state->x.data()), inputrec, cr, &(fcd->orires),
                state);

    if (inputrecDeform(inputrec))
    {
        /* Store the deform reference box before reading the checkpoint */
        if (SIMMASTER(cr))
        {
            copy_mat(state->box, box);
        }
        if (PAR(cr))
        {
            gmx_bcast(sizeof(box), box, cr);
        }
        /* Because we do not have the update struct available yet
         * in which the reference values should be stored,
         * we store them temporarily in static variables.
         * This should be thread safe, since they are only written once
         * and with identical values.
         */
        tMPI_Thread_mutex_lock(&deform_init_box_mutex);
        deform_init_init_step_tpx = inputrec->init_step;
        copy_mat(box, deform_init_box_tpx);
        tMPI_Thread_mutex_unlock(&deform_init_box_mutex);
    }

    ObservablesHistory observablesHistory = {};

    if (Flags & MD_STARTFROMCPT)
    {
        /* Check if checkpoint file exists before doing continuation.
         * This way we can use identical input options for the first and subsequent runs...
         */
        gmx_bool bReadEkin;

        load_checkpoint(opt2fn_master("-cpi", nfile, fnm, cr), &fplog,
                        cr, ddxyz, &npme,
                        inputrec, state, &bReadEkin, &observablesHistory,
                        (Flags & MD_APPENDFILES),
                        (Flags & MD_APPENDFILESSET),
                        (Flags & MD_REPRODUCIBLE));

        if (bReadEkin)
        {
            Flags |= MD_READ_EKIN;
        }
    }

    if (SIMMASTER(cr) && (Flags & MD_APPENDFILES))
    {
        gmx_log_open(ftp2fn(efLOG, nfile, fnm), cr,
                     Flags, &fplog);
        logOwner = buildLogger(fplog, nullptr);
        mdlog    = logOwner.logger();
    }

    /* override nsteps with value from cmdline */
    override_nsteps_cmdline(mdlog, nsteps_cmdline, inputrec);

    if (SIMMASTER(cr))
    {
        copy_mat(state->box, box);
    }

    if (PAR(cr))
    {
        gmx_bcast(sizeof(box), box, cr);
    }

    if (PAR(cr) && !(EI_TPI(inputrec->eI) ||
                     inputrec->eI == eiNM))
    {
        gmx::DomDecParams ddParams {
            Flags,
            npme,
            dd_rank_order,
            rdd,
            rconstr,
            dddlb_opt,
            dlb_scale,
            ddcsx,
            ddcsy,
            ddcsz
        };
        cr->dd = init_domain_decomposition(fplog,
                                           cr,
                                           ddParams,
                                           ddxyz,
                                           mtop, inputrec,
                                           box, as_rvec_array(state->x.data()),
                                           &ddbox, &npme_major, &npme_minor);
    }
    else
    {
        /* PME, if used, is done on all nodes with 1D decomposition */
        cr->npmenodes = 0;
        cr->duty      = (DUTY_PP | DUTY_PME);
        npme_major    = 1;
        npme_minor    = 1;

        if (inputrec->ePBC == epbcSCREW)
        {
            gmx_fatal(FARGS,
                      "pbc=%s is only implemented with domain decomposition",
                      epbc_names[inputrec->ePBC]);
        }
    }

    if (PAR(cr))
    {
        /* After possible communicator splitting in make_dd_communicators.
         * we can set up the intra/inter node communication.
         */
        gmx_setup_nodecomm(fplog, cr);
    }

    /* Initialize per-physical-node MPI process/thread ID and counters. */
    gmx_init_intranode_counters(cr);
#if GMX_MPI
    if (MULTISIM(cr))
    {
        GMX_LOG(mdlog.warning).asParagraph().appendTextFormatted(
                "This is simulation %d out of %d running as a composite GROMACS\n"
                "multi-simulation job. Setup for this simulation:\n",
                cr->ms->sim, cr->ms->nsim);
    }
    GMX_LOG(mdlog.warning).appendTextFormatted(
            "Using %d MPI %s\n",
            cr->nnodes,
#if GMX_THREAD_MPI
            cr->nnodes == 1 ? "thread" : "threads"
#else
            cr->nnodes == 1 ? "process" : "processes"
#endif
            );
    fflush(stderr);
#endif

    /* Check and update hw_opt for the cut-off scheme */
    check_and_update_hw_opt_2(hw_opt, inputrec->cutoff_scheme);

    /* Check and update hw_opt for the number of MPI ranks */
    check_and_update_hw_opt_3(hw_opt);

    gmx_omp_nthreads_init(mdlog, cr,
                          hwinfo->nthreads_hw_avail,
                          hw_opt->nthreads_omp,
                          hw_opt->nthreads_omp_pme,
                          (cr->duty & DUTY_PP) == 0,
                          inputrec->cutoff_scheme == ecutsVERLET);

#ifndef NDEBUG
    if (EI_TPI(inputrec->eI) &&
        inputrec->cutoff_scheme == ecutsVERLET)
    {
        gmx_feenableexcept();
    }
#endif

    bool userSetGpuIds = hasUserSetGpuIds(&hw_opt->gpu_opt);

    if (bUseGPU)
    {
        /* Select GPU id's to use */
        gmx_select_rank_gpu_ids(mdlog, cr, &hwinfo->gpu_info, bForceUseGPU,
                                userSetGpuIds, &hw_opt->gpu_opt);
    }
    else
    {
        /* Ignore (potentially) manually selected GPUs */
        hw_opt->gpu_opt.n_dev_use = 0;
    }

    /* check consistency across ranks of things like SIMD
     * support and number of GPUs selected */
    gmx_check_hw_runconf_consistency(mdlog, hwinfo, cr, hw_opt, userSetGpuIds, bUseGPU);

    /* Now that we know the setup is consistent, check for efficiency */
    check_resource_division_efficiency(hwinfo, hw_opt, hw_opt->gpu_opt.n_dev_use, Flags & MD_NTOMPSET,
                                       cr, mdlog);

    if (DOMAINDECOMP(cr))
    {
        /* When we share GPUs over ranks, we need to know this for the DLB */
        dd_setup_dlb_resource_sharing(cr, hwinfo, hw_opt);
    }

    /* getting number of PP/PME threads
       PME: env variable should be read only on one node to make sure it is
       identical everywhere;
     */
    nthreads_pme = gmx_omp_nthreads_get(emntPME);

    wcycle = wallcycle_init(fplog, resetstep, cr);

    if (PAR(cr))
    {
        /* Master synchronizes its value of reset_counters with all nodes
         * including PME only nodes */
        reset_counters = wcycle_get_reset_counters(wcycle);
        gmx_bcast_sim(sizeof(reset_counters), &reset_counters, cr);
        wcycle_set_reset_counters(wcycle, reset_counters);
    }

    // Membrane embedding must be initialized before we call init_forcerec()
    if (doMembed)
    {
        if (MASTER(cr))
        {
            fprintf(stderr, "Initializing membed");
        }
        /* Note that membed cannot work in parallel because mtop is
         * changed here. Fix this if we ever want to make it run with
         * multiple ranks. */
        membed = init_membed(fplog, nfile, fnm, mtop, inputrec, state, cr, &cpt_period);
    }

    snew(nrnb, 1);
    if (cr->duty & DUTY_PP)
    {
        bcast_state(cr, state);

        /* Initiate forcerecord */
        fr          = mk_forcerec();
        fr->hwinfo  = hwinfo;
        fr->gpu_opt = &hw_opt->gpu_opt;
        init_forcerec(fplog, mdlog, fr, fcd, mdModules.forceProvider(),
                      inputrec, mtop, cr, box,
                      opt2fn("-table", nfile, fnm),
                      opt2fn("-tablep", nfile, fnm),
                      getFilenm("-tableb", nfile, fnm),
                      nbpu_opt,
                      FALSE,
                      pforce);

        /* Initialize QM-MM */
        if (fr->bQMMM)
        {
            init_QMMMrec(cr, mtop, inputrec, fr);
        }

        /* Initialize the mdatoms structure.
         * mdatoms is not filled with atom data,
         * as this can not be done now with domain decomposition.
         */
        mdatoms = init_mdatoms(fplog, mtop, inputrec->efep != efepNO);

        /* Initialize the virtual site communication */
        vsite = init_vsite(mtop, cr, FALSE);

        calc_shifts(box, fr->shift_vec);

        /* With periodic molecules the charge groups should be whole at start up
         * and the virtual sites should not be far from their proper positions.
         */
        if (!inputrec->bContinuation && MASTER(cr) &&
            !(inputrec->ePBC != epbcNONE && inputrec->bPeriodicMols))
        {
            /* Make molecules whole at start of run */
            if (fr->ePBC != epbcNONE)
            {
                do_pbc_first_mtop(fplog, inputrec->ePBC, box, mtop, as_rvec_array(state->x.data()));
            }
            if (vsite)
            {
                /* Correct initial vsite positions are required
                 * for the initial distribution in the domain decomposition
                 * and for the initial shell prediction.
                 */
                construct_vsites_mtop(vsite, mtop, as_rvec_array(state->x.data()));
            }
        }

        if (EEL_PME(fr->eeltype) || EVDW_PME(fr->vdwtype))
        {
            ewaldcoeff_q  = fr->ewaldcoeff_q;
            ewaldcoeff_lj = fr->ewaldcoeff_lj;
            pmedata       = &fr->pmedata;
        }
        else
        {
            pmedata = nullptr;
        }
    }
    else
    {
        /* This is a PME only node */

        /* We don't need the state */
        stateInstance.reset();
        state         = nullptr;

        ewaldcoeff_q  = calc_ewaldcoeff_q(inputrec->rcoulomb, inputrec->ewald_rtol);
        ewaldcoeff_lj = calc_ewaldcoeff_lj(inputrec->rvdw, inputrec->ewald_rtol_lj);
        snew(pmedata, 1);
    }

    if (hw_opt->thread_affinity != threadaffOFF)
    {
        /* Before setting affinity, check whether the affinity has changed
         * - which indicates that probably the OpenMP library has changed it
         * since we first checked).
         */
        gmx_check_thread_affinity_set(mdlog, cr,
                                      hw_opt, hwinfo->nthreads_hw_avail, TRUE);

        int nthread_local;
        /* threads on this MPI process or TMPI thread */
        if (cr->duty & DUTY_PP)
        {
            nthread_local = gmx_omp_nthreads_get(emntNonbonded);
        }
        else
        {
            nthread_local = gmx_omp_nthreads_get(emntPME);
        }

        /* Set the CPU affinity */
        gmx_set_thread_affinity(mdlog, cr, hw_opt, *hwinfo->hardwareTopology,
                                nthread_local, nullptr);
    }

    /* Initiate PME if necessary,
     * either on all nodes or on dedicated PME nodes only. */
    if (EEL_PME(inputrec->coulombtype) || EVDW_PME(inputrec->vdwtype))
    {
        if (mdatoms)
        {
            nChargePerturbed = mdatoms->nChargePerturbed;
            if (EVDW_PME(inputrec->vdwtype))
            {
                nTypePerturbed   = mdatoms->nTypePerturbed;
            }
        }
        if (cr->npmenodes > 0)
        {
            /* The PME only nodes need to know nChargePerturbed(FEP on Q) and nTypePerturbed(FEP on LJ)*/
            gmx_bcast_sim(sizeof(nChargePerturbed), &nChargePerturbed, cr);
            gmx_bcast_sim(sizeof(nTypePerturbed), &nTypePerturbed, cr);
        }

        if (cr->duty & DUTY_PME)
        {
            try
            {
                status = gmx_pme_init(pmedata, cr, npme_major, npme_minor, inputrec,
                                      mtop ? mtop->natoms : 0, nChargePerturbed, nTypePerturbed,
                                      (Flags & MD_REPRODUCIBLE),
                                      ewaldcoeff_q, ewaldcoeff_lj,
                                      nthreads_pme);
            }
            GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
            if (status != 0)
            {
                gmx_fatal(FARGS, "Error %d initializing PME", status);
            }
        }
    }


    if (EI_DYNAMICS(inputrec->eI))
    {
        /* Turn on signal handling on all nodes */
        /*
         * (A user signal from the PME nodes (if any)
         * is communicated to the PP nodes.
         */
        signal_handler_install();
    }

    if (cr->duty & DUTY_PP)
    {
        /* Assumes uniform use of the number of OpenMP threads */
        walltime_accounting = walltime_accounting_init(gmx_omp_nthreads_get(emntDefault));

        if (inputrec->bPull)
        {
            /* Initialize pull code */
            inputrec->pull_work =
                init_pull(fplog, inputrec->pull, inputrec, nfile, fnm,
                          mtop, cr, oenv, inputrec->fepvals->init_lambda,
                          EI_DYNAMICS(inputrec->eI) && MASTER(cr), Flags);
        }

        if (inputrec->bRot)
        {
            /* Initialize enforced rotation code */
            init_rot(fplog, inputrec, nfile, fnm, cr, as_rvec_array(state->x.data()), state->box, mtop, oenv,
                     bVerbose, Flags);
        }

        /* Let init_constraints know whether we have essential dynamics constraints.
         * TODO: inputrec should tell us whether we use an algorithm, not a file option or the checkpoint
         */
        bool doEdsam = (opt2fn_null("-ei", nfile, fnm) != nullptr || observablesHistory.edsamHistory);

        constr = init_constraints(fplog, mtop, inputrec, doEdsam, cr);

        if (DOMAINDECOMP(cr))
        {
            GMX_RELEASE_ASSERT(fr, "fr was NULL while cr->duty was DUTY_PP");
            /* This call is not included in init_domain_decomposition mainly
             * because fr->cginfo_mb is set later.
             */
            dd_init_bondeds(fplog, cr->dd, mtop, vsite, inputrec,
                            Flags & MD_DDBONDCHECK, fr->cginfo_mb);
        }

        /* Now do whatever the user wants us to do (how flexible...) */
        my_integrator(inputrec->eI) (fplog, cr, mdlog, nfile, fnm,
                                     oenv, bVerbose,
                                     nstglobalcomm,
                                     vsite, constr,
                                     nstepout, mdModules.outputProvider(),
                                     inputrec, mtop,
                                     fcd, state, &observablesHistory,
                                     mdatoms, nrnb, wcycle, fr,
                                     replExParams,
                                     membed,
                                     cpt_period, max_hours,
                                     imdport,
                                     Flags,
                                     walltime_accounting);

        if (inputrec->bRot)
        {
            finish_rot(inputrec->rot);
        }

        if (inputrec->bPull)
        {
            finish_pull(inputrec->pull_work);
        }

    }
    else
    {
        GMX_RELEASE_ASSERT(pmedata, "pmedata was NULL while cr->duty was not DUTY_PP");
        /* do PME only */
        walltime_accounting = walltime_accounting_init(gmx_omp_nthreads_get(emntPME));
        gmx_pmeonly(*pmedata, cr, nrnb, wcycle, walltime_accounting, ewaldcoeff_q, ewaldcoeff_lj, inputrec);
    }

    wallcycle_stop(wcycle, ewcRUN);

    /* Finish up, write some stuff
     * if rerunMD, don't write last frame again
     */
    finish_run(fplog, mdlog, cr,
               inputrec, nrnb, wcycle, walltime_accounting,
               fr ? fr->nbv : nullptr,
               EI_DYNAMICS(inputrec->eI) && !MULTISIM(cr));

    // Free PME data
    if (pmedata)
    {
        gmx_pme_destroy(*pmedata); // TODO: pmedata is always a single element list, refactor
        pmedata = nullptr;
    }

    /* Free GPU memory and context */
    free_gpu_resources(fr, cr, &hwinfo->gpu_info, fr ? fr->gpu_opt : nullptr);

    if (doMembed)
    {
        free_membed(membed);
    }

    gmx_hardware_info_free(hwinfo);

    /* Does what it says */
    print_date_and_time(fplog, cr->nodeid, "Finished mdrun", gmx_gettime());
    walltime_accounting_destroy(walltime_accounting);

    /* Close logfile already here if we were appending to it */
    if (MASTER(cr) && (Flags & MD_APPENDFILES))
    {
        gmx_log_close(fplog);
    }

    rc = (int)gmx_get_stop_condition();

#if GMX_THREAD_MPI
    /* we need to join all threads. The sub-threads join when they
       exit this function, but the master thread needs to be told to
       wait for that. */
    if (PAR(cr) && MASTER(cr))
    {
        tMPI_Finalize();
    }
#endif

    return rc;
};
/// \endcond

RunnerImpl::RunnerImpl() :
    flags_ {0UL},
input_ {
    nullptr
},
fplog_ {
    nullptr
},
replExParams_ {
    gmx::compat::make_unique<ReplicaExchangeParameters>()
},
logOwner_ {
    nullptr
},
hardwareInfo_ {
    nullptr
},
pmeData_ {
    nullptr
},
commRec_ {
    nullptr
},
doAppendFiles_ {
    0
},
doMembed_ {
    false
},
ewaldcoeff_q_ {
    real(0)
},
ewaldcoeff_lj_ {
    real(0)
},
mdModules_ {},
oenv_ {
    nullptr
},                  // allocated by parse_common_args
verbose_ {
    false
},
nstglobalcomm_ {
    0
},
ddParams_ {
    gmx::compat::make_unique<gmx::DomDecParams>()
},
observablesHistory_ {
    gmx::compat::make_unique<ObservablesHistory>()
},
imdport_ {
    8888
},                  /* can be almost anything, 8888 is easy to remember */
nstepout_ {
    100
},
cpt_period_ {
    real(15.0)
},
mdAtoms_ {
    nullptr
},
vSite_ {
    nullptr
},
forceCalcData_ {
    gmx::compat::make_unique<t_fcdata>()
},
nrNonBonded_ {
    gmx::compat::make_unique<t_nrnb>()
},
wallCycle_ {
    nullptr
},
forceRecord_ {
    nullptr
},
membed_ {
    nullptr
},
maxHours_ {
    real(-1)
},
    // Allowing the destructor to handle a non-initialized walltimeAccounting is too much work...
walltimeAccounting_ {
    walltime_accounting_init(gmx_omp_nthreads_get(emntDefault))
},
initialized_{false}
{}

RunnerImpl::RunnerImpl(const std::string &filename) :
    RunnerImpl {}
{
    initialized_ = true;
    input_ = gmxapi::MDInput::from_tpr_file(filename);
    auto &input = *input_;

    int   argc = 3;
    char* argv[3];
    // Put these on the stack so I don't have to worry about it.
    char  argv0[] = "";   // the parser skips what it assumes is the command name
    char  argv1[] = "-s"; // the topology / TPR file
    char  filename_stack_string[filename.length()+1];
    filename.copy(filename_stack_string, filename.length());
    filename_stack_string[filename.length()] = '\0';

    argv[0] = argv0;
    argv[1] = argv1;
    argv[2] = filename_stack_string;

    // //! Implements C-style main function for mdrun
    // int gmx_mdrun(int argc, char *argv[])
    // {
    const char   *desc[] = {
        ""  // don't need this right now, but will need to deal with it later...
    };
    auto        * fnm   = fnm_.data();
    const int     nfile = fnm_.size();

    /* Command line option parameters, with their default values */
    gmx_bool          bDDBondCheck  = TRUE;
    gmx_bool          bDDBondComm   = TRUE;
    gmx_bool          bTunePME      = TRUE;
    bool             &bVerbose      = verbose_;
    gmx_bool          bRerunVSite   = FALSE;
    gmx_bool          bConfout      = TRUE;
    gmx_bool          bReproducible = FALSE;
    gmx_bool          bIMDwait      = FALSE;
    gmx_bool          bIMDterm      = FALSE;
    gmx_bool          bIMDpull      = FALSE;

    int               nstlist       = 0;
    int               nmultisim     = 0;
    int               nstglobalcomm = -1;
    int               resetstep     = -1;
    gmx_int64_t       nsteps        = -2;   /* the value -2 means that the mdp option will be used */

    /* Special algorithms section */
    ReplicaExchangeParameters replExParams;

    /* Command line options */
    rvec              realddxyz                   = {0, 0, 0};
    const char       *ddrank_opt[ddrankorderNR+1] =
    { nullptr, "interleave", "pp_pme", "cartesian", nullptr };
    const char       *dddlb_opt[] =
    { nullptr, "auto", "no", "yes", nullptr };
    const char       *thread_aff_opt[threadaffNR+1] =
    { nullptr, "auto", "on", "off", nullptr };
    const char       *nbpu_opt[] =
    { nullptr, "auto", "cpu", "gpu", "gpu_cpu", nullptr };
    real              pforce                = -1;
    gmx_bool          bTryToAppendFiles     = TRUE;
    gmx_bool          bKeepAndNumCPT        = FALSE;
    gmx_bool          bResetCountersHalfWay = FALSE;

    /* Non transparent initialization of a complex gmx_hw_opt_t struct.
     * But unfortunately we are not allowed to call a function here,
     * since declarations follow below.
     */
    gmx_hw_opt_t    hw_opt = {
        0, 0, 0, 0, threadaffSEL, 0, 0,
        { nullptr, 0, nullptr }
    };

    t_pargs         pa[] = {

        { "-dd",      FALSE, etRVEC, {&realddxyz},
          "Domain decomposition grid, 0 is optimize" },
        { "-ddorder", FALSE, etENUM, {ddrank_opt},
          "DD rank order" },
        { "-npme",    FALSE, etINT, {&ddParams_->nPmeRanks},
          "Number of separate ranks to be used for PME, -1 is guess" },
        { "-nt",      FALSE, etINT, {&hw_opt.nthreads_tot},
          "Total number of threads to start (0 is guess)" },
        { "-ntmpi",   FALSE, etINT, {&hw_opt.nthreads_tmpi},
          "Number of thread-MPI threads to start (0 is guess)" },
        { "-ntomp",   FALSE, etINT, {&hw_opt.nthreads_omp},
          "Number of OpenMP threads per MPI rank to start (0 is guess)" },
        { "-ntomp_pme", FALSE, etINT, {&hw_opt.nthreads_omp_pme},
          "Number of OpenMP threads per MPI rank to start (0 is -ntomp)" },
        { "-pin",     FALSE, etENUM, {thread_aff_opt},
          "Whether mdrun should try to set thread affinities" },
        { "-pinoffset", FALSE, etINT, {&hw_opt.core_pinning_offset},
          "The lowest logical core number to which mdrun should pin the first thread" },
        { "-pinstride", FALSE, etINT, {&hw_opt.core_pinning_stride},
          "Pinning distance in logical cores for threads, use 0 to minimize the number of threads per physical core" },
        { "-gpu_id",  FALSE, etSTR, {&hw_opt.gpu_opt.gpu_id},
          "List of GPU device id-s to use, specifies the per-node PP rank to GPU mapping" },
        { "-ddcheck", FALSE, etBOOL, {&bDDBondCheck},
          "Check for all bonded interactions with DD" },
        { "-ddbondcomm", FALSE, etBOOL, {&bDDBondComm},
          "HIDDENUse special bonded atom communication when [TT]-rdd[tt] > cut-off" },
        { "-rdd",     FALSE, etREAL, {&ddParams_->commDistanceMin},
          "The maximum distance for bonded interactions with DD (nm), 0 is determine from initial coordinates" },
        { "-rcon",    FALSE, etREAL, {&ddParams_->rConstraints},
          "Maximum distance for P-LINCS (nm), 0 is estimate" },
        { "-dlb",     FALSE, etENUM, {dddlb_opt},
          "Dynamic load balancing (with DD)" },
        { "-dds",     FALSE, etREAL, {&ddParams_->dlbScale},
          "Fraction in (0,1) by whose reciprocal the initial DD cell size will be increased in order to "
          "provide a margin in which dynamic load balancing can act while preserving the minimum cell size." },
        { "-ddcsx",   FALSE, etSTR, {&ddParams_->sizeX},
          "HIDDENA string containing a vector of the relative sizes in the x "
          "direction of the corresponding DD cells. Only effective with static "
          "load balancing." },
        { "-ddcsy",   FALSE, etSTR, {&ddParams_->sizeY},
          "HIDDENA string containing a vector of the relative sizes in the y "
          "direction of the corresponding DD cells. Only effective with static "
          "load balancing." },
        { "-ddcsz",   FALSE, etSTR, {&ddParams_->sizeZ},
          "HIDDENA string containing a vector of the relative sizes in the z "
          "direction of the corresponding DD cells. Only effective with static "
          "load balancing." },
        { "-gcom",    FALSE, etINT, {&nstglobalcomm},
          "Global communication frequency" },
        { "-nb",      FALSE, etENUM, {&nbpu_opt},
          "Calculate non-bonded interactions on" },
        { "-nstlist", FALSE, etINT, {&nstlist},
          "Set nstlist when using a Verlet buffer tolerance (0 is guess)" },
        { "-tunepme", FALSE, etBOOL, {&bTunePME},
          "Optimize PME load between PP/PME ranks or GPU/CPU" },
        { "-v",       FALSE, etBOOL, {&bVerbose},
          "Be loud and noisy" },
        { "-pforce",  FALSE, etREAL, {&pforce},
          "Print all forces larger than this (kJ/mol nm)" },
        { "-reprod",  FALSE, etBOOL, {&bReproducible},
          "Try to avoid optimizations that affect binary reproducibility" },
        { "-cpt",     FALSE, etREAL, {&cpt_period_},
          "Checkpoint interval (minutes)" },
        { "-cpnum",   FALSE, etBOOL, {&bKeepAndNumCPT},
          "Keep and number checkpoint files" },
        { "-append",  FALSE, etBOOL, {&bTryToAppendFiles},
          "Append to previous output files when continuing from checkpoint instead of adding the simulation part number to all file names" },
        { "-nsteps",  FALSE, etINT64, {&nsteps},
          "Run this number of steps, overrides .mdp file option (-1 means infinite, -2 means use mdp option, smaller is invalid)" },
        { "-maxh",   FALSE, etREAL, {&maxHours_},
          "Terminate after 0.99 times this time (hours)" },
        { "-multi",   FALSE, etINT, {&nmultisim},
          "Do multiple simulations in parallel" },
        { "-replex",  FALSE, etINT, {&replExParams.exchangeInterval},
          "Attempt replica exchange periodically with this period (steps)" },
        { "-nex",  FALSE, etINT, {&replExParams.numExchanges},
          "Number of random exchanges to carry out each exchange interval (N^3 is one suggestion).  -nex zero or not specified gives neighbor replica exchange." },
        { "-reseed",  FALSE, etINT, {&replExParams.randomSeed},
          "Seed for replica exchange, -1 is generate a seed" },
        { "-imdport",    FALSE, etINT, {&imdport_},
          "HIDDENIMD listening port" },
        { "-imdwait",  FALSE, etBOOL, {&bIMDwait},
          "HIDDENPause the simulation while no IMD client is connected" },
        { "-imdterm",  FALSE, etBOOL, {&bIMDterm},
          "HIDDENAllow termination of the simulation from IMD client" },
        { "-imdpull",  FALSE, etBOOL, {&bIMDpull},
          "HIDDENAllow pulling in the simulation from IMD client" },
        { "-rerunvsite", FALSE, etBOOL, {&bRerunVSite},
          "HIDDENRecalculate virtual site coordinates with [TT]-rerun[tt]" },
        { "-confout", FALSE, etBOOL, {&bConfout},
          "HIDDENWrite the last configuration with [TT]-c[tt] and force checkpointing at the last step" },
        { "-stepout", FALSE, etINT, {&nstepout_},
          "HIDDENFrequency of writing the remaining wall clock time for the run" },
        { "-resetstep", FALSE, etINT, {&resetstep},
          "HIDDENReset cycle counters after these many time steps" },
        { "-resethway", FALSE, etBOOL, {&bResetCountersHalfWay},
          "HIDDENReset the cycle counters after half the number of steps or halfway [TT]-maxh[tt]" }
    };
    unsigned long   Flags;
    ivec            ddxyz;
    gmx_bool        bStartFromCpt;
    char          **multidir = nullptr;

    commRec_ = init_commrec();

    unsigned long PCA_Flags = PCA_CAN_SET_DEFFNM;
    // With -multi or -multidir, the file names are going to get processed
    // further (or the working directory changed), so we can't check for their
    // existence during parsing.  It isn't useful to do any completion based on
    // file system contents, either.
    if (is_multisim_option_set(argc, argv))
    {
        PCA_Flags |= PCA_DISABLE_INPUT_FILE_CHECKING;
    }

    /* Comment this in to do fexist calls only on master
     * works not with rerun or tables at the moment
     * also comment out the version of init_forcerec in md.c
     * with NULL instead of opt2fn
     */
    /*
       if (!MASTER(cr))
       {
       PCA_Flags |= PCA_NOT_READ_NODE;
       }
     */

    if (!parse_common_args(&argc, argv, PCA_Flags, nfile, fnm, asize(pa), pa,
                           asize(desc), desc, 0, nullptr, &oenv_))
    {
        sfree(commRec_);
        return;
    }
    ddParams_->ddRankOrder = nenum(ddrank_opt);

    // Handle option that parses GPU ids, which could be in an
    // environment variable, so that there is a way to customize it
    // when using MPI in heterogeneous contexts.
    {
        char *env = getenv("GMX_GPU_ID");
        if (env != nullptr && hw_opt.gpu_opt.gpu_id != nullptr)
        {
            gmx_fatal(FARGS, "GMX_GPU_ID and -gpu_id can not be used at the same time");
        }
        if (env != nullptr)
        {
            hw_opt.gpu_opt.gpu_id = env;
        }
    }

    hw_opt.thread_affinity = nenum(thread_aff_opt);

    /* now check the -multi and -multidir option */
    if (opt2bSet("-multidir", nfile, fnm))
    {
        if (nmultisim > 0)
        {
            gmx_fatal(FARGS, "mdrun -multi and -multidir options are mutually exclusive.");
        }
        nmultisim = opt2fns(&multidir, "-multidir", nfile, fnm);
    }


    if (replExParams.exchangeInterval != 0 && nmultisim < 2)
    {
        gmx_fatal(FARGS, "Need at least two replicas for replica exchange (option -multi)");
    }

    if (replExParams.numExchanges < 0)
    {
        gmx_fatal(FARGS, "Replica exchange number of exchanges needs to be positive");
    }

    if (nmultisim >= 1)
    {
#if !GMX_THREAD_MPI
        init_multisystem(commRec_, nmultisim, multidir, nfile, fnm);
#else
        gmx_fatal(FARGS, "mdrun -multi or -multidir are not supported with the thread-MPI library. "
                  "Please compile GROMACS with a proper external MPI library.");
#endif
    }

    if (!opt2bSet("-cpi", nfile, fnm))
    {
        // If we are not starting from a checkpoint we never allow files to be appended
        // to, since that has caused a ton of strange behaviour and bugs in the past.
        if (opt2parg_bSet("-append", asize(pa), pa))
        {
            // If the user explicitly used the -append option, explain that it is not possible.
            gmx_fatal(FARGS, "GROMACS can only append to files when restarting from a checkpoint.");
        }
        else
        {
            // If the user did not say anything explicit, just disable appending.
            bTryToAppendFiles = FALSE;
        }
    }

    handleRestart(commRec_, bTryToAppendFiles, nfile, fnm, &doAppendFiles_, &bStartFromCpt);

    Flags = opt2bSet("-rerun", nfile, fnm) ? MD_RERUN : 0;
    Flags = Flags | (bDDBondCheck  ? MD_DDBONDCHECK  : 0);
    Flags = Flags | (bDDBondComm   ? MD_DDBONDCOMM   : 0);
    Flags = Flags | (bTunePME      ? MD_TUNEPME      : 0);
    Flags = Flags | (bConfout      ? MD_CONFOUT      : 0);
    Flags = Flags | (bRerunVSite   ? MD_RERUN_VSITE  : 0);
    Flags = Flags | (bReproducible ? MD_REPRODUCIBLE : 0);
    Flags = Flags | (doAppendFiles_  ? MD_APPENDFILES  : 0);
    Flags = Flags | (opt2parg_bSet("-append", asize(pa), pa) ? MD_APPENDFILESSET : 0);
    Flags = Flags | (bKeepAndNumCPT ? MD_KEEPANDNUMCPT : 0);
    Flags = Flags | (bStartFromCpt ? MD_STARTFROMCPT : 0);
    Flags = Flags | (bResetCountersHalfWay ? MD_RESETCOUNTERSHALFWAY : 0);
    Flags = Flags | (opt2parg_bSet("-ntomp", asize(pa), pa) ? MD_NTOMPSET : 0);
    Flags = Flags | (bIMDwait      ? MD_IMDWAIT      : 0);
    Flags = Flags | (bIMDterm      ? MD_IMDTERM      : 0);
    Flags = Flags | (bIMDpull      ? MD_IMDPULL      : 0);

    /* We postpone opening the log file if we are appending, so we can
       first truncate the old log file and append to the correct position
       there instead.  */
    if (MASTER(commRec_) && !doAppendFiles_)
    {
        gmx_log_open(ftp2fn(efLOG, nfile, fnm), commRec_,
                     Flags & MD_APPENDFILES, &fplog_);
    }
    else
    {
        fplog_ = nullptr;
    }

    ddxyz[XX] = (int)(realddxyz[XX] + 0.5);
    ddxyz[YY] = (int)(realddxyz[YY] + 0.5);
    ddxyz[ZZ] = (int)(realddxyz[ZZ] + 0.5);

    // auto args = make_mdrunner_arglist(&hw_opt, fplog, cr, nfile, fnm, oenv, bVerbose,
    //                                        nstglobalcomm, ddxyz, dd_rank_order, npme, rdd, rconstr,
    //                                        dddlb_opt[0], dlb_scale, ddcsx, ddcsy, ddcsz,
    //                                        nbpu_opt[0], nstlist,
    //                                        nsteps, nstepout, resetstep,
    //                                        nmultisim, replExParams,
    //                                        pforce, cpt_period, max_hours, imdport, Flags);

//
// int mdrunner(gmx_hw_opt_t *hw_opt,
//              FILE *fplog, t_commrec *cr, int nfile,
//              const t_filenm fnm[], const gmx_output_env_t *oenv, gmx_bool bVerbose,
//              int nstglobalcomm,
//              ivec ddxyz, int dd_rank_order, int npme, real rdd, real rconstr,
//              const char *dddlb_opt, real dlb_scale,
//              const char *ddcsx, const char *ddcsy, const char *ddcsz,
//              const char *nbpu_opt, int nstlist_cmdline,
//              gmx_int64_t nsteps_cmdline, int nstepout, int resetstep,
//              int gmx_unused nmultisim,
//              const ReplicaExchangeParameters &replExParams,
//              real pforce, real cpt_period, real max_hours,
//              int imdport, unsigned long Flags)
// {
    gmx_bool                  bForceUseGPU, bTryUseGPU, bRerunMD;
    matrix                    box;
    gmx_ddbox_t               ddbox = {0};
    int                       npme_major, npme_minor;
    real                      ewaldcoeff_q     = 0;
    real                      ewaldcoeff_lj    = 0;
    struct gmx_pme_t        **pmedata          = nullptr;
    int                       nChargePerturbed = -1, nTypePerturbed = 0, status;
    gmx_int64_t               reset_counters;
    int                       nthreads_pme = 1;
    //gmx_hw_info_t            *hwinfo       = nullptr;
    /* The master rank decides early on bUseGPU and broadcasts this later */
    gmx_bool                  bUseGPU            = FALSE;

    /* CAUTION: threads may be started later on in this function, so
       cr doesn't reflect the final parallel state right now */
    // gmx::MDModules mdModules;
    // t_inputrec     inputrecInstance;
    // t_inputrec    *inputrec = &inputrecInstance;
    t_inputrec* inputrec = input.inputRecord_.get();

    if (Flags & MD_APPENDFILES)
    {
        fplog_ = nullptr;
    }

    doMembed_    = opt2bSet("-membed", nfile, fnm);
    bRerunMD     = (Flags & MD_RERUN);

    /* Handle GPU-related user options. Later, we check consistency
     * with things like whether support is compiled, or tMPI thread
     * count. */
    bForceUseGPU = (strncmp(nbpu_opt[0], "gpu", 3) == 0);
    bTryUseGPU   = (strncmp(nbpu_opt[0], "auto", 4) == 0) || bForceUseGPU;
    gmx_parse_gpu_ids(&hw_opt.gpu_opt);

    // Here we assume that SIMMASTER(cr) does not change even after the
    // threads are started.
    logOwner_ = gmx::compat::make_unique<gmx::LoggerOwner>(buildLogger(fplog_, commRec_));
    auto mdlog = logOwner_->logger();

    /* Detect hardware, gather information. This is an operation that is
     * global for this process (MPI rank). */
    hardwareInfo_.reset(gmx_detect_hardware(mdlog, commRec_, bTryUseGPU));
    // hwinfo is now a non-owning pointer that will get passed along.
    // TODO: resources like these could be owned by the Context object that
    // will be guaranteed to outlive the execution.
    auto hwinfo = hardwareInfo_.get();

    gmx_print_detected_hardware(fplog_, commRec_, mdlog, hwinfo);

    if (fplog_ != nullptr)
    {
        /* Print references after all software/hardware printing */
        please_cite(fplog_, "Abraham2015");
        please_cite(fplog_, "Pall2015");
        please_cite(fplog_, "Pronk2013");
        please_cite(fplog_, "Hess2008b");
        please_cite(fplog_, "Spoel2005a");
        please_cite(fplog_, "Lindahl2001a");
        please_cite(fplog_, "Berendsen95a");
    }

    //std::unique_ptr<t_state> stateInstance = std::unique_ptr<t_state>(new t_state);
    t_state *state = input.state_.get();

    auto     mtop = input.topology_.get();

    if (SIMMASTER(commRec_))
    {
        /* Read (nearly) all data required for the simulation */
        //read_tpx_state(ftp2fn(efTPR, nfile, fnm), inputrec, state, mtop);

        if (inputrec->cutoff_scheme == ecutsVERLET)
        {
            /* Here the master rank decides if all ranks will use GPUs */
            bUseGPU = (hwinfo->gpu_info.n_dev_compatible > 0 ||
                       getenv("GMX_EMULATE_GPU") != nullptr);

            /* TODO add GPU kernels for this and replace this check by:
             * (bUseGPU && (ir->vdwtype == evdwPME &&
             *               ir->ljpme_combination_rule == eljpmeLB))
             * update the message text and the content of nbnxn_acceleration_supported.
             */
            if (bUseGPU &&
                !nbnxn_gpu_acceleration_supported(mdlog, inputrec, bRerunMD))
            {
                /* Fallback message printed by nbnxn_acceleration_supported */
                if (bForceUseGPU)
                {
                    gmx_fatal(FARGS, "GPU acceleration requested, but not supported with the given input settings");
                }
                bUseGPU = FALSE;
            }

            prepare_verlet_scheme(fplog_, commRec_,
                                  inputrec, nstlist, mtop, state->box,
                                  bUseGPU, *hwinfo->cpuInfo);
        }
        else
        {
            if (nstlist > 0)
            {
                gmx_fatal(FARGS, "Can not set nstlist with the group cut-off scheme");
            }

            if (hwinfo->gpu_info.n_dev_compatible > 0)
            {
                GMX_LOG(mdlog.warning).asParagraph().appendText(
                        "NOTE: GPU(s) found, but the current simulation can not use GPUs\n"
                        "      To use a GPU, set the mdp option: cutoff-scheme = Verlet");
            }

            if (bForceUseGPU)
            {
                gmx_fatal(FARGS, "GPU requested, but can't be used without cutoff-scheme=Verlet");
            }

#if GMX_TARGET_BGQ
            md_print_warn(commRec_, fplog_,
                          "NOTE: There is no SIMD implementation of the group scheme kernels on\n"
                          "      BlueGene/Q. You will observe better performance from using the\n"
                          "      Verlet cut-off scheme.\n");
#endif
        }
    }

    /* Check and update the hardware options for internal consistency */
    check_and_update_hw_opt_1(&hw_opt, commRec_, ddParams_->nPmeRanks);

    /* Early check for externally set process affinity. */
    gmx_check_thread_affinity_set(mdlog, commRec_,
                                  &hw_opt, hwinfo->nthreads_hw_avail, FALSE);

#if GMX_THREAD_MPI
    if (SIMMASTER(commRec_))
    {
        if (ddParams_->nPmeRanks > 0 && hw_opt.nthreads_tmpi <= 0)
        {
            gmx_fatal(FARGS, "You need to explicitly specify the number of MPI threads (-ntmpi) when using separate PME ranks");
        }

        /* Since the master knows the cut-off scheme, update hw_opt for this.
         * This is done later for normal MPI and also once more with tMPI
         * for all tMPI ranks.
         */
        check_and_update_hw_opt_2(&hw_opt, inputrec->cutoff_scheme);

        /* NOW the threads will be started: */
        hw_opt.nthreads_tmpi = get_nthreads_mpi(hwinfo,
                                                &hw_opt,
                                                inputrec, mtop,
                                                mdlog, bUseGPU,
                                                doMembed_);

        if (hw_opt.nthreads_tmpi > 1)
        {
            t_commrec *cr_old       = commRec_;
            /* now start the threads. */
            commRec_ = mdrunner_start_threads(&hw_opt,
                                              fplog_,
                                              cr_old,
                                              nfile,
                                              fnm,
                                              oenv_,
                                              bVerbose,
                                              nstglobalcomm,
                                              ddxyz,
                                              ddParams_->ddRankOrder,
                                              ddParams_->nPmeRanks,
                                              ddParams_->commDistanceMin,
                                              ddParams_->rConstraints,
                                              dddlb_opt[0],
                                              ddParams_->dlbScale,
                                              ddParams_->sizeX.c_str(),
                                              ddParams_->sizeY.c_str(),
                                              ddParams_->sizeZ.c_str(),
                                              nbpu_opt[0],
                                              nstlist,
                                              nsteps,
                                              nstepout_,
                                              resetstep,
                                              nmultisim,
                                              *replExParams_,
                                              pforce,
                                              cpt_period_,
                                              maxHours_,
                                              flags_);
            /* the main thread continues here with a new cr. We don't deallocate
               the old cr because other threads may still be reading it. */
            if (commRec_ == nullptr)
            {
                gmx_comm("Failed to spawn threads");
            }
        }
    }
#endif
    /* END OF CAUTION: cr is now reliable */

    if (PAR(commRec_))
    {
        /* now broadcast everything to the non-master nodes/threads: */
        init_parallel(commRec_, inputrec, mtop);

        /* The master rank decided on the use of GPUs,
         * broadcast this information to all ranks.
         */
        gmx_bcast_sim(sizeof(bUseGPU), &bUseGPU, commRec_);
    }
    // TODO: Error handling
    mdModules_.assignOptionsToModules(*inputrec->params, nullptr);

    if (fplog_ != nullptr)
    {
        pr_inputrec(fplog_, 0, "Input Parameters", inputrec, FALSE);
        fprintf(fplog_, "\n");
    }

    /* now make sure the state is initialized and propagated */
    set_state_entries(state, inputrec);

    /* A parallel command line option consistency check that we can
       only do after any threads have started. */
    if (!PAR(commRec_) &&
        (ddxyz[XX] > 1 || ddxyz[YY] > 1 || ddxyz[ZZ] > 1 || ddParams_->nPmeRanks > 0))
    {
        gmx_fatal(FARGS,
                  "The -dd or -npme option request a parallel simulation, "
#if !GMX_MPI
                  "but %s was compiled without threads or MPI enabled"
#else
#if GMX_THREAD_MPI
                  "but the number of MPI-threads (option -ntmpi) is not set or is 1"
#else
                  "but %s was not started through mpirun/mpiexec or only one rank was requested through mpirun/mpiexec"
#endif
#endif
                  , output_env_get_program_display_name(oenv_)
                  );
    }

    if (bRerunMD &&
        (EI_ENERGY_MINIMIZATION(inputrec->eI) || eiNM == inputrec->eI))
    {
        gmx_fatal(FARGS, "The .mdp file specified an energy mininization or normal mode algorithm, and these are not compatible with mdrun -rerun");
    }

    if (can_use_allvsall(inputrec, TRUE, commRec_, fplog_) && DOMAINDECOMP(commRec_))
    {
        gmx_fatal(FARGS, "All-vs-all loops do not work with domain decomposition, use a single MPI rank");
    }

    if (!(EEL_PME(inputrec->coulombtype) || EVDW_PME(inputrec->vdwtype)))
    {
        if (ddParams_->nPmeRanks > 0)
        {
            gmx_fatal_collective(FARGS, commRec_->mpi_comm_mysim, MASTER(commRec_),
                                 "PME-only ranks are requested, but the system does not use PME for electrostatics or LJ");
        }

        ddParams_->nPmeRanks = 0;
    }

    if (bUseGPU && ddParams_->nPmeRanks < 0)
    {
        /* With GPUs we don't automatically use PME-only ranks. PME ranks can
         * improve performance with many threads per GPU, since our OpenMP
         * scaling is bad, but it's difficult to automate the setup.
         */
        ddParams_->nPmeRanks = 0;
    }

#ifdef GMX_FAHCORE
    if (MASTER(commRec_))
    {
        fcRegisterSteps(inputrec->nsteps, inputrec->init_step);
    }
#endif

    /* NMR restraints must be initialized before load_checkpoint,
     * since with time averaging the history is added to t_state.
     * For proper consistency check we therefore need to extend
     * t_state here.
     * So the PME-only nodes (if present) will also initialize
     * the distance restraints.
     */

    /* This needs to be called before read_checkpoint to extend the state */

    init_disres(fplog_, mtop, inputrec, commRec_, forceCalcData_.get(), state, replExParams.exchangeInterval > 0);

    init_orires(fplog_, mtop, as_rvec_array(state->x.data()), inputrec, commRec_, &(forceCalcData_->orires),
                state);

    if (inputrecDeform(inputrec))
    {
        /* Store the deform reference box before reading the checkpoint */
        if (SIMMASTER(commRec_))
        {
            copy_mat(state->box, box);
        }
        if (PAR(commRec_))
        {
            gmx_bcast(sizeof(box), box, commRec_);
        }
        /* Because we do not have the update struct available yet
         * in which the reference values should be stored,
         * we store them temporarily in static variables.
         * This should be thread safe, since they are only written once
         * and with identical values.
         */
        tMPI_Thread_mutex_lock(&deform_init_box_mutex);
        deform_init_init_step_tpx = inputrec->init_step;
        copy_mat(box, deform_init_box_tpx);
        tMPI_Thread_mutex_unlock(&deform_init_box_mutex);
    }

    if (Flags & MD_STARTFROMCPT)
    {
        /* Check if checkpoint file exists before doing continuation.
         * This way we can use identical input options for the first and subsequent runs...
         */
        gmx_bool bReadEkin;

        load_checkpoint(opt2fn_master("-cpi", nfile, fnm, commRec_), &fplog_,
                        commRec_, ddxyz, &ddParams_->nPmeRanks,
                        inputrec, state, &bReadEkin, observablesHistory_.get(),
                        (Flags & MD_APPENDFILES),
                        (Flags & MD_APPENDFILESSET),
                        (Flags & MD_REPRODUCIBLE));

        if (bReadEkin)
        {
            Flags |= MD_READ_EKIN;
        }
    }

    if (SIMMASTER(commRec_) && (Flags & MD_APPENDFILES))
    {
        gmx_log_open(ftp2fn(efLOG, nfile, fnm), commRec_,
                     Flags, &fplog_);
        logOwner_ = gmx::compat::make_unique<gmx::LoggerOwner>(buildLogger(fplog_, nullptr));
        mdlog     = logOwner_->logger();
    }

    /* override nsteps with value from cmdline */
    override_nsteps_cmdline(mdlog, nsteps, inputrec);

    if (SIMMASTER(commRec_))
    {
        copy_mat(state->box, box);
    }

    if (PAR(commRec_))
    {
        gmx_bcast(sizeof(box), box, commRec_);
    }

    if (PAR(commRec_) && !(EI_TPI(inputrec->eI) ||
                           inputrec->eI == eiNM))
    {
        commRec_->dd = init_domain_decomposition(fplog_,
                                                 commRec_,
                                                 *ddParams_,
                                                 ddxyz,
                                                 mtop, inputrec,
                                                 box, as_rvec_array(state->x.data()),
                                                 &ddbox, &npme_major, &npme_minor);
    }
    else
    {
        /* PME, if used, is done on all nodes with 1D decomposition */
        commRec_->npmenodes = 0;
        commRec_->duty      = (DUTY_PP | DUTY_PME);
        npme_major          = 1;
        npme_minor          = 1;

        if (inputrec->ePBC == epbcSCREW)
        {
            gmx_fatal(FARGS,
                      "pbc=%s is only implemented with domain decomposition",
                      epbc_names[inputrec->ePBC]);
        }
    }

    if (PAR(commRec_))
    {
        /* After possible communicator splitting in make_dd_communicators.
         * we can set up the intra/inter node communication.
         */
        gmx_setup_nodecomm(fplog_, commRec_);
    }

    /* Initialize per-physical-node MPI process/thread ID and counters. */
    gmx_init_intranode_counters(commRec_);
#if GMX_MPI
    if (MULTISIM(commRec_))
    {
        GMX_LOG(mdlog.warning).asParagraph().appendTextFormatted(
                "This is simulation %d out of %d running as a composite GROMACS\n"
                "multi-simulation job. Setup for this simulation:\n",
                commRec_->ms->sim, commRec_->ms->nsim);
    }
    GMX_LOG(mdlog.warning).appendTextFormatted(
            "Using %d MPI %s\n",
            commRec_->nnodes,
#if GMX_THREAD_MPI
            commRec_->nnodes == 1 ? "thread" : "threads"
#else
            commRec_->nnodes == 1 ? "process" : "processes"
#endif
            );
    fflush(stderr);
#endif

    /* Check and update hw_opt for the cut-off scheme */
    check_and_update_hw_opt_2(&hw_opt, inputrec->cutoff_scheme);

    /* Check and update hw_opt for the number of MPI ranks */
    check_and_update_hw_opt_3(&hw_opt);

    gmx_omp_nthreads_init(mdlog, commRec_,
                          hwinfo->nthreads_hw_avail,
                          hw_opt.nthreads_omp,
                          hw_opt.nthreads_omp_pme,
                          (commRec_->duty & DUTY_PP) == 0,
                          inputrec->cutoff_scheme == ecutsVERLET);

#ifndef NDEBUG
    if (EI_TPI(inputrec->eI) &&
        inputrec->cutoff_scheme == ecutsVERLET)
    {
        gmx_feenableexcept();
    }
#endif

    // TODO: does hasUserSetGpuIds need to write?
    bool userSetGpuIds = hasUserSetGpuIds(&hw_opt.gpu_opt);

    if (bUseGPU)
    {
        /* Select GPU id's to use */
        // TODO: does this function need to write to these?
        gmx_select_rank_gpu_ids(mdlog, commRec_, &hwinfo->gpu_info, bForceUseGPU,
                                userSetGpuIds, &hw_opt.gpu_opt);
    }
    else
    {
        /* Ignore (potentially) manually selected GPUs */
        hw_opt.gpu_opt.n_dev_use = 0;
    }

    /* check consistency across ranks of things like SIMD
     * support and number of GPUs selected */
    gmx_check_hw_runconf_consistency(mdlog, hwinfo, commRec_, &hw_opt, userSetGpuIds, bUseGPU);

    /* Now that we know the setup is consistent, check for efficiency */
    check_resource_division_efficiency(hwinfo, &hw_opt, hw_opt.gpu_opt.n_dev_use, Flags & MD_NTOMPSET,
                                       commRec_, mdlog);

    if (DOMAINDECOMP(commRec_))
    {
        /* When we share GPUs over ranks, we need to know this for the DLB */
        dd_setup_dlb_resource_sharing(commRec_, hwinfo, &hw_opt);
    }

    /* getting number of PP/PME threads
       PME: env variable should be read only on one node to make sure it is
       identical everywhere;
     */
    nthreads_pme = gmx_omp_nthreads_get(emntPME);

    wallCycle_ = wallcycle_init(fplog_, resetstep, commRec_);

    if (PAR(commRec_))
    {
        /* Master synchronizes its value of reset_counters with all nodes
         * including PME only nodes */
        reset_counters = wcycle_get_reset_counters(wallCycle_);
        gmx_bcast_sim(sizeof(reset_counters), &reset_counters, commRec_);
        wcycle_set_reset_counters(wallCycle_, reset_counters);
    }

    // Membrane embedding must be initialized before we call init_forcerec()
    if (doMembed_)
    {
        if (MASTER(commRec_))
        {
            fprintf(stderr, "Initializing membed");
        }
        /* Note that membed cannot work in parallel because mtop is
         * changed here. Fix this if we ever want to make it run with
         * multiple ranks. */
        membed_ = std::unique_ptr < gmx_membed_t, std::function < void(gmx_membed_t*)>>(init_membed(fplog_, nfile, fnm, mtop, inputrec, state, commRec_, &cpt_period_), free_membed);
    }

    if (commRec_->duty & DUTY_PP)
    {
        bcast_state(commRec_, state);

        /* Initiate forcerecord */
        forceRecord_          = std::unique_ptr<t_forcerec>(mk_forcerec());
        forceRecord_->hwinfo  = hwinfo;
        //
        // Note! fr now has a non-owning pointer to a member of hw_opt.
        // TODO: Check lifetime assumptions.
        //
        forceRecord_->gpu_opt = &hw_opt.gpu_opt;
        init_forcerec(fplog_, mdlog, forceRecord_.get(), forceCalcData_.get(), mdModules_.forceProvider(),
                      inputrec, mtop, commRec_, box,
                      opt2fn("-table", nfile, fnm),
                      opt2fn("-tablep", nfile, fnm),
                      getFilenm("-tableb", nfile, fnm),
                      nbpu_opt[0],
                      FALSE,
                      pforce);

        /* Initialize QM-MM */
        if (forceRecord_->bQMMM)
        {
            init_QMMMrec(commRec_, mtop, inputrec, forceRecord_.get());
        }

        /* Initialize the mdatoms structure.
         * mdatoms is not filled with atom data,
         * as this can not be done now with domain decomposition.
         */
        mdAtoms_.reset(init_mdatoms(fplog_, mtop, inputrec->efep != efepNO));

        /* Initialize the virtual site communication */
        vSite_.reset(init_vsite(mtop, commRec_, FALSE));

        calc_shifts(box, forceRecord_->shift_vec);

        /* With periodic molecules the charge groups should be whole at start up
         * and the virtual sites should not be far from their proper positions.
         */
        if (!inputrec->bContinuation && MASTER(commRec_) &&
            !(inputrec->ePBC != epbcNONE && inputrec->bPeriodicMols))
        {
            /* Make molecules whole at start of run */
            if (forceRecord_->ePBC != epbcNONE)
            {
                do_pbc_first_mtop(fplog_, inputrec->ePBC, box, mtop, as_rvec_array(state->x.data()));
            }
            if (vSite_ != nullptr)
            {
                /* Correct initial vsite positions are required
                 * for the initial distribution in the domain decomposition
                 * and for the initial shell prediction.
                 */
                construct_vsites_mtop(vSite_.get(), mtop, as_rvec_array(state->x.data()));
            }
        }

        if (EEL_PME(forceRecord_->eeltype) || EVDW_PME(forceRecord_->vdwtype))
        {
            ewaldcoeff_q  = forceRecord_->ewaldcoeff_q;
            ewaldcoeff_lj = forceRecord_->ewaldcoeff_lj;
            pmedata       = &forceRecord_->pmedata;
        }
        else
        {
            pmedata = nullptr;
        }
    }
    else
    {
        /* This is a PME only node */

        /* We don't need the state */
        input.state_.reset();
        state         = nullptr;

        ewaldcoeff_q  = calc_ewaldcoeff_q(inputrec->rcoulomb, inputrec->ewald_rtol);
        ewaldcoeff_lj = calc_ewaldcoeff_lj(inputrec->rvdw, inputrec->ewald_rtol_lj);
        snew(pmedata, 1);
    }

    if (hw_opt.thread_affinity != threadaffOFF)
    {
        /* Before setting affinity, check whether the affinity has changed
         * - which indicates that probably the OpenMP library has changed it
         * since we first checked).
         */
        gmx_check_thread_affinity_set(mdlog, commRec_,
                                      &hw_opt, hwinfo->nthreads_hw_avail, TRUE);

        int nthread_local;
        /* threads on this MPI process or TMPI thread */
        if (commRec_->duty & DUTY_PP)
        {
            nthread_local = gmx_omp_nthreads_get(emntNonbonded);
        }
        else
        {
            nthread_local = gmx_omp_nthreads_get(emntPME);
        }

        /* Set the CPU affinity */
        // TODO: Note some of these functions take pointers to const. E.g. the
        // addresses are not output params...
        gmx_set_thread_affinity(mdlog, commRec_, &hw_opt, *hwinfo->hardwareTopology,
                                nthread_local, nullptr);
    }

    /* Initiate PME if necessary,
     * either on all nodes or on dedicated PME nodes only. */
    if (EEL_PME(inputrec->coulombtype) || EVDW_PME(inputrec->vdwtype))
    {
        if (mdAtoms_ != nullptr)
        {
            nChargePerturbed = mdAtoms_->nChargePerturbed;
            if (EVDW_PME(inputrec->vdwtype))
            {
                nTypePerturbed   = mdAtoms_->nTypePerturbed;
            }
        }
        if (commRec_->npmenodes > 0)
        {
            /* The PME only nodes need to know nChargePerturbed(FEP on Q) and nTypePerturbed(FEP on LJ)*/
            gmx_bcast_sim(sizeof(nChargePerturbed), &nChargePerturbed, commRec_);
            gmx_bcast_sim(sizeof(nTypePerturbed), &nTypePerturbed, commRec_);
        }

        if (commRec_->duty & DUTY_PME)
        {
            try
            {
                status = gmx_pme_init(pmedata, commRec_, npme_major, npme_minor, inputrec,
                                      mtop ? mtop->natoms : 0, nChargePerturbed, nTypePerturbed,
                                      (Flags & MD_REPRODUCIBLE),
                                      ewaldcoeff_q, ewaldcoeff_lj,
                                      nthreads_pme);
            }
            GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
            if (status != 0)
            {
                gmx_fatal(FARGS, "Error %d initializing PME", status);
            }
        }
    }

    // Stash electrostatics parameters.
    ewaldcoeff_q_  = ewaldcoeff_q;
    ewaldcoeff_lj_ = ewaldcoeff_lj;
    // Capture pmedata if it exists
    if (pmedata != nullptr)
    {
        pmeData_ = std::unique_ptr<gmx_pme_t, std::function<void(gmx_pme_t*)> >(*pmedata, &gmx_pme_destroy);
    }

    if (EI_DYNAMICS(inputrec->eI))
    {
        /* Turn on signal handling on all nodes */
        /*
         * (A user signal from the PME nodes (if any)
         * is communicated to the PP nodes.
         */
        signal_handler_install();
    }

    nstglobalcomm_ = nstglobalcomm;
}

int RunnerImpl::run()
{
    int  rc    = 0; // initialize return code
    auto mdlog = logOwner_->logger();
    auto mtop  = input_->topology_.get();

    if (commRec_->duty & DUTY_PP)
    {
        /* Assumes uniform use of the number of OpenMP threads */
        walltimeAccounting_ = walltime_accounting_init(gmx_omp_nthreads_get(emntDefault));

        if (input_->inputRecord_->bPull)
        {
            /* Initialize pull code */
            input_->inputRecord_->pull_work =
                init_pull(fplog_,
                          input_->inputRecord_->pull,
                          input_->inputRecord_.get(),
                          fnm_.size(),
                          fnm_.data(),
                          mtop,
                          commRec_,
                          oenv_,
                          input_->inputRecord_->fepvals->init_lambda,
                          EI_DYNAMICS(input_->inputRecord_->eI) && MASTER(commRec_), flags_);
        }

        if (input_->inputRecord_->bRot)
        {
            /* Initialize enforced rotation code */
            init_rot(fplog_,
                     input_->inputRecord_.get(),
                     fnm_.size(),
                     fnm_.data(),
                     commRec_,
                     as_rvec_array(input_->state_->x.data()),
                     input_->state_->box,
                     mtop,
                     oenv_,
                     verbose_,
                     flags_);
        }


        /* Let init_constraints know whether we have essential dynamics constraints.
         * TODO: inputrec should tell us whether we use an algorithm, not a file option or the checkpoint
         */
        bool doEdsam = (opt2fn_null("-ei", fnm_.size(), fnm_.data()) != nullptr || observablesHistory_->edsamHistory);

        auto constraints = init_constraints(fplog_,
                                            input_->topology_.get(),
                                            input_->inputRecord_.get(),
                                            doEdsam,
                                            commRec_);
        if (DOMAINDECOMP(commRec_))
        {
            GMX_RELEASE_ASSERT(forceRecord_.get(), "fr was NULL while cr->duty was DUTY_PP");
            /* This call is not included in init_domain_decomposition mainly
             * because forceRecord_->cginfo_mb is set later.
             */
            dd_init_bondeds(fplog_, commRec_->dd, mtop, vSite_.get(), input_->inputRecord_.get(),
                            flags_ & MD_DDBONDCHECK, forceRecord_->cginfo_mb);
        }

        auto mdEngine = my_integrator(input_->inputRecord_->eI);

        rc =  mdEngine(
                    fplog_,
                    commRec_,
                    mdlog,
                    fnm_.size(),
                    fnm_.data(),
                    oenv_,
                    verbose_,
                    nstglobalcomm_,
                    vSite_.get(),
                    constraints,
                    nstepout_,
                    mdModules_.outputProvider(),
                    input_->inputRecord_.get(),
                    input_->topology_.get(),
                    forceCalcData_.get(),
                    input_->state_.get(),
                    observablesHistory_.get(),
                    mdAtoms_.get(),
                    nrNonBonded_.get(),
                    wallCycle_,
                    forceRecord_.get(),
                    *replExParams_,
                    membed_.get(),
                    cpt_period_,
                    maxHours_,
                    imdport_,
                    flags_,
                    walltimeAccounting_);
        if (input_->inputRecord_.get()->bRot)
        {
            finish_rot(input_->inputRecord_.get()->rot);
        }

        if (input_->inputRecord_.get()->bPull)
        {
            finish_pull(input_->inputRecord_.get()->pull_work);
        }
    }
    else // not an active particle-pair rank:
    {
        GMX_RELEASE_ASSERT(pmeData_ != nullptr, "pmedata was NULL while cr->duty was not DUTY_PP");
        /* do PME only */
        gmx_pmeonly(pmeData_.get(),
                    commRec_,
                    nrNonBonded_.get(),
                    wallCycle_,
                    walltimeAccounting_,
                    ewaldcoeff_q_,
                    ewaldcoeff_lj_,
                    input_->inputRecord_.get());
    }

    wallcycle_stop(wallCycle_, ewcRUN);

    // Use initial step MD input to track current step. Assume we took the steps we said we would.
    input_->inputRecord_->init_step += input_->inputRecord_->nsteps;

    return rc;
}

int RunnerImpl::run(unsigned int numSteps)
{
    input_->inputRecord_->nsteps = numSteps;
    return run();
}

std::shared_ptr < std::vector < std::array<real, 3>>> RunnerImpl::getX() const
{
    // t_state.x is PaddedRVecVector, which is ultimately a std::vector<real[3]>, but at this point let's
    // emphasize safety over performance and not assume layout or lifetime.
    // Note that state is nullptr on PME nodes.
    if (input_->state_ == nullptr)
    {
        return nullptr;
    }
    auto        &x         = input_->state_->x;
    auto         length    = x.size();
    auto         positions = std::make_shared < std::vector < std::array<real, 3>>>(length);
    unsigned int i {
        0
    };
    for (auto &position : *positions)
    {
        position[0] = x[i][0];
        position[1] = x[i][1];
        position[2] = x[i][2];
        i++;
    }
    return positions;
}


int RunnerImpl::close()
{
    if (!initialized_)
    {
        return 0;
    }
    std::cout << "Shutting down runner" << std::endl;
    initialized_ = false;
    auto              mdlog = logOwner_->logger();

    auto              hwinfo  = hardwareInfo_.release(); // hwinfo is explicitly deleted elsewhere.
    struct gmx_pme_t *pmedata = pmeData_.release();      // *pmedata is explicitly deleted elsewhere.

    /* Finish up, write some stuff
     * if rerunMD, don't write last frame again
     */
    finish_run(fplog_,
               mdlog,
               commRec_,
               input_->inputRecord_.get(),
               nrNonBonded_.get(),
               wallCycle_, // presumably the memory is freed here...
               walltimeAccounting_,
               forceRecord_ != nullptr ? forceRecord_->nbv : nullptr,
               EI_DYNAMICS(input_->inputRecord_.get()->eI) && !MULTISIM(commRec_));

    // Free PME data
    if (pmedata != nullptr)
    {
        // gmx_pme_destroy _does_ do sfree(pme)
        gmx_pme_destroy(pmedata); // TODO: pmedata is always a single element list, refactor
        pmedata = nullptr;
    }

    /* Free GPU memory and context */
    free_gpu_resources(forceRecord_.get(),
                       commRec_,
                       &hwinfo->gpu_info,
                       forceRecord_ != nullptr ? forceRecord_->gpu_opt : nullptr);

    // Not sure if the order matters, but this is where the deleter is called in the earlier code.
    membed_.reset();

    // Frees hwinfo_g after count of hwinfo pointers goes to 0.
    gmx_hardware_info_free(hwinfo);

    /* Does what it says */
    print_date_and_time(fplog_, commRec_->nodeid, "Finished mdrun", gmx_gettime());
    walltime_accounting_destroy(walltimeAccounting_);

    /* Close logfile already here if we were appending to it */
    if (MASTER(commRec_) && (flags_ & MD_APPENDFILES))
    {
        gmx_log_close(fplog_);
    }

    auto rc = (int)gmx_get_stop_condition();

    #if GMX_THREAD_MPI
    /* we need to join all threads. The sub-threads join when they
       exit this function, but the master thread needs to be told to
       wait for that. */
    if (PAR(commRec_) && MASTER(commRec_))
    {
        tMPI_Finalize();
    }
    #endif

    /* Log file has to be closed in mdrunner if we are appending to it
       (fplog not set here) */
    if (MASTER(commRec_) && !doAppendFiles_)
    {
        gmx_log_close(fplog_);
    }
    return rc;
}

RunnerImpl::~RunnerImpl()
{
    if (initialized_)
    {
        this->close();
    }
};

} // end namespace gmxapi
