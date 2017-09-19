/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */

/*! \libinternal \file
 *
 *
 * \brief
 * This file contains datatypes and function declarations for internal
   use in the pull code.
 *
 * \author Berk Hess
 *
 * \inlibraryapi
 */

#ifndef GMX_PULLING_PULL_INTERNAL_H
#define GMX_PULLING_PULL_INTERNAL_H

#include "config.h"

#include "gromacs/mdtypes/pull-params.h"
#include "gromacs/utility/gmxmpi.h"

/*! \cond INTERNAL */

/*! \brief Determines up to what local atom count a pull group gets processed single-threaded.
 *
 * We set this limit to 1 with debug to catch bugs.
 * On Haswell with GCC 5 the cross-over point is around 400 atoms,
 * independent of thread count and hyper-threading.
 */
#ifdef NDEBUG
static const int c_pullMaxNumLocalAtomsSingleThreaded = 100;
#else
static const int c_pullMaxNumLocalAtomsSingleThreaded = 1;
#endif

enum {
    epgrppbcNONE, epgrppbcREFAT, epgrppbcCOS
};

typedef struct
{
    t_pull_group  params;

    gmx_bool      bCalcCOM;   /* Calculate COM? Not if only used as cylinder group */
    int           epgrppbc;   /* The type of pbc for this pull group, see enum above */

    int           nat_loc;    /* Number of local pull atoms */
    int           nalloc_loc; /* Allocation size for ind_loc and weight_loc */
    int          *ind_loc;    /* Local pull indices */
    real         *weight_loc; /* Weights for the local indices */

    real          mwscale;    /* mass*weight scaling factor 1/sum w m */
    real          wscale;     /* scaling factor for the weights: sum w m/sum w w m */
    real          invtm;      /* inverse total mass of the group: 1/wscale sum w m */
    dvec         *mdw;        /* mass*gradient(weight) for atoms */
    double       *dv;         /* distance to the other group along vec */
    dvec          x;          /* center of mass before update */
    dvec          xp;         /* center of mass after update before constraining */
}
pull_group_work_t;

/*!
 * \brief State of pull coords and work
 */
typedef struct
{
    t_pull_coord  params;     /*!< \brief Pull coordinate (constant) parameters */

    double        value_ref;  /*!< \brief The reference value, usually init+rate*t, units of nm or rad */
    double        value;      /*!< \brief The current value of the coordinate, units of nm or rad */
    dvec          dr01;       /*!< \brief The direction vector of group 1 relative to group 0 */
    dvec          dr23;       /*!< \brief The direction vector of group 3 relative to group 2 */
    dvec          dr45;       /*!< \brief The direction vector of group 5 relative to group 4 */
    dvec          vec;        /*!< \brief The pull direction */
    double        vec_len;    /*!< \brief Length of vec for direction-relative */
    dvec          ffrad;      /*!< \brief conversion factor from vec to radial force */
    double        cyl_dev;    /*!< \brief The deviation from the reference position */
    double        f_scal;     /*!< \brief Scalar force for directional pulling */
    dvec          f01;        /*!< \brief Force due to the pulling/constraining for groups 0, 1 */
    dvec          f23;        /*!< \brief Force for groups 2 and 3 */
    dvec          f45;        /*!< \brief Force for groups 4 and 5 */
    dvec          planevec_m; /*!< \brief Normal of plane for groups 0, 1, 2, 3 for geometry dihedral */
    dvec          planevec_n; /*!< \brief Normal of plane for groups 2, 3, 4, 5 for geometry dihedral */

    /*! \brief For external-potential coordinates only, for checking if a provider has been registered */
    bool          bExternalPotentialProviderHasBeenRegistered;
}
pull_coord_work_t;

/*! \brief Struct for sums over (local) atoms in a pull group */
struct pull_sum_com_t {
    /* For normal weighting */
    double sum_wm;    /*!< \brief Sum of weight*mass        */
    double sum_wwm;   /*!< \brief Sum of weight*weight*mass */
    dvec   sum_wmx;   /*!< \brief Sum of weight*mass*x      */
    dvec   sum_wmxp;  /*!< \brief Sum of weight*mass*xp     */

    /*!< \brief For cosine weighting */
    double sum_cm;    /*!< \brief Sum of cos(x)*mass          */
    double sum_sm;    /*!< \brief Sum of sin(x)*mass          */
    double sum_ccm;   /*!< \brief Sum of cos(x)*cos(x)*mass   */
    double sum_csm;   /*!< \brief Sum of cos(x)*sin(x)*mass   */
    double sum_ssm;   /*!< \brief Sum of sin(x)*sin(x)*mass   */
    double sum_cmp;   /*!< \brief Sum of cos(xp)*sin(xp)*mass */
    double sum_smp;   /*!< \brief Sum of sin(xp)*sin(xp)*mass */

    /*! \brief Assure cache line size
     * 
     * Dummy data to ensure adjacent elements in an array are separated
     * by a cache line size, max 128 bytes.
     * TODO: Replace this by some automated mechanism.
     */
    int    dummy[32];
};

/*!
 * \brief communications state for pulling code
 */
typedef struct {
    gmx_bool    bParticipateAll; /*!< \brief Do all ranks always participate in pulling? */
    gmx_bool    bParticipate;    /*!< \brief Does our rank participate in pulling? */
#if GMX_MPI
    MPI_Comm    mpi_comm_com;    /*!< \brief Communicator for pulling */
#endif
    int         nparticipate;    /*!< \brief The number of ranks participating */

    gmx_int64_t setup_count;     /*!< \brief The number of decomposition calls */
    gmx_int64_t must_count;      /*!< \brief The last count our rank needed to be part */

    rvec       *rbuf;            /*!< \brief COM calculation buffer */
    dvec       *dbuf;            /*!< \brief COM calculation buffer */
    double     *dbuf_cyl;        /*!< \brief cylinder ref. groups calculation buffer */
}
pull_comm_t;

/*!
 * \brief Pull work struct
 */
struct pull_t
{
    pull_params_t      params;       /*!< \brief The pull parameters, from inputrec */

    gmx_bool           bPotential;   /*!< \brief Are there coordinates with potential? */
    gmx_bool           bConstraint;  /*!< \brief Are there constrained coordinates? */
    gmx_bool           bAngle;       /*!< \brief Are there angle geometry coordinates? */

    int                ePBC;         /*!< \brief the boundary conditions */
    int                npbcdim;      /*!< \brief do pbc in dims 0 <= dim < npbcdim */
    gmx_bool           bRefAt;       /*!< \brief do we need reference atoms for a group COM ? */
    int                cosdim;       /*!< \brief dimension for cosine weighting, -1 if none */

    int                ngroup;       /*!< \brief Number of pull groups */
    int                ncoord;       /*!< \brief Number of pull coordinates */
    pull_group_work_t *group;        /*!< \brief The pull group param and work data */
    pull_group_work_t *dyna;         /*!< \brief Dynamic groups for geom=cylinder */
    pull_coord_work_t *coord;        /*!< \brief The pull group param and work data */

    gmx_bool           bCylinder;    /*!< \brief Is group 0 a cylinder group? */

    gmx_bool           bSetPBCatoms; /*!< \brief Do we need to set x_pbc for the groups? */

    int                nthreads;     /*!< \brief Number of threads used by the pull code */
    pull_sum_com_t    *sum_com;      /*!< \brief Work array for summing for COM, 1 entry per thread */

    pull_comm_t        comm;         /*!< \brief Communication parameters, communicator and buffers */

    FILE              *out_x;        /*!< \brief Output file for pull data */
    FILE              *out_f;        /*!< \brief Output file for pull data */

    /*! \brief The number of coordinates using an external potential */
    int                numCoordinatesWithExternalPotential;
    /*! \brief Counter for checking external potential registration */
    int                numUnregisteredExternalPotentials;
    /*! \brief counter... */
    int                numExternalPotentialsStillToBeAppliedThisStep;
};

/*! \endcond */

#endif
