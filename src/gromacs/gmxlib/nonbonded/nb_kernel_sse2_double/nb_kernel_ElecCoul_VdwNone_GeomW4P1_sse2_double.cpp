/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2017,2018, by the GROMACS development team, led by
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
/*
 * Note: this file was generated by the GROMACS sse2_double kernel generator.
 */
#include "gmxpre.h"

#include "config.h"

#include <math.h>

#include "../nb_kernel.h"
#include "gromacs/gmxlib/nrnb.h"

#include "kernelutil_x86_sse2_double.h"

/*
 * Gromacs nonbonded kernel:   nb_kernel_ElecCoul_VdwNone_GeomW4P1_VF_sse2_double
 * Electrostatics interaction: Coulomb
 * VdW interaction:            None
 * Geometry:                   Water4-Particle
 * Calculate force/pot:        PotentialAndForce
 */
void
nb_kernel_ElecCoul_VdwNone_GeomW4P1_VF_sse2_double
                    (t_nblist                    * gmx_restrict       nlist,
                     rvec                        * gmx_restrict          xx,
                     rvec                        * gmx_restrict          ff,
                     struct t_forcerec           * gmx_restrict          fr,
                     t_mdatoms                   * gmx_restrict     mdatoms,
                     nb_kernel_data_t gmx_unused * gmx_restrict kernel_data,
                     t_nrnb                      * gmx_restrict        nrnb)
{
    /* Suffixes 0,1,2,3 refer to particle indices for waters in the inner or outer loop, or
     * just 0 for non-waters.
     * Suffixes A,B refer to j loop unrolling done with SSE double precision, e.g. for the two different
     * jnr indices corresponding to data put in the four positions in the SIMD register.
     */
    int              i_shift_offset,i_coord_offset,outeriter,inneriter;
    int              j_index_start,j_index_end,jidx,nri,inr,ggid,iidx;
    int              jnrA,jnrB;
    int              j_coord_offsetA,j_coord_offsetB;
    int              *iinr,*jindex,*jjnr,*shiftidx,*gid;
    real             rcutoff_scalar;
    real             *shiftvec,*fshift,*x,*f;
    __m128d          tx,ty,tz,fscal,rcutoff,rcutoff2,jidxall;
    int              vdwioffset1;
    __m128d          ix1,iy1,iz1,fix1,fiy1,fiz1,iq1,isai1;
    int              vdwioffset2;
    __m128d          ix2,iy2,iz2,fix2,fiy2,fiz2,iq2,isai2;
    int              vdwioffset3;
    __m128d          ix3,iy3,iz3,fix3,fiy3,fiz3,iq3,isai3;
    int              vdwjidx0A,vdwjidx0B;
    __m128d          jx0,jy0,jz0,fjx0,fjy0,fjz0,jq0,isaj0;
    __m128d          dx10,dy10,dz10,rsq10,rinv10,rinvsq10,r10,qq10,c6_10,c12_10;
    __m128d          dx20,dy20,dz20,rsq20,rinv20,rinvsq20,r20,qq20,c6_20,c12_20;
    __m128d          dx30,dy30,dz30,rsq30,rinv30,rinvsq30,r30,qq30,c6_30,c12_30;
    __m128d          velec,felec,velecsum,facel,crf,krf,krf2;
    real             *charge;
    __m128d          dummy_mask,cutoff_mask;
    __m128d          signbit   = gmx_mm_castsi128_pd( _mm_set_epi32(0x80000000,0x00000000,0x80000000,0x00000000) );
    __m128d          one     = _mm_set1_pd(1.0);
    __m128d          two     = _mm_set1_pd(2.0);
    x                = xx[0];
    f                = ff[0];

    nri              = nlist->nri;
    iinr             = nlist->iinr;
    jindex           = nlist->jindex;
    jjnr             = nlist->jjnr;
    shiftidx         = nlist->shift;
    gid              = nlist->gid;
    shiftvec         = fr->shift_vec[0];
    fshift           = fr->fshift[0];
    facel            = _mm_set1_pd(fr->ic->epsfac);
    charge           = mdatoms->chargeA;

    /* Setup water-specific parameters */
    inr              = nlist->iinr[0];
    iq1              = _mm_mul_pd(facel,_mm_set1_pd(charge[inr+1]));
    iq2              = _mm_mul_pd(facel,_mm_set1_pd(charge[inr+2]));
    iq3              = _mm_mul_pd(facel,_mm_set1_pd(charge[inr+3]));

    /* Avoid stupid compiler warnings */
    jnrA = jnrB = 0;
    j_coord_offsetA = 0;
    j_coord_offsetB = 0;

    outeriter        = 0;
    inneriter        = 0;

    /* Start outer loop over neighborlists */
    for(iidx=0; iidx<nri; iidx++)
    {
        /* Load shift vector for this list */
        i_shift_offset   = DIM*shiftidx[iidx];

        /* Load limits for loop over neighbors */
        j_index_start    = jindex[iidx];
        j_index_end      = jindex[iidx+1];

        /* Get outer coordinate index */
        inr              = iinr[iidx];
        i_coord_offset   = DIM*inr;

        /* Load i particle coords and add shift vector */
        gmx_mm_load_shift_and_3rvec_broadcast_pd(shiftvec+i_shift_offset,x+i_coord_offset+DIM,
                                                 &ix1,&iy1,&iz1,&ix2,&iy2,&iz2,&ix3,&iy3,&iz3);

        fix1             = _mm_setzero_pd();
        fiy1             = _mm_setzero_pd();
        fiz1             = _mm_setzero_pd();
        fix2             = _mm_setzero_pd();
        fiy2             = _mm_setzero_pd();
        fiz2             = _mm_setzero_pd();
        fix3             = _mm_setzero_pd();
        fiy3             = _mm_setzero_pd();
        fiz3             = _mm_setzero_pd();

        /* Reset potential sums */
        velecsum         = _mm_setzero_pd();

        /* Start inner kernel loop */
        for(jidx=j_index_start; jidx<j_index_end-1; jidx+=2)
        {

            /* Get j neighbor index, and coordinate index */
            jnrA             = jjnr[jidx];
            jnrB             = jjnr[jidx+1];
            j_coord_offsetA  = DIM*jnrA;
            j_coord_offsetB  = DIM*jnrB;

            /* load j atom coordinates */
            gmx_mm_load_1rvec_2ptr_swizzle_pd(x+j_coord_offsetA,x+j_coord_offsetB,
                                              &jx0,&jy0,&jz0);

            /* Calculate displacement vector */
            dx10             = _mm_sub_pd(ix1,jx0);
            dy10             = _mm_sub_pd(iy1,jy0);
            dz10             = _mm_sub_pd(iz1,jz0);
            dx20             = _mm_sub_pd(ix2,jx0);
            dy20             = _mm_sub_pd(iy2,jy0);
            dz20             = _mm_sub_pd(iz2,jz0);
            dx30             = _mm_sub_pd(ix3,jx0);
            dy30             = _mm_sub_pd(iy3,jy0);
            dz30             = _mm_sub_pd(iz3,jz0);

            /* Calculate squared distance and things based on it */
            rsq10            = gmx_mm_calc_rsq_pd(dx10,dy10,dz10);
            rsq20            = gmx_mm_calc_rsq_pd(dx20,dy20,dz20);
            rsq30            = gmx_mm_calc_rsq_pd(dx30,dy30,dz30);

            rinv10           = sse2_invsqrt_d(rsq10);
            rinv20           = sse2_invsqrt_d(rsq20);
            rinv30           = sse2_invsqrt_d(rsq30);

            rinvsq10         = _mm_mul_pd(rinv10,rinv10);
            rinvsq20         = _mm_mul_pd(rinv20,rinv20);
            rinvsq30         = _mm_mul_pd(rinv30,rinv30);

            /* Load parameters for j particles */
            jq0              = gmx_mm_load_2real_swizzle_pd(charge+jnrA+0,charge+jnrB+0);

            fjx0             = _mm_setzero_pd();
            fjy0             = _mm_setzero_pd();
            fjz0             = _mm_setzero_pd();

            /**************************
             * CALCULATE INTERACTIONS *
             **************************/

            /* Compute parameters for interactions between i and j atoms */
            qq10             = _mm_mul_pd(iq1,jq0);

            /* COULOMB ELECTROSTATICS */
            velec            = _mm_mul_pd(qq10,rinv10);
            felec            = _mm_mul_pd(velec,rinvsq10);

            /* Update potential sum for this i atom from the interaction with this j atom. */
            velecsum         = _mm_add_pd(velecsum,velec);

            fscal            = felec;

            /* Calculate temporary vectorial force */
            tx               = _mm_mul_pd(fscal,dx10);
            ty               = _mm_mul_pd(fscal,dy10);
            tz               = _mm_mul_pd(fscal,dz10);

            /* Update vectorial force */
            fix1             = _mm_add_pd(fix1,tx);
            fiy1             = _mm_add_pd(fiy1,ty);
            fiz1             = _mm_add_pd(fiz1,tz);

            fjx0             = _mm_add_pd(fjx0,tx);
            fjy0             = _mm_add_pd(fjy0,ty);
            fjz0             = _mm_add_pd(fjz0,tz);

            /**************************
             * CALCULATE INTERACTIONS *
             **************************/

            /* Compute parameters for interactions between i and j atoms */
            qq20             = _mm_mul_pd(iq2,jq0);

            /* COULOMB ELECTROSTATICS */
            velec            = _mm_mul_pd(qq20,rinv20);
            felec            = _mm_mul_pd(velec,rinvsq20);

            /* Update potential sum for this i atom from the interaction with this j atom. */
            velecsum         = _mm_add_pd(velecsum,velec);

            fscal            = felec;

            /* Calculate temporary vectorial force */
            tx               = _mm_mul_pd(fscal,dx20);
            ty               = _mm_mul_pd(fscal,dy20);
            tz               = _mm_mul_pd(fscal,dz20);

            /* Update vectorial force */
            fix2             = _mm_add_pd(fix2,tx);
            fiy2             = _mm_add_pd(fiy2,ty);
            fiz2             = _mm_add_pd(fiz2,tz);

            fjx0             = _mm_add_pd(fjx0,tx);
            fjy0             = _mm_add_pd(fjy0,ty);
            fjz0             = _mm_add_pd(fjz0,tz);

            /**************************
             * CALCULATE INTERACTIONS *
             **************************/

            /* Compute parameters for interactions between i and j atoms */
            qq30             = _mm_mul_pd(iq3,jq0);

            /* COULOMB ELECTROSTATICS */
            velec            = _mm_mul_pd(qq30,rinv30);
            felec            = _mm_mul_pd(velec,rinvsq30);

            /* Update potential sum for this i atom from the interaction with this j atom. */
            velecsum         = _mm_add_pd(velecsum,velec);

            fscal            = felec;

            /* Calculate temporary vectorial force */
            tx               = _mm_mul_pd(fscal,dx30);
            ty               = _mm_mul_pd(fscal,dy30);
            tz               = _mm_mul_pd(fscal,dz30);

            /* Update vectorial force */
            fix3             = _mm_add_pd(fix3,tx);
            fiy3             = _mm_add_pd(fiy3,ty);
            fiz3             = _mm_add_pd(fiz3,tz);

            fjx0             = _mm_add_pd(fjx0,tx);
            fjy0             = _mm_add_pd(fjy0,ty);
            fjz0             = _mm_add_pd(fjz0,tz);

            gmx_mm_decrement_1rvec_2ptr_swizzle_pd(f+j_coord_offsetA,f+j_coord_offsetB,fjx0,fjy0,fjz0);

            /* Inner loop uses 87 flops */
        }

        if(jidx<j_index_end)
        {

            jnrA             = jjnr[jidx];
            j_coord_offsetA  = DIM*jnrA;

            /* load j atom coordinates */
            gmx_mm_load_1rvec_1ptr_swizzle_pd(x+j_coord_offsetA,
                                              &jx0,&jy0,&jz0);

            /* Calculate displacement vector */
            dx10             = _mm_sub_pd(ix1,jx0);
            dy10             = _mm_sub_pd(iy1,jy0);
            dz10             = _mm_sub_pd(iz1,jz0);
            dx20             = _mm_sub_pd(ix2,jx0);
            dy20             = _mm_sub_pd(iy2,jy0);
            dz20             = _mm_sub_pd(iz2,jz0);
            dx30             = _mm_sub_pd(ix3,jx0);
            dy30             = _mm_sub_pd(iy3,jy0);
            dz30             = _mm_sub_pd(iz3,jz0);

            /* Calculate squared distance and things based on it */
            rsq10            = gmx_mm_calc_rsq_pd(dx10,dy10,dz10);
            rsq20            = gmx_mm_calc_rsq_pd(dx20,dy20,dz20);
            rsq30            = gmx_mm_calc_rsq_pd(dx30,dy30,dz30);

            rinv10           = sse2_invsqrt_d(rsq10);
            rinv20           = sse2_invsqrt_d(rsq20);
            rinv30           = sse2_invsqrt_d(rsq30);

            rinvsq10         = _mm_mul_pd(rinv10,rinv10);
            rinvsq20         = _mm_mul_pd(rinv20,rinv20);
            rinvsq30         = _mm_mul_pd(rinv30,rinv30);

            /* Load parameters for j particles */
            jq0              = _mm_load_sd(charge+jnrA+0);

            fjx0             = _mm_setzero_pd();
            fjy0             = _mm_setzero_pd();
            fjz0             = _mm_setzero_pd();

            /**************************
             * CALCULATE INTERACTIONS *
             **************************/

            /* Compute parameters for interactions between i and j atoms */
            qq10             = _mm_mul_pd(iq1,jq0);

            /* COULOMB ELECTROSTATICS */
            velec            = _mm_mul_pd(qq10,rinv10);
            felec            = _mm_mul_pd(velec,rinvsq10);

            /* Update potential sum for this i atom from the interaction with this j atom. */
            velec            = _mm_unpacklo_pd(velec,_mm_setzero_pd());
            velecsum         = _mm_add_pd(velecsum,velec);

            fscal            = felec;

            fscal            = _mm_unpacklo_pd(fscal,_mm_setzero_pd());

            /* Calculate temporary vectorial force */
            tx               = _mm_mul_pd(fscal,dx10);
            ty               = _mm_mul_pd(fscal,dy10);
            tz               = _mm_mul_pd(fscal,dz10);

            /* Update vectorial force */
            fix1             = _mm_add_pd(fix1,tx);
            fiy1             = _mm_add_pd(fiy1,ty);
            fiz1             = _mm_add_pd(fiz1,tz);

            fjx0             = _mm_add_pd(fjx0,tx);
            fjy0             = _mm_add_pd(fjy0,ty);
            fjz0             = _mm_add_pd(fjz0,tz);

            /**************************
             * CALCULATE INTERACTIONS *
             **************************/

            /* Compute parameters for interactions between i and j atoms */
            qq20             = _mm_mul_pd(iq2,jq0);

            /* COULOMB ELECTROSTATICS */
            velec            = _mm_mul_pd(qq20,rinv20);
            felec            = _mm_mul_pd(velec,rinvsq20);

            /* Update potential sum for this i atom from the interaction with this j atom. */
            velec            = _mm_unpacklo_pd(velec,_mm_setzero_pd());
            velecsum         = _mm_add_pd(velecsum,velec);

            fscal            = felec;

            fscal            = _mm_unpacklo_pd(fscal,_mm_setzero_pd());

            /* Calculate temporary vectorial force */
            tx               = _mm_mul_pd(fscal,dx20);
            ty               = _mm_mul_pd(fscal,dy20);
            tz               = _mm_mul_pd(fscal,dz20);

            /* Update vectorial force */
            fix2             = _mm_add_pd(fix2,tx);
            fiy2             = _mm_add_pd(fiy2,ty);
            fiz2             = _mm_add_pd(fiz2,tz);

            fjx0             = _mm_add_pd(fjx0,tx);
            fjy0             = _mm_add_pd(fjy0,ty);
            fjz0             = _mm_add_pd(fjz0,tz);

            /**************************
             * CALCULATE INTERACTIONS *
             **************************/

            /* Compute parameters for interactions between i and j atoms */
            qq30             = _mm_mul_pd(iq3,jq0);

            /* COULOMB ELECTROSTATICS */
            velec            = _mm_mul_pd(qq30,rinv30);
            felec            = _mm_mul_pd(velec,rinvsq30);

            /* Update potential sum for this i atom from the interaction with this j atom. */
            velec            = _mm_unpacklo_pd(velec,_mm_setzero_pd());
            velecsum         = _mm_add_pd(velecsum,velec);

            fscal            = felec;

            fscal            = _mm_unpacklo_pd(fscal,_mm_setzero_pd());

            /* Calculate temporary vectorial force */
            tx               = _mm_mul_pd(fscal,dx30);
            ty               = _mm_mul_pd(fscal,dy30);
            tz               = _mm_mul_pd(fscal,dz30);

            /* Update vectorial force */
            fix3             = _mm_add_pd(fix3,tx);
            fiy3             = _mm_add_pd(fiy3,ty);
            fiz3             = _mm_add_pd(fiz3,tz);

            fjx0             = _mm_add_pd(fjx0,tx);
            fjy0             = _mm_add_pd(fjy0,ty);
            fjz0             = _mm_add_pd(fjz0,tz);

            gmx_mm_decrement_1rvec_1ptr_swizzle_pd(f+j_coord_offsetA,fjx0,fjy0,fjz0);

            /* Inner loop uses 87 flops */
        }

        /* End of innermost loop */

        gmx_mm_update_iforce_3atom_swizzle_pd(fix1,fiy1,fiz1,fix2,fiy2,fiz2,fix3,fiy3,fiz3,
                                              f+i_coord_offset+DIM,fshift+i_shift_offset);

        ggid                        = gid[iidx];
        /* Update potential energies */
        gmx_mm_update_1pot_pd(velecsum,kernel_data->energygrp_elec+ggid);

        /* Increment number of inner iterations */
        inneriter                  += j_index_end - j_index_start;

        /* Outer loop uses 19 flops */
    }

    /* Increment number of outer iterations */
    outeriter        += nri;

    /* Update outer/inner flops */

    inc_nrnb(nrnb,eNR_NBKERNEL_ELEC_W4_VF,outeriter*19 + inneriter*87);
}
/*
 * Gromacs nonbonded kernel:   nb_kernel_ElecCoul_VdwNone_GeomW4P1_F_sse2_double
 * Electrostatics interaction: Coulomb
 * VdW interaction:            None
 * Geometry:                   Water4-Particle
 * Calculate force/pot:        Force
 */
void
nb_kernel_ElecCoul_VdwNone_GeomW4P1_F_sse2_double
                    (t_nblist                    * gmx_restrict       nlist,
                     rvec                        * gmx_restrict          xx,
                     rvec                        * gmx_restrict          ff,
                     struct t_forcerec           * gmx_restrict          fr,
                     t_mdatoms                   * gmx_restrict     mdatoms,
                     nb_kernel_data_t gmx_unused * gmx_restrict kernel_data,
                     t_nrnb                      * gmx_restrict        nrnb)
{
    /* Suffixes 0,1,2,3 refer to particle indices for waters in the inner or outer loop, or
     * just 0 for non-waters.
     * Suffixes A,B refer to j loop unrolling done with SSE double precision, e.g. for the two different
     * jnr indices corresponding to data put in the four positions in the SIMD register.
     */
    int              i_shift_offset,i_coord_offset,outeriter,inneriter;
    int              j_index_start,j_index_end,jidx,nri,inr,ggid,iidx;
    int              jnrA,jnrB;
    int              j_coord_offsetA,j_coord_offsetB;
    int              *iinr,*jindex,*jjnr,*shiftidx,*gid;
    real             rcutoff_scalar;
    real             *shiftvec,*fshift,*x,*f;
    __m128d          tx,ty,tz,fscal,rcutoff,rcutoff2,jidxall;
    int              vdwioffset1;
    __m128d          ix1,iy1,iz1,fix1,fiy1,fiz1,iq1,isai1;
    int              vdwioffset2;
    __m128d          ix2,iy2,iz2,fix2,fiy2,fiz2,iq2,isai2;
    int              vdwioffset3;
    __m128d          ix3,iy3,iz3,fix3,fiy3,fiz3,iq3,isai3;
    int              vdwjidx0A,vdwjidx0B;
    __m128d          jx0,jy0,jz0,fjx0,fjy0,fjz0,jq0,isaj0;
    __m128d          dx10,dy10,dz10,rsq10,rinv10,rinvsq10,r10,qq10,c6_10,c12_10;
    __m128d          dx20,dy20,dz20,rsq20,rinv20,rinvsq20,r20,qq20,c6_20,c12_20;
    __m128d          dx30,dy30,dz30,rsq30,rinv30,rinvsq30,r30,qq30,c6_30,c12_30;
    __m128d          velec,felec,velecsum,facel,crf,krf,krf2;
    real             *charge;
    __m128d          dummy_mask,cutoff_mask;
    __m128d          signbit   = gmx_mm_castsi128_pd( _mm_set_epi32(0x80000000,0x00000000,0x80000000,0x00000000) );
    __m128d          one     = _mm_set1_pd(1.0);
    __m128d          two     = _mm_set1_pd(2.0);
    x                = xx[0];
    f                = ff[0];

    nri              = nlist->nri;
    iinr             = nlist->iinr;
    jindex           = nlist->jindex;
    jjnr             = nlist->jjnr;
    shiftidx         = nlist->shift;
    gid              = nlist->gid;
    shiftvec         = fr->shift_vec[0];
    fshift           = fr->fshift[0];
    facel            = _mm_set1_pd(fr->ic->epsfac);
    charge           = mdatoms->chargeA;

    /* Setup water-specific parameters */
    inr              = nlist->iinr[0];
    iq1              = _mm_mul_pd(facel,_mm_set1_pd(charge[inr+1]));
    iq2              = _mm_mul_pd(facel,_mm_set1_pd(charge[inr+2]));
    iq3              = _mm_mul_pd(facel,_mm_set1_pd(charge[inr+3]));

    /* Avoid stupid compiler warnings */
    jnrA = jnrB = 0;
    j_coord_offsetA = 0;
    j_coord_offsetB = 0;

    outeriter        = 0;
    inneriter        = 0;

    /* Start outer loop over neighborlists */
    for(iidx=0; iidx<nri; iidx++)
    {
        /* Load shift vector for this list */
        i_shift_offset   = DIM*shiftidx[iidx];

        /* Load limits for loop over neighbors */
        j_index_start    = jindex[iidx];
        j_index_end      = jindex[iidx+1];

        /* Get outer coordinate index */
        inr              = iinr[iidx];
        i_coord_offset   = DIM*inr;

        /* Load i particle coords and add shift vector */
        gmx_mm_load_shift_and_3rvec_broadcast_pd(shiftvec+i_shift_offset,x+i_coord_offset+DIM,
                                                 &ix1,&iy1,&iz1,&ix2,&iy2,&iz2,&ix3,&iy3,&iz3);

        fix1             = _mm_setzero_pd();
        fiy1             = _mm_setzero_pd();
        fiz1             = _mm_setzero_pd();
        fix2             = _mm_setzero_pd();
        fiy2             = _mm_setzero_pd();
        fiz2             = _mm_setzero_pd();
        fix3             = _mm_setzero_pd();
        fiy3             = _mm_setzero_pd();
        fiz3             = _mm_setzero_pd();

        /* Start inner kernel loop */
        for(jidx=j_index_start; jidx<j_index_end-1; jidx+=2)
        {

            /* Get j neighbor index, and coordinate index */
            jnrA             = jjnr[jidx];
            jnrB             = jjnr[jidx+1];
            j_coord_offsetA  = DIM*jnrA;
            j_coord_offsetB  = DIM*jnrB;

            /* load j atom coordinates */
            gmx_mm_load_1rvec_2ptr_swizzle_pd(x+j_coord_offsetA,x+j_coord_offsetB,
                                              &jx0,&jy0,&jz0);

            /* Calculate displacement vector */
            dx10             = _mm_sub_pd(ix1,jx0);
            dy10             = _mm_sub_pd(iy1,jy0);
            dz10             = _mm_sub_pd(iz1,jz0);
            dx20             = _mm_sub_pd(ix2,jx0);
            dy20             = _mm_sub_pd(iy2,jy0);
            dz20             = _mm_sub_pd(iz2,jz0);
            dx30             = _mm_sub_pd(ix3,jx0);
            dy30             = _mm_sub_pd(iy3,jy0);
            dz30             = _mm_sub_pd(iz3,jz0);

            /* Calculate squared distance and things based on it */
            rsq10            = gmx_mm_calc_rsq_pd(dx10,dy10,dz10);
            rsq20            = gmx_mm_calc_rsq_pd(dx20,dy20,dz20);
            rsq30            = gmx_mm_calc_rsq_pd(dx30,dy30,dz30);

            rinv10           = sse2_invsqrt_d(rsq10);
            rinv20           = sse2_invsqrt_d(rsq20);
            rinv30           = sse2_invsqrt_d(rsq30);

            rinvsq10         = _mm_mul_pd(rinv10,rinv10);
            rinvsq20         = _mm_mul_pd(rinv20,rinv20);
            rinvsq30         = _mm_mul_pd(rinv30,rinv30);

            /* Load parameters for j particles */
            jq0              = gmx_mm_load_2real_swizzle_pd(charge+jnrA+0,charge+jnrB+0);

            fjx0             = _mm_setzero_pd();
            fjy0             = _mm_setzero_pd();
            fjz0             = _mm_setzero_pd();

            /**************************
             * CALCULATE INTERACTIONS *
             **************************/

            /* Compute parameters for interactions between i and j atoms */
            qq10             = _mm_mul_pd(iq1,jq0);

            /* COULOMB ELECTROSTATICS */
            velec            = _mm_mul_pd(qq10,rinv10);
            felec            = _mm_mul_pd(velec,rinvsq10);

            fscal            = felec;

            /* Calculate temporary vectorial force */
            tx               = _mm_mul_pd(fscal,dx10);
            ty               = _mm_mul_pd(fscal,dy10);
            tz               = _mm_mul_pd(fscal,dz10);

            /* Update vectorial force */
            fix1             = _mm_add_pd(fix1,tx);
            fiy1             = _mm_add_pd(fiy1,ty);
            fiz1             = _mm_add_pd(fiz1,tz);

            fjx0             = _mm_add_pd(fjx0,tx);
            fjy0             = _mm_add_pd(fjy0,ty);
            fjz0             = _mm_add_pd(fjz0,tz);

            /**************************
             * CALCULATE INTERACTIONS *
             **************************/

            /* Compute parameters for interactions between i and j atoms */
            qq20             = _mm_mul_pd(iq2,jq0);

            /* COULOMB ELECTROSTATICS */
            velec            = _mm_mul_pd(qq20,rinv20);
            felec            = _mm_mul_pd(velec,rinvsq20);

            fscal            = felec;

            /* Calculate temporary vectorial force */
            tx               = _mm_mul_pd(fscal,dx20);
            ty               = _mm_mul_pd(fscal,dy20);
            tz               = _mm_mul_pd(fscal,dz20);

            /* Update vectorial force */
            fix2             = _mm_add_pd(fix2,tx);
            fiy2             = _mm_add_pd(fiy2,ty);
            fiz2             = _mm_add_pd(fiz2,tz);

            fjx0             = _mm_add_pd(fjx0,tx);
            fjy0             = _mm_add_pd(fjy0,ty);
            fjz0             = _mm_add_pd(fjz0,tz);

            /**************************
             * CALCULATE INTERACTIONS *
             **************************/

            /* Compute parameters for interactions between i and j atoms */
            qq30             = _mm_mul_pd(iq3,jq0);

            /* COULOMB ELECTROSTATICS */
            velec            = _mm_mul_pd(qq30,rinv30);
            felec            = _mm_mul_pd(velec,rinvsq30);

            fscal            = felec;

            /* Calculate temporary vectorial force */
            tx               = _mm_mul_pd(fscal,dx30);
            ty               = _mm_mul_pd(fscal,dy30);
            tz               = _mm_mul_pd(fscal,dz30);

            /* Update vectorial force */
            fix3             = _mm_add_pd(fix3,tx);
            fiy3             = _mm_add_pd(fiy3,ty);
            fiz3             = _mm_add_pd(fiz3,tz);

            fjx0             = _mm_add_pd(fjx0,tx);
            fjy0             = _mm_add_pd(fjy0,ty);
            fjz0             = _mm_add_pd(fjz0,tz);

            gmx_mm_decrement_1rvec_2ptr_swizzle_pd(f+j_coord_offsetA,f+j_coord_offsetB,fjx0,fjy0,fjz0);

            /* Inner loop uses 84 flops */
        }

        if(jidx<j_index_end)
        {

            jnrA             = jjnr[jidx];
            j_coord_offsetA  = DIM*jnrA;

            /* load j atom coordinates */
            gmx_mm_load_1rvec_1ptr_swizzle_pd(x+j_coord_offsetA,
                                              &jx0,&jy0,&jz0);

            /* Calculate displacement vector */
            dx10             = _mm_sub_pd(ix1,jx0);
            dy10             = _mm_sub_pd(iy1,jy0);
            dz10             = _mm_sub_pd(iz1,jz0);
            dx20             = _mm_sub_pd(ix2,jx0);
            dy20             = _mm_sub_pd(iy2,jy0);
            dz20             = _mm_sub_pd(iz2,jz0);
            dx30             = _mm_sub_pd(ix3,jx0);
            dy30             = _mm_sub_pd(iy3,jy0);
            dz30             = _mm_sub_pd(iz3,jz0);

            /* Calculate squared distance and things based on it */
            rsq10            = gmx_mm_calc_rsq_pd(dx10,dy10,dz10);
            rsq20            = gmx_mm_calc_rsq_pd(dx20,dy20,dz20);
            rsq30            = gmx_mm_calc_rsq_pd(dx30,dy30,dz30);

            rinv10           = sse2_invsqrt_d(rsq10);
            rinv20           = sse2_invsqrt_d(rsq20);
            rinv30           = sse2_invsqrt_d(rsq30);

            rinvsq10         = _mm_mul_pd(rinv10,rinv10);
            rinvsq20         = _mm_mul_pd(rinv20,rinv20);
            rinvsq30         = _mm_mul_pd(rinv30,rinv30);

            /* Load parameters for j particles */
            jq0              = _mm_load_sd(charge+jnrA+0);

            fjx0             = _mm_setzero_pd();
            fjy0             = _mm_setzero_pd();
            fjz0             = _mm_setzero_pd();

            /**************************
             * CALCULATE INTERACTIONS *
             **************************/

            /* Compute parameters for interactions between i and j atoms */
            qq10             = _mm_mul_pd(iq1,jq0);

            /* COULOMB ELECTROSTATICS */
            velec            = _mm_mul_pd(qq10,rinv10);
            felec            = _mm_mul_pd(velec,rinvsq10);

            fscal            = felec;

            fscal            = _mm_unpacklo_pd(fscal,_mm_setzero_pd());

            /* Calculate temporary vectorial force */
            tx               = _mm_mul_pd(fscal,dx10);
            ty               = _mm_mul_pd(fscal,dy10);
            tz               = _mm_mul_pd(fscal,dz10);

            /* Update vectorial force */
            fix1             = _mm_add_pd(fix1,tx);
            fiy1             = _mm_add_pd(fiy1,ty);
            fiz1             = _mm_add_pd(fiz1,tz);

            fjx0             = _mm_add_pd(fjx0,tx);
            fjy0             = _mm_add_pd(fjy0,ty);
            fjz0             = _mm_add_pd(fjz0,tz);

            /**************************
             * CALCULATE INTERACTIONS *
             **************************/

            /* Compute parameters for interactions between i and j atoms */
            qq20             = _mm_mul_pd(iq2,jq0);

            /* COULOMB ELECTROSTATICS */
            velec            = _mm_mul_pd(qq20,rinv20);
            felec            = _mm_mul_pd(velec,rinvsq20);

            fscal            = felec;

            fscal            = _mm_unpacklo_pd(fscal,_mm_setzero_pd());

            /* Calculate temporary vectorial force */
            tx               = _mm_mul_pd(fscal,dx20);
            ty               = _mm_mul_pd(fscal,dy20);
            tz               = _mm_mul_pd(fscal,dz20);

            /* Update vectorial force */
            fix2             = _mm_add_pd(fix2,tx);
            fiy2             = _mm_add_pd(fiy2,ty);
            fiz2             = _mm_add_pd(fiz2,tz);

            fjx0             = _mm_add_pd(fjx0,tx);
            fjy0             = _mm_add_pd(fjy0,ty);
            fjz0             = _mm_add_pd(fjz0,tz);

            /**************************
             * CALCULATE INTERACTIONS *
             **************************/

            /* Compute parameters for interactions between i and j atoms */
            qq30             = _mm_mul_pd(iq3,jq0);

            /* COULOMB ELECTROSTATICS */
            velec            = _mm_mul_pd(qq30,rinv30);
            felec            = _mm_mul_pd(velec,rinvsq30);

            fscal            = felec;

            fscal            = _mm_unpacklo_pd(fscal,_mm_setzero_pd());

            /* Calculate temporary vectorial force */
            tx               = _mm_mul_pd(fscal,dx30);
            ty               = _mm_mul_pd(fscal,dy30);
            tz               = _mm_mul_pd(fscal,dz30);

            /* Update vectorial force */
            fix3             = _mm_add_pd(fix3,tx);
            fiy3             = _mm_add_pd(fiy3,ty);
            fiz3             = _mm_add_pd(fiz3,tz);

            fjx0             = _mm_add_pd(fjx0,tx);
            fjy0             = _mm_add_pd(fjy0,ty);
            fjz0             = _mm_add_pd(fjz0,tz);

            gmx_mm_decrement_1rvec_1ptr_swizzle_pd(f+j_coord_offsetA,fjx0,fjy0,fjz0);

            /* Inner loop uses 84 flops */
        }

        /* End of innermost loop */

        gmx_mm_update_iforce_3atom_swizzle_pd(fix1,fiy1,fiz1,fix2,fiy2,fiz2,fix3,fiy3,fiz3,
                                              f+i_coord_offset+DIM,fshift+i_shift_offset);

        /* Increment number of inner iterations */
        inneriter                  += j_index_end - j_index_start;

        /* Outer loop uses 18 flops */
    }

    /* Increment number of outer iterations */
    outeriter        += nri;

    /* Update outer/inner flops */

    inc_nrnb(nrnb,eNR_NBKERNEL_ELEC_W4_F,outeriter*18 + inneriter*84);
}
