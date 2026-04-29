/* ============================================================================

   Copyright (C) 2019  Konrad Bernloehr

   This file is part of the eventio/hessio library and other packages.

   The eventio/hessio library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with this library. If not, see <http://www.gnu.org/licenses/>.

============================================================================ */

/** @file mc_atmprof.h 
 *  @short A data structure shared between io_simtel.c and atmo.c -
 *  which is used by both sim_telarray and the CORSIKA IACT/atmo package.
 *  Filling the structure from text format tables is handled by atmo.c
 *  while EventIO input and output is handled by io_simtel.c.
 *  The purpose of the structure is for keeping track of the profile
 *  actually used. Evaluating/interpolating it is handled elsewhere.
 *  In addition to the tabulated profiles, it can also keep track of the
 *  5-layer parametrization as hard-wired into the CORSIKA EGS part.
 *
 *  @author  Konrad Bernloehr 
 *  @date    2019
 */
/* ========================================================= */

#ifndef HAVE_MC_ATMPROF

# define HAVE_MC_ATMPROF 1

#ifdef __cplusplus
extern "C" {
#endif

/** Atmospheric profile as stored in atmprof*.dat files - the actually used columns only. */

struct atmospheric_profile
{
   int atmprof_id;     /**< Profile ID number ('atmprof<i>.dat') or 99 */
   char *atmprof_fname;/**< Original name of atmospheric profile loaded */
   double obslev;      /**< Observation level [cm], a.s.l., as used in CORSIKA. */
   unsigned n_alt;     /**< Number of altitude levels */
   double *alt_km;     /**< Altitude a.s.l. [km] at each level */
   double *rho;        /**< Density [g/cm^3] at each level */
   double *thick;      /**< Vertical column density from space to given level [g/cm^2] */
   double *refidx_m1;  /**< Index of refraction minus one (n-1) at given level */
   int have_lay5_param;/**< Is 1 if the 5-layer CORSIKA built-in parametrization is known, 0 if not. */
   double hlay[6];     /**< Layer bounderies a.s.l. [cm]; see ATMLAY CORSIKA inputs card */
   double aatm[5];     /**< See ATMA CORSIKA inputs card */
   double batm[5];     /**< See ATMB CORSIKA inputs card */
   double catm[5];     /**< See ATMC CORSIKA inputs card */
   double datm[5];     /**< Inverse of catm values (if non-zero) */
   double thickl[6];   /**< Atmospheric thickness at given hlay heights */
   double htoa;        /**< Height (a.s.l.) at top of atmosphere [cm] */
};
typedef struct atmospheric_profile AtmProf;

/* Make the common profile available */
AtmProf *get_common_atmprof(void);
/* Set the common profile from a separate copy. */
void set_common_atmprof(AtmProf *atmprof);
void show_atmprof (AtmProf *atmprof);

void atmegs_(int *nlay, double *hlay,  double *aatm, 
   double *batm, double *catm, double *datm, double *htoa);
void atmegs_default(void);

/** C-called functions equivalent to the CORSIKA-built-in functions
   to evaluate the 5-layer parametrization.
   Assumes that these parameters have been set before.
   Where the numerical table is available it should be used once to
   initialize the atmospheric profile and then use the corresponding
   rhofx_(), ... functions for the evaulation instead.
 */

double rhofc (double *height);
double thickc (double *height);
double refidc (double *height);
double refim1c (double *height);
double heighc (double *thick);

#ifdef __cplusplus
}
#endif

#endif

