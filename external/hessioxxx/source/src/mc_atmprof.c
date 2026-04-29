/* ============================================================================

   Copyright (C) 2019, 2020  Konrad Bernloehr

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

/** @file mc_atmprof.c
 *  @short Interface to the atmospheric profile structure.
 *  @author Konrad Bernloehr
 *  @date   2019, 2020
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mc_atmprof.h"

#ifdef __GLIBC__
# ifndef _POSIX_C_SOURCE
/* Workaround for gcc in c99 language standard mode (too strict
   but glibc providing the function we need here). */
char *strdup(const char *s);
# endif
#endif

static AtmProf common_atmprof; /**< Keep track of atmospheric profiles loaded from text tables */


/* ------------------- get_common_atmprof -------------- */
/** 
 *  Make this copy of the atmospheric profile available elsewhere.
 */

AtmProf *get_common_atmprof(void)
{
   return &common_atmprof;
}

/* ---------------- set_common_atmprof --------------- */
/**
 *  Set the common profile from a separate copy.
 */

void set_common_atmprof(AtmProf *aprof)
{
   int i;

   if ( aprof == &common_atmprof )
      return; /* Don't try to copy it into itself. */

   /* Copy tabulated profiles to the extent they are known. */
   common_atmprof.atmprof_id = aprof->atmprof_id;
   if ( common_atmprof.atmprof_fname != NULL )
      free(common_atmprof.atmprof_fname);
   common_atmprof.atmprof_fname = NULL;
   if ( aprof->atmprof_fname != NULL )
      common_atmprof.atmprof_fname = strdup(aprof->atmprof_fname);
   common_atmprof.obslev = aprof->obslev;

   if ( common_atmprof.alt_km != NULL )
      free(common_atmprof.alt_km);
   common_atmprof.alt_km = NULL;
   if ( aprof->alt_km != NULL )
   {
      common_atmprof.alt_km = (double *) malloc(aprof->n_alt*sizeof(double));
      memcpy(common_atmprof.alt_km,aprof->alt_km,aprof->n_alt*sizeof(double));
   }
   if ( common_atmprof.rho != NULL )
      free(common_atmprof.rho);
   common_atmprof.rho = NULL;
   if ( aprof->rho != NULL )
   {
      common_atmprof.rho = (double *) malloc(aprof->n_alt*sizeof(double));
      memcpy(common_atmprof.rho,aprof->rho,aprof->n_alt*sizeof(double));
   }
   if ( common_atmprof.thick != NULL )
      free(common_atmprof.thick);
   common_atmprof.thick = NULL;
   if ( aprof->thick != NULL )
   {
      common_atmprof.thick = (double *) malloc(aprof->n_alt*sizeof(double));
      memcpy(common_atmprof.thick,aprof->thick,aprof->n_alt*sizeof(double));
   }
   if ( common_atmprof.refidx_m1 != NULL )
      free(common_atmprof.refidx_m1);
   common_atmprof.refidx_m1 = NULL;
   if ( aprof->refidx_m1 != NULL )
   {
      common_atmprof.refidx_m1 = (double *) malloc(aprof->n_alt*sizeof(double));
      memcpy(common_atmprof.refidx_m1,aprof->refidx_m1,aprof->n_alt*sizeof(double));
   }
   common_atmprof.n_alt = aprof->n_alt;

   /* Also copy the 5-layer CORSIKA EGS internal representation. */
   common_atmprof.have_lay5_param = aprof->have_lay5_param;
   for ( i=0; i<5; i++ )
   {
      common_atmprof.hlay[i] = aprof->hlay[i];
      common_atmprof.aatm[i] = aprof->aatm[i];
      common_atmprof.batm[i] = aprof->batm[i];
      common_atmprof.catm[i] = aprof->catm[i];
      common_atmprof.datm[i] = aprof->datm[i];
      common_atmprof.htoa = common_atmprof.hlay[5] = aprof->htoa;
   }
   /* Thickness at layer boundaries gets re-evaluated rather than copied. */
   for ( i=0; i<5; i++ )
      common_atmprof.thickl[i] = thickc(&common_atmprof.hlay[i]);
   common_atmprof.thickl[5] = 0.;
}

/* ---------------------- show_atmprof ---------------- */
/**
 *  Show a readable version of the tabulated atmospheric profile
 *  (basically like in the original tables, except for comments
 *  and extra unused columns), plus the 5-layer parametrization,
 *  if available.
 */

void show_atmprof(AtmProf *aprof)
{
   unsigned i;

   if ( aprof == NULL )
      return;

   if ( aprof->atmprof_fname != NULL )
      printf("\nAtmospheric profile %d (originally from file %s):\n", 
         aprof->atmprof_id, aprof->atmprof_fname);
   else
      printf("\nAtmospheric profile %d:\n", 
         aprof->atmprof_id);
   printf("Observation level is at %6.0f cm (%5.3f km) a.s.l.\n\n", 
      aprof->obslev, aprof->obslev*1e-5);
   printf(" Col. #1          #2           #3            #4\n");
   printf(" Altitude [km]  rho[g/cm^3]  thick[g/cm^2]  n-1\n");
   for ( i=0; i<aprof->n_alt; i++ )
   {
      printf("%10.3f    %11.5E  %11.5E  %11.5E\n", 
         aprof->alt_km[i], aprof->rho[i], aprof->thick[i], aprof->refidx_m1[i]);
   }
   if ( aprof->have_lay5_param )
   {
      printf("\nTop of atmosphere at %8.0f cm (%5.3f km).\n", aprof->htoa, aprof->htoa*1e-5);
      printf("Values for CORSIKA 5-layer parametrization:\n");
      printf("      HLAY        AATM         BATM          CATM       DATM=1/CATM    THICKL\n");
      for ( i=0; i<5; i++ )
      {
         printf("   %8.0f  %12.5E  %12.5E  %12.5E  %12.5E  %12.5E\n", 
            aprof->hlay[i], aprof->aatm[i], aprof->batm[i], aprof->catm[i], 
            aprof->datm[i], aprof->thickl[i]);
      }
   }
   else
      printf("\nNo values included for CORSIKA 5-layer parametrization.\n");
}

/* ------------------------- atmegs_ ---------------------- */
/**
 *  Fill the 5-layer parameters into the common atmospheric
 *  profile structure for keeping track of that together with
 *  the tabular input.
 */

void atmegs_(int *nlay, double *hlay,  double *aatm, double *batm, double *catm, double *datm, double *htoa)
{
   int i, n;
   double top = 120e5;

   if ( nlay == NULL || hlay == NULL || aatm == NULL || batm == NULL || catm == NULL )
   {
      fprintf(stderr,"Invalid or missing parameters for atmegs_()\n");
      return;
   }
   if ( htoa != NULL )
      top = *htoa;

   n = *nlay;
   if ( n < 0 )
      n = -n;

   if ( n == 5 )
   {
      for ( i=0; i<n; i++ )
      {
         common_atmprof.hlay[i] = hlay[i];
         common_atmprof.aatm[i] = aatm[i];
         common_atmprof.batm[i] = batm[i];
         common_atmprof.catm[i] = catm[i];
         if ( datm != NULL )
            common_atmprof.datm[i] = datm[i];
         else if ( catm[i] == 0. )
            common_atmprof.datm[i] = 0.;
         else
            common_atmprof.datm[i] = 1./catm[i];
      }
      common_atmprof.have_lay5_param = 1;
      for ( i=0; i<n; i++ )
         common_atmprof.thickl[i] = thickc(&common_atmprof.hlay[i]);
      common_atmprof.thickl[5] = 0.;
      if ( htoa == NULL )
         top = common_atmprof.aatm[4]*common_atmprof.catm[4];
      common_atmprof.htoa = common_atmprof.hlay[5] = top;
   }
   else
      fprintf(stderr,"Invalid number of atmospheric layers for atmegs_()\n");
}

void atmegs_default(void)
{
   int hl5 = common_atmprof.have_lay5_param; /* Preserve this value */
   /* Set parameters like in CORSIKA hard-coded default to enable for subsequent fit. */
   double hlay[] = { -5779.5e2, 4e5, 10e5, 40e5, 100e5, 112.8292e5 };
   /* DATA AATM / -186.555306D0,-94.919D0,  0.61289D0,0.D0,.01128292D0 / */
   /* DATA BATM / 1222.6562D0,1144.9069D0,1305.5948D0,540.1778D0,1.D0  / */
   /* DATA CATM / 994186.38D0,878153.55D0,636143.04D0,772170.16D0,1.D9 / */
   double aatm[] = { -186.555306, -94.919, 0.61289, 0., .01128292 };
   double batm[] = { 1222.6562, 1144.9069, 1305.5948, 540.1778, 1.0 };
   double catm[] = { 994186.38, 878153.55, 636143.04, 772170.16, 1.e9 };
   int n5 = 5;
   atmegs_(&n5,hlay,aatm,batm,catm,NULL,hlay+5);
   /* Restore value */
   common_atmprof.have_lay5_param = hl5;
}

double rhofc (double *height)
{
   if ( (*height) < common_atmprof.hlay[1] )
      return common_atmprof.batm[0] * common_atmprof.datm[0] * exp(-(*height)*common_atmprof.datm[0]);
   else if ( (*height) < common_atmprof.hlay[2] )
      return common_atmprof.batm[1] * common_atmprof.datm[1] * exp(-(*height)*common_atmprof.datm[1]);
   else if ( (*height) < common_atmprof.hlay[3] )
      return common_atmprof.batm[2] * common_atmprof.datm[2] * exp(-(*height)*common_atmprof.datm[2]);
   else if ( (*height) < common_atmprof.hlay[4] )
      return common_atmprof.batm[3] * common_atmprof.datm[3] * exp(-(*height)*common_atmprof.datm[3]);
   else
      return common_atmprof.datm[4];
}

double thickc (double *height)
{
   if ( (*height) < common_atmprof.hlay[1] )
      return common_atmprof.aatm[0] + common_atmprof.batm[0] * exp(-(*height)*common_atmprof.datm[0]);
   else if ( (*height) < common_atmprof.hlay[2] )
      return common_atmprof.aatm[1] + common_atmprof.batm[1] * exp(-(*height)*common_atmprof.datm[1]);
   else if ( (*height) < common_atmprof.hlay[3] )
      return common_atmprof.aatm[2] + common_atmprof.batm[2] * exp(-(*height)*common_atmprof.datm[2]);
   else if ( (*height) < common_atmprof.hlay[4] )
      return common_atmprof.aatm[3] + common_atmprof.batm[3] * exp(-(*height)*common_atmprof.datm[3]);
   else
      return common_atmprof.aatm[4] - (*height) * common_atmprof.datm[4];
}

/* This value (inc contrast to etadsn as used in atmo.c) remains at the
   value as used in CORSIKA and does not get adapted to the actual ratio of
   refractivity over density. */
static double etadsn0 = 0.000283 * 994186.38 / 1222.656; /* CORSIKA default */;

double refidc (double *height)
{
   return 1.0 + etadsn0 * rhofc(height);
}

double refim1c (double *height)
{
   return etadsn0 * rhofc(height);
}

double heighc (double *thick)
{
   if ( (*thick) > common_atmprof.thickl[1] )
      return common_atmprof.catm[0] * log(common_atmprof.batm[0] / ((*thick) - common_atmprof.aatm[0]));
   else if ( (*thick) > common_atmprof.thickl[2] )
      return common_atmprof.catm[1] * log(common_atmprof.batm[1] / ((*thick) - common_atmprof.aatm[1]));
   else if ( (*thick) > common_atmprof.thickl[3] )
      return common_atmprof.catm[2] * log(common_atmprof.batm[2] / ((*thick) - common_atmprof.aatm[2]));
   else if ( (*thick) > common_atmprof.thickl[4] )
      return common_atmprof.catm[3] * log(common_atmprof.batm[3] / ((*thick) - common_atmprof.aatm[3]));
   else
      return (common_atmprof.aatm[4] - (*thick)) * common_atmprof.catm[4];
}

